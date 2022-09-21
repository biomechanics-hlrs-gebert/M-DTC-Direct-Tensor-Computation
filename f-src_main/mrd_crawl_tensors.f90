!------------------------------------------------------------------------------
!> MeRaDat - Crawl Tensors from the Direct Tensor Computation
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    19.03.2022
!> LastMod: 20.03.2022
!
!> @brief:
!> Program to colect the results of the Direct Tensor Computation. 
!
!> @description:
!> This program can be included in struct_process.f90, which prevents
!> redundancies. 
!> However, the transparency, modularity etc. are better this way.
!------------------------------------------------------------------------------
PROGRAM mrd_crawl_tensors

USE global_std
USE puredat
USE meta
USE user_interaction
USE formatted_plain
USE mechanical

IMPLICIT NONE

CHARACTER(*), PARAMETER :: mrd_dbg_lvl = "PRODUCTION" ! "DEBUG"

TYPE(tLeaf), POINTER :: leaf_domain_no, leaf_density, leaf_tensors, leaf_pos
TYPE(tBranch)        :: rank_data
TYPE(tensor_2nd_rank_R66), DIMENSION(3) :: tensor

INTEGER(ik), DIMENSION(:,:), ALLOCATABLE :: Domain_stats
INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0

INTEGER(ik) :: alloc_stat, parts, fh_covo, fh_cr1, fh_cr2, fh_crdiag, fh_tens, &
    domains_crawled = 0_ik, nn_comm, par_domains, activity_size, size_mpi, &
    ii, jj, kk, ll, mm, tt, handle, &
    aun, domains_per_comm, rank_mpi, No_of_domains, local_domain_no, last_domain_rank

REAL(rk) :: local_domain_density, sym, start, end
REAL(rk), DIMENSION(3)   :: local_domain_opt_pos
REAL(rk), DIMENSION(6,6) :: local_domain_tensor

INTEGER(8), DIMENSION(:), ALLOCATABLE :: dat_domains
REAL   (8), DIMENSION(:), ALLOCATABLE :: dat_densities, dat_tensors, dat_pos

CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(mcl) :: cmd_arg_history='', target_path, file_head, binary, activity_file, stat=""
CHARACTER(20)  :: dmn_char
CHARACTER(7)   :: crdiag_suf = ".crdiag"
CHARACTER(5)   :: covo_suf = ".covo"
CHARACTER(4)   :: cr1_suf = ".cr1", cr2_suf = ".cr2"
CHARACTER(1)   :: tt_char, restart_cmd_arg='U' ! U = 'undefined'

LOGICAL :: success=.FALSE., stat_exists=.FALSE., opened=.FALSE., last_domain=.FALSE.
LOGICAL :: write_to_diag = .FALSE., abrt = .FALSE.

CALL CPU_TIME(start)

!------------------------------------------------------------------------------
! Parse the command arguments
!------------------------------------------------------------------------------
CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history)

IF (in%full=='') THEN
    CALL usage(binary, [ &
        "Give the Output Meta File of the tensors, computed by DTC.", &
        "The *.covo file printed acts as the input file for M-ROT. ", &
        "Only empty fields are computed in M-ROT.                  " ])    

    !------------------------------------------------------------------------------
    ! On std_out since file of std_out is not spawned
    !------------------------------------------------------------------------------
    CALL print_err_stop(6, "No input file given", 1)
END IF

!------------------------------------------------------------------------------
! Open the given meta file and parse its basename
!------------------------------------------------------------------------------
CALL meta_invoke(m_rry)

!------------------------------------------------------------------------------
! Redirect std_out into a file in case std_out is not useful by environment.
!------------------------------------------------------------------------------
std_out = 6 ! determine_stout()

!------------------------------------------------------------------------------
! Spawn standard out after(!) the basename is known
!------------------------------------------------------------------------------
IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "])

!------------------------------------------------------------------------------
! Parse input
!------------------------------------------------------------------------------
CALL meta_check_contains_program ('TENSOR-COMPUTATION', m_rry, success)

IF (.NOT. success) THEN
    CALL print_trimmed(std_out, &
        "The program 'TENSOR-COMPUTATION' apparently did not run successfully. &
        &Maybe it crashed or was stopped by purpose. However, it can also point &
        &to an incorrect implementation.", &
        FMT_WRN)    
    WRITE(std_out, FMT_WRN_SEP)
END IF

CALL meta_read('LO_BNDS_DMN_RANGE' , m_rry, xa_d, stat); IF(stat/="") STOP
CALL meta_read('UP_BNDS_DMN_RANGE' , m_rry, xe_d, stat); IF(stat/="") STOP
CALL meta_read('MESH_PER_SUB_DMN'  , m_rry, parts, stat); IF(stat/="") STOP
CALL meta_read('PROCESSORS'        , m_rry, size_mpi, stat); IF(stat/="") STOP

INQUIRE(UNIT=fhmei, OPENED=opened)
IF(opened) CLOSE (fhmei)

!------------------------------------------------------------------------------
! Open the file to store the tensors positioned in their global coordinates
!------------------------------------------------------------------------------
fh_covo   = give_new_unit(); CALL meta_start_ascii(fh_covo, covo_suf)
fh_cr1    = give_new_unit(); CALL meta_start_ascii(fh_cr1, cr1_suf)
fh_cr2    = give_new_unit(); CALL meta_start_ascii(fh_cr2, cr2_suf)
fh_crdiag = give_new_unit(); CALL meta_start_ascii(fh_crdiag, crdiag_suf)

CALL write_tensor_2nd_rank_R66_header(fh_covo)
CALL write_tensor_2nd_rank_R66_header(fh_cr1)
CALL write_tensor_2nd_rank_R66_header(fh_cr2)
CALL write_tensor_2nd_rank_R66_header(fh_crdiag)

!------------------------------------------------------------------------------
! Calculate the number of domains per communicator of the DTC
! 
! IMPORTANT: This calculations must be the exact same like in the 
! struct_process.f90 main program! Otherwise, the data may be garbage.
!------------------------------------------------------------------------------
No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

par_domains = (size_mpi-1) / parts

domains_per_comm = CEILING(REAL(No_of_domains, rk) / REAL(par_domains, rk), ik)

domains_per_comm = domains_per_comm * 2_ik

!------------------------------------------------------------------------------
! Allocate memory. Only necessary for one single PETSc sub comm of the DTC.
! Directories/Ranks get crawled one by one.
!------------------------------------------------------------------------------
ALLOCATE(Domain_stats(No_of_domains, 2_ik))

Domain_stats(:,:) = 0_ik

!------------------------------------------------------------------------------
! Read activity/status file of the DTC computation to get the list of domains
! computed.
!------------------------------------------------------------------------------
aun = give_new_unit()

activity_file = TRIM(in%p_n_bsnm)//".stat"
INQUIRE(FILE = TRIM(activity_file), EXIST=stat_exists, SIZE=activity_size)

IF(.NOT. stat_exists) THEN
    mssg = "No status file found. Please check the input and/or computation."
    CALL print_err_stop(std_out, mssg, 1_ik) 
END IF 

IF(activity_size /= No_of_domains * ik) THEN
    mssg = "File size and number of computed domains does not match. &
        &Please check the boundaries, hexdump the status file and have a look, &
        &if the kind of the integers of this file is correct (INT64/8)."
    CALL print_err_stop(std_out, mssg, 1_ik)
END IF

OPEN (UNIT=aun, FILE=TRIM(activity_file), ACCESS="STREAM")

READ (aun) Domain_stats(:, 1)

CLOSE(aun)

!------------------------------------------------------------------------------
! Crawl data
!------------------------------------------------------------------------------
DO rank_mpi = 1, size_mpi-1, parts
    nn_comm = 1_ik

    target_path = TRIM(in%p_n_bsnm)//"/"
    file_head = "results"

    pro_path = ""
    pro_name = ""

    WRITE(pro_path,'(A,A,I7.7,A)') TRIM(target_path), "Rank_", rank_mpi, "/"
    WRITE(pro_name,'(A,A,I7.7)')   TRIM(file_head), "_", rank_mpi

    !------------------------------------------------------------------------------
    ! Read PureDat tree
    !------------------------------------------------------------------------------
    rank_data = read_tree()

    CALL open_stream_files(rank_data, "read" , "old")

    !------------------------------------------------------------------------------
    ! Read domain numbers
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 2_pd_ik, leaf_domain_no, success)

    IF(ALLOCATED(dat_domains)) DEALLOCATE(dat_domains)
    ALLOCATE(dat_domains(leaf_domain_no%dat_no), stat=alloc_stat)
    
    CALL print_err_stop(std_out, "Allocating 'dat_domains' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, leaf_domain_no, dat_domains)

    last_domain_rank = dat_domains(leaf_domain_no%dat_no)

    !------------------------------------------------------------------------------
    ! Read effective density
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 20_pd_ik, leaf_density, success)
    
    IF(ALLOCATED(dat_densities)) DEALLOCATE(dat_densities)
    ALLOCATE(dat_densities(leaf_density%dat_no), stat=alloc_stat)
    
    CALL print_err_stop(std_out, "Allocating 'dat_densities' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, leaf_density, dat_densities)

    DO tt = 1, 3
        !------------------------------------------------------------------------------
        ! Reset steering variables
        !------------------------------------------------------------------------------
        Domain_stats(:, 2) = 0_ik
        last_domain = .FALSE.

        !------------------------------------------------------------------------------
        ! Read effective stiffness and euler (!) positions
        !------------------------------------------------------------------------------
        SELECT CASE(tt)
            CASE(1)
                CALL get_leaf_with_num(rank_data, 8_pd_ik, leaf_tensors, success)
                local_domain_opt_pos = 0._rk
                fh_tens     = fh_covo 
            CASE(2)
                CALL get_leaf_with_num(rank_data, 15_pd_ik, leaf_tensors, success)
                CALL get_leaf_with_num(rank_data, 13_pd_ik, leaf_pos    , success)
                fh_tens = fh_cr1
            CASE(3)
                CALL get_leaf_with_num(rank_data, 19_pd_ik, leaf_tensors, success)
                CALL get_leaf_with_num(rank_data, 17_pd_ik, leaf_pos    , success)
                fh_tens = fh_cr2
        END SELECT

        !------------------------------------------------------------------------------
        ! Required for distinction between regular domains and domains which appear
        ! twice in the dataset.
        !------------------------------------------------------------------------------
        handle = fh_tens

        IF(ALLOCATED(dat_tensors)) DEALLOCATE(dat_tensors)
        ALLOCATE(dat_tensors(leaf_tensors%dat_no), stat=alloc_stat)
    
        WRITE(tt_char, '(I0)') tt
        CALL print_err_stop(std_out, "Allocating 'dat_tensors("//tt_char//")' failed.", alloc_stat)
        CALL pd_read_leaf(rank_data%streams, leaf_tensors, dat_tensors)

        !------------------------------------------------------------------------------
        ! Rotation Vector CR_1
        !------------------------------------------------------------------------------
        IF(tt/= 1) THEN
            IF(ALLOCATED(dat_pos)) DEALLOCATE(dat_pos)
            ALLOCATE(dat_pos(leaf_pos%dat_no), stat=alloc_stat)
            CALL print_err_stop(std_out, "Allocating 'dat_pos("//tt_char//")' failed.", alloc_stat)
            CALL pd_read_leaf(rank_data%streams, leaf_pos, dat_pos)
        END IF 

        !------------------------------------------------------------------------------
        ! Begin searching for the domain to extract
        !------------------------------------------------------------------------------
        DO ii = 1, domains_Per_comm - 1_ik

            !------------------------------------------------------------------------------
            ! Check if it is the last domain. Exit loop if it is not, but the last
            ! iteration was.
            !------------------------------------------------------------------------------
            IF (dat_domains(ii) .EQ. last_domain_rank) THEN
                last_domain = .TRUE.
            ELSE 
                IF (last_domain) EXIT
            END IF 

            !------------------------------------------------------------------------------
            ! Compare the current domain with list of computed domains during DTC 
            ! If it is contained - the tensor can be read. If not, issue an error, since
            ! the domain given apparently is wrong
            !------------------------------------------------------------------------------
            DO jj = 1, No_of_domains
                !------------------------------------------------------------------------------
                ! Reset the information stored before proceeding to the next one.
                !------------------------------------------------------------------------------

                IF(Domain_stats(jj, 1) == dat_domains(ii)) THEN

                    !------------------------------------------------------------------------------
                    ! Check whether the domain was already found.
                    !------------------------------------------------------------------------------
                    IF((Domain_stats(jj, 2) == 0_ik) .OR. (last_domain)) THEN
                        Domain_stats(jj, 2) = 1_ik
                    ELSE
!                       WRITE(std_out, FMT_DBG_AxI0) "Domain_stats(:): ", Domain_stats(:,1)   
!                       WRITE(std_out, FMT_DBG_AxI0) "Domain_stats(:): ", Domain_stats(:,2)   

                        WRITE(std_out, FMT_DBG_xAL)  "last_domain: ", last_domain
                        WRITE(std_out, FMT_DBG_AxI0) "Domain_stats(jj, 1)", Domain_stats(jj, 1) 
                        WRITE(std_out, FMT_DBG_AxI0) "Domain_stats(jj, 2)", Domain_stats(jj, 2) 
                        WRITE(std_out, FMT_DBG_AxI0) "last_domain_rank", last_domain_rank

                        mssg = "The domain was already found. The second appearance of a domain &
                            &number in a 'valid' state is a hint for a corrupted data set. &
                            &One would expect to find domain numbers only once in all leafs of &
                            &domain numbers, scattered in the rank-directories, as long as the &
                            &first entries are monotonously increasing."
                        ! CALL print_err_stop(std_out, mssg, 1_ik)
                        CALL print_trimmed(std_out, TRIM(mssg), FMT_WRN)
                        
                        WRITE(std_out, FMT_WRN_SEP)

                        write_to_diag = .TRUE.
                    END IF

                    !------------------------------------------------------------------------------
                    ! Effective domain density, which is crawled at the moment
                    !------------------------------------------------------------------------------
                    local_domain_density = dat_densities(ii)

                    !------------------------------------------------------------------------------
                    ! Get the position of the optimized tensor
                    !------------------------------------------------------------------------------
                    IF(tt/= 1) local_domain_opt_pos = dat_pos(ii:ii + 2_ik)

                    !------------------------------------------------------------------------------
                    ! Local tensor
                    !------------------------------------------------------------------------------
                    local_domain_tensor = 0._rk
                    mm = (ii-1) * 36 + 1
                    DO kk = 1, 6
                    DO ll = 1, 6
                        local_domain_tensor(kk, ll) = dat_tensors(mm)
                        mm = mm + 1
                    END DO 
                    END DO

                    !------------------------------------------------------------------------------
                    ! Domain number, which is crawled at the moment
                    !------------------------------------------------------------------------------
                    local_domain_no = dat_domains(ii)

                    domains_crawled = domains_crawled + 1_ik

                    IF(mrd_dbg_lvl == "DEBUG") THEN
                        WRITE(dmn_char, '(I0)') local_domain_no
                        CALL write_matrix(std_out, "Effective stiffness of domain "//dmn_char, &
                            local_domain_tensor, fmti='std', unit='MPa')
                    END IF 

                    !------------------------------------------------------------------------------
                    ! Check symmetry of the tensor
                    !------------------------------------------------------------------------------
                    sym = check_sym (local_domain_tensor)

                    !------------------------------------------------------------------------------
                    ! Write the data to the file.
                    !------------------------------------------------------------------------------
                    tensor(tt)%dmn  = dat_domains(ii)
                    tensor(tt)%bvtv = dat_densities(ii)
                    tensor(tt)%sym  = sym
                    tensor(tt)%mat  = local_domain_tensor
                    tensor(tt)%pos  = local_domain_opt_pos

                    SELECT CASE (tt)
                        CASE(1); tensor(tt)%opt_crit = "covo" 
                        CASE(2); tensor(tt)%opt_crit = "cr1" 
                        CASE(3); tensor(tt)%opt_crit = "cr2" 
                    END SELECT

                    IF(write_to_diag) handle = fh_crdiag

                    CALL write_tensor_2nd_rank_R66_row(handle, tensor(tt))

                    !------------------------------------------------------------------------------
                    ! Reset diagnostic stuff
                    !------------------------------------------------------------------------------
                    write_to_diag = .FALSE.
                    handle = fh_tens

                END IF ! (Domain_stats(jj, 1) == leaf_domain_no%p_int8(ii)) THEN
            END DO

        END DO ! ii = 1, domains_Per_comm - 1_ik
    END DO ! tt = 1, 3
    !------------------------------------------------------------------------------
    ! Close the files of this particular directory (Rank_000xyz)
    ! basically a reset of the tree which was read before.
    !------------------------------------------------------------------------------
    CALL close_stream_files(rank_data)

    CALL nullify_pd_root()

    !------------------------------------------------------------------------------
    ! Increment the domain per PETSc sub-comm counter
    !------------------------------------------------------------------------------
    nn_comm = nn_comm + 1_ik 
END DO

CALL CPU_TIME(end)

!------------------------------------------------------------------------------
! Finish program
!------------------------------------------------------------------------------
CALL meta_stop_ascii(fh_covo, covo_suf)
CALL meta_stop_ascii(fh_cr1, cr1_suf)
CALL meta_stop_ascii(fh_cr2, cr2_suf)
CALL meta_stop_ascii(fh_crdiag, crdiag_suf)

WRITE(std_out, FMT_TXT_xAF0) "Program finished successfully in ", end-start, " seconds."
WRITE(std_out, FMT_TXT_SEP)

IF(std_out /= 6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END PROGRAM mrd_crawl_tensors
