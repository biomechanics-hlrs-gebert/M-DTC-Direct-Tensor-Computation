!------------------------------------------------------------------------------
!> MeRaDat - Crawl Tensors from the Direct Tensor Computation
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    19.03.2022
!> LastMod: 09.11.2022
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
USE ser_binary

IMPLICIT NONE

CHARACTER(*), PARAMETER :: mrd_dbg_lvl =  "PRODUCTION" ! "DEBUG"

TYPE(tLeaf), POINTER :: domain_no, eff_num_stiffness, density, tensors, leaf_pos, &
    no_elems, no_nodes, t_start, t_duration, collected_logs
TYPE(tBranch)        :: rank_data
TYPE(domain_data), DIMENSION(3) :: tensor

INTEGER(ik), DIMENSION(:,:), ALLOCATABLE :: Domain_status
INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0, vdim=0, grid=0

INTEGER(ik) :: alloc_stat, parts, fh_covo, fh_cr1, fh_cr2, fh_tens, dat_no_highscore, &
    domains_crawled = 0_ik, nn_comm, activity_size, size_mpi, fh_covo_num, dat_no, &
    ii, kk, ll, mm, tt, xx, par_dmn_number, jj, current_domain, pntr, no_lc, ma_el_order, &
    aun, domains_per_comm = 0, rank_mpi, No_of_domains, local_domain_no, last_domain_rank

REAL(rk) :: local_domain_density, sym, start, end
REAL(rk), DIMENSION(3)   :: local_domain_opt_pos, spcng, dmn_size
REAL(rk), DIMENSION(6,6) :: local_domain_tensor
REAL(rk), DIMENSION(:,:), ALLOCATABLE :: local_num_tensor

INTEGER(8), DIMENSION(:), ALLOCATABLE :: dat_domains, dat_no_elems, dat_no_nodes, dat_collected_logs
REAL(8), DIMENSION(:), ALLOCATABLE :: dat_densities, dat_eff_num_stiffnesses, dat_tensors, &
    dat_pos, dat_t_start, dat_t_duration

CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(mcl) :: cmd_arg_history='', target_path, file_head, binary, activity_file, stat=""
CHARACTER(scl) :: mi_el_type

CHARACTER(100*mcl) :: string

CHARACTER(20)  :: dmn_char
CHARACTER(9)   :: covo_num_suf = ".covo.num"
CHARACTER(5)   :: covo_suf = ".covo"
CHARACTER(4)   :: cr1_suf = ".cr1", cr2_suf = ".cr2"
CHARACTER(1)   :: tt_char, restart_cmd_arg='U' ! U = 'undefined'

LOGICAL :: success=.FALSE., stat_exists=.FALSE., opened=.FALSE., last_domain=.FALSE., ex=.FALSE.

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

    CALL print_err_stop(6, "No input file given", 1)
END IF

!------------------------------------------------------------------------------
! Open the given meta file and parse its basename
!------------------------------------------------------------------------------
CALL meta_invoke(m_rry)
! Normally, one would call meta_append which wraps meta_invoke and meta_continue.
! However, meta_continue provides the assignment out=in which is required to set 
! the proper basename for meta_start_ascii.
out = in

CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "])

!------------------------------------------------------------------------------
! Parse input
!------------------------------------------------------------------------------
CALL meta_check_contains_program ('TENSOR_COMPUTATION', m_rry, success)

IF (.NOT. success) THEN
    CALL print_trimmed(6, &
        "The program 'TENSOR_COMPUTATION' probably did not run successfully. &
        &Maybe it crashed or was stopped by purpose. However, it can also point &
        &to an incorrect implementation.", &
        FMT_WRN)    
    WRITE(6, FMT_WRN_SEP)
END IF

CALL meta_read('LO_BNDS_DMN_RANGE',m_rry,xa_d,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
CALL meta_read('UP_BNDS_DMN_RANGE',m_rry,xe_d,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)

!------------------------------------------------------------------------------
! Macro Element order via string for more flexibility
!------------------------------------------------------------------------------
CALL meta_read('MACRO_ELMNT_ORDER',m_rry,ma_el_order,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
CALL meta_read('MICRO_ELMNT_TYPE' ,m_rry,mi_el_type,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)

CALL meta_read('PROCESSORS' ,m_rry,size_mpi,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
CALL meta_read('SIZE_DOMAIN',m_rry,dmn_size,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
CALL meta_read('SPACING'    ,m_rry,spcng,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
CALL meta_read('DIMENSIONS' ,m_rry,vdim,stat);IF(stat/="") WRITE(6,FMT_ERR) TRIM(stat)
IF(stat/="") STOP

! Optional:
CALL meta_read('MESH_PER_SUB_DMN' ,m_rry,parts,stat)
IF(stat/="") THEN
    WRITE(6, FMT_TXT) "Keyword "//TRIM(stat)//" not active."
    WRITE(6, FMT_WRN_SEP)
END IF 

IF (ma_el_order == 1) THEN
    no_lc = 24
ELSE IF (ma_el_order == 2) THEN
    no_lc = 60
ELSE
    CALL print_err_stop(6, "Macro element order not recognized. Chosse '1' or '2'.", 1)
END IF 

ALLOCATE(local_num_tensor(no_lc, no_lc))

INQUIRE(UNIT=fhmei, OPENED=opened)
IF(opened) CLOSE (fhmei)

!------------------------------------------------------------------------------
! Open the file to store the tensors positioned in their global coordinates
!------------------------------------------------------------------------------
fh_covo     = give_new_unit(); CALL meta_start_ascii(fh_covo, covo_suf)
fh_cr1      = give_new_unit(); CALL meta_start_ascii(fh_cr1, cr1_suf)
fh_cr2      = give_new_unit(); CALL meta_start_ascii(fh_cr2, cr2_suf)
fh_covo_num = give_new_unit(); CALL meta_start_ascii(fh_covo_num, covo_num_suf)

string = write_tensor_2nd_rank_R66_header()
WRITE(fh_covo, '(A)') TRIM(string)
WRITE(fh_cr1, '(A)') TRIM(string)
WRITE(fh_cr2, '(A)') TRIM(string)

string = write_eff_num_stiff_header(no_lc)
WRITE(fh_covo_num, '(A)') TRIM(string)

!------------------------------------------------------------------------------
! Calculate the number of domains per communicator of the DTC
! 
! IMPORTANT: This calculations must be the exact same like in the 
! struct_process.f90 main program! Otherwise, the data may be garbage.
!------------------------------------------------------------------------------
No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

!------------------------------------------------------------------------------
! Allocate memory. Only necessary for one single PETSc sub comm of the DTC.
! Directories/Ranks get crawled one by one.
!------------------------------------------------------------------------------
ALLOCATE(Domain_status(No_of_domains,4))

Domain_status = 0_ik

!------------------------------------------------------------------------------
! Read activity/status file of the DTC computation to get the list of domains
! computed.
!------------------------------------------------------------------------------
aun = give_new_unit()

activity_file = TRIM(in%p_n_bsnm)//".status"
INQUIRE(FILE = TRIM(activity_file), EXIST=stat_exists, SIZE=activity_size)

IF(.NOT. stat_exists) THEN
    mssg = "No status file found. Please check the input and/or computation."
    CALL print_err_stop(6, mssg, 1_ik) 
END IF 

IF(activity_size /= No_of_domains * ik) THEN
    mssg = "File size and number of computed domains does not match. &
        &Please check the boundaries, hexdump the status file and have a look, &
        &if the kind of the integers of this file is correct (INT64/8)."
    CALL print_err_stop(6, mssg, 1_ik)
END IF

OPEN (UNIT=aun, FILE=TRIM(activity_file), ACCESS="STREAM")

READ (aun) Domain_status(:,1)

CLOSE(aun)

!------------------------------------------------------------------------------
! User feedback
!------------------------------------------------------------------------------
WRITE(6, FMT_TXT_AxI0) "size_mpi:      ", size_mpi
WRITE(6, FMT_TXT_AxI0) "No_of_domains: ", No_of_domains
WRITE(6, FMT_SEP)

!------------------------------------------------------------------------------
! Crawl data
!------------------------------------------------------------------------------
par_dmn_number=0
DO rank_mpi = 1, size_mpi-1

    target_path = TRIM(in%p_n_bsnm)//"/"
    file_head = "results"

    pro_path = ""
    pro_name = ""

    WRITE(pro_path,'(A,A,I7.7,A)') TRIM(target_path), "Rank_", rank_mpi, "/"
    WRITE(pro_name,'(A,A,I7.7)')   TRIM(file_head), "_", rank_mpi

    ! Search for valid  directories
    INQUIRE(FILE=pro_path, EXIST=ex)
    IF(.NOT. ex) CYCLE

    dat_no_highscore=0

    par_dmn_number = par_dmn_number + 1

    nn_comm = 1_ik

    !------------------------------------------------------------------------------
    ! Read PureDat tree
    !------------------------------------------------------------------------------
    rank_data = read_tree()

    CALL open_stream_files(rank_data, "read" , "old")

    !------------------------------------------------------------------------------
    ! Read domain numbers
    ! Check that with e.g. 
    ! ../bin/pd_dump_leaf_x86_64 $PWD/FH01-2_mu_Dev_dtc_Tensors/Rank_0002041/ results_0002041 2
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 2_pd_ik, domain_no, success)

    IF(ALLOCATED(dat_domains)) DEALLOCATE(dat_domains)
    IF(dat_no_highscore<domain_no%dat_no) dat_no_highscore = domain_no%dat_no
    ALLOCATE(dat_domains(domain_no%dat_no), stat=alloc_stat)

    dat_no = domain_no%dat_no

    CALL print_err_stop(6, "Allocating 'dat_domains' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, domain_no, dat_domains)
    last_domain_rank = MAXVAL(dat_domains)

    !------------------------------------------------------------------------------
    ! Number of elements
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 3_pd_ik, no_elems, success)

    IF(ALLOCATED(dat_no_elems)) DEALLOCATE(dat_no_elems)
    IF(dat_no_highscore<no_elems%dat_no) dat_no_highscore = no_elems%dat_no
    ALLOCATE(dat_no_elems(no_elems%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_no_elems' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, no_elems, dat_no_elems)

    !------------------------------------------------------------------------------
    ! Number of Nodes
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 4_pd_ik, no_nodes, success)

    IF(ALLOCATED(dat_no_nodes)) DEALLOCATE(dat_no_nodes)
    IF(dat_no_highscore<no_nodes%dat_no) dat_no_highscore = no_nodes%dat_no
    ALLOCATE(dat_no_nodes(no_nodes%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_no_nodes' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, no_nodes, dat_no_nodes)

    !------------------------------------------------------------------------------
    ! Collected logs
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 5_pd_ik, collected_logs, success)

    IF(ALLOCATED(dat_collected_logs)) DEALLOCATE(dat_collected_logs)
    IF(dat_no_highscore<no_nodes%dat_no) dat_no_highscore = no_nodes%dat_no
    ALLOCATE(dat_collected_logs(collected_logs%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_collected_logs' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, collected_logs, dat_collected_logs)

    !------------------------------------------------------------------------------
    ! Start time of the domain
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 6_pd_ik, t_start, success)

    IF(ALLOCATED(dat_t_start)) DEALLOCATE(dat_t_start)
    IF(dat_no_highscore<t_start%dat_no) dat_no_highscore = t_start%dat_no
    ALLOCATE(dat_t_start(t_start%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_t_start' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, t_start, dat_t_start)

    !------------------------------------------------------------------------------
    ! Duration of the domain
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 7_pd_ik, t_duration, success)

    IF(ALLOCATED(dat_t_duration)) DEALLOCATE(dat_t_duration)
    IF(dat_no_highscore<t_duration%dat_no) dat_no_highscore = t_duration%dat_no
    ALLOCATE(dat_t_duration(t_duration%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_t_duration' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, t_duration, dat_t_duration)

    !------------------------------------------------------------------------------
    ! Read effective numerical stiffness
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 9_pd_ik, eff_num_stiffness, success)
        
    IF(ALLOCATED(dat_eff_num_stiffnesses)) DEALLOCATE(dat_eff_num_stiffnesses)
    IF(dat_no_highscore<eff_num_stiffness%dat_no) dat_no_highscore = eff_num_stiffness%dat_no
    ALLOCATE(dat_eff_num_stiffnesses(eff_num_stiffness%dat_no), stat=alloc_stat)

    CALL print_err_stop(6, "Allocating 'dat_eff_num_stiffnesses' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, eff_num_stiffness, dat_eff_num_stiffnesses)

    !------------------------------------------------------------------------------
    ! Read effective density
    !------------------------------------------------------------------------------
    CALL get_leaf_with_num(rank_data, 25_pd_ik, density, success)
    
    IF(ALLOCATED(dat_densities)) DEALLOCATE(dat_densities)
    IF(dat_no_highscore<density%dat_no) dat_no_highscore = density%dat_no
    ALLOCATE(dat_densities(density%dat_no), stat=alloc_stat)
    
    CALL print_err_stop(6, "Allocating 'dat_densities' failed.", alloc_stat)
    CALL pd_read_leaf(rank_data%streams, density, dat_densities)

    DO tt = 1, 3

        !------------------------------------------------------------------------------
        ! Read effective stiffness and euler (!) positions
        !------------------------------------------------------------------------------
        SELECT CASE(tt)
            CASE(1)
                CALL get_leaf_with_num(rank_data, 13_pd_ik, tensors, success)
                local_domain_opt_pos = 0._rk
                fh_tens = fh_covo 
                tensor(tt)%opt_crit = "covo" 
            CASE(2)
                CALL get_leaf_with_num(rank_data, 20_pd_ik, tensors, success)
                CALL get_leaf_with_num(rank_data, 18_pd_ik, leaf_pos, success)
                fh_tens = fh_cr1
                tensor(tt)%opt_crit = "cr1" 
            CASE(3)
                CALL get_leaf_with_num(rank_data, 24_pd_ik, tensors, success)
                CALL get_leaf_with_num(rank_data, 22_pd_ik, leaf_pos, success)
                fh_tens = fh_cr2
                tensor(tt)%opt_crit = "cr2" 
        END SELECT


        !------------------------------------------------------------------------------
        ! User feedback
        !------------------------------------------------------------------------------
        WRITE(6, FMT_TXT_xAI0) &
            "Crawling "//TRIM(tensor(tt)%opt_crit)//" of DTC rank ", &
            rank_mpi, "       ", par_dmn_number

        IF(ALLOCATED(dat_tensors)) DEALLOCATE(dat_tensors)
        ALLOCATE(dat_tensors(tensors%dat_no), stat=alloc_stat)
    
        WRITE(tt_char, '(I0)') tt
        CALL print_err_stop(6, "Allocating 'dat_tensors("//tt_char//")' failed.", alloc_stat)
        CALL pd_read_leaf(rank_data%streams, tensors, dat_tensors)

        !------------------------------------------------------------------------------
        ! Rotation Vector CR_1
        !------------------------------------------------------------------------------
        IF(tt/= 1) THEN
            IF(ALLOCATED(dat_pos)) DEALLOCATE(dat_pos)
            ALLOCATE(dat_pos(leaf_pos%dat_no), stat=alloc_stat)
            CALL print_err_stop(6, "Allocating 'dat_pos("//tt_char//")' failed.", alloc_stat)
            CALL pd_read_leaf(rank_data%streams, leaf_pos, dat_pos)
        END IF 

        !------------------------------------------------------------------------------
        ! Begin searching for the domain to extract
        !------------------------------------------------------------------------------
        last_domain = .FALSE.
        DO ii = 1, dat_no 

            !------------------------------------------------------------------------------
            ! Quite a naive implementation
            !------------------------------------------------------------------------------
            ! DO jj =1, No_of_domains
            !------------------------------------------------------------------------------
            ! Reset the information stored before proceeding to the next one.
            ! Iterating over jj while searching for the actual domain number stored in 
            ! dat_domains implicitly sorts the data. 
            !------------------------------------------------------------------------------
            tensor(tt)%dmn            = 0_ik
            tensor(tt)%no_elems       = 0_ik
            tensor(tt)%no_nodes       = 0_ik
            tensor(tt)%collected_logs = 0_ik
            tensor(tt)%bvtv           = 0._rk
            tensor(tt)%doa_zener      = 0._rk
            tensor(tt)%doa_gebert     = 0._rk
            tensor(tt)%sym            = 0._rk
            tensor(tt)%mat            = 0._rk
            tensor(tt)%num8           = 0._rk
            tensor(tt)%num20          = 0._rk
            tensor(tt)%pos            = 0._rk
            tensor(tt)%sym            = 0._rk
            tensor(tt)%mps            = 0._rk
            tensor(tt)%t_start        = 0._rk
            tensor(tt)%t_duration     = 0._rk

            tensor(tt)%mi_el_type = "" 
            tensor(tt)%ma_el_order = 0_ik
            !------------------------------------------------------------------------------
            ! This value is hardcoded in ./f-src/mod_struct_calcmat.f90
            !------------------------------------------------------------------------------
            tensor(tt)%opt_res = 1._rk 

            !------------------------------------------------------------------------------
            ! Search for the current domain
            !------------------------------------------------------------------------------
            current_domain = dat_domains(ii)

            pntr = -1
            DO jj = 1, No_of_domains
                IF (Domain_status(jj,1) == current_domain) pntr = jj
            END DO

            !------------------------------------------------------------------------------
            ! If, e.g. 500 domains are computed in parallel, but only 800 are requested
            ! for computation, then cycle the last 200.
            !------------------------------------------------------------------------------
            IF (pntr == -1) CYCLE

            !------------------------------------------------------------------------------
            ! Check whether the domain was already found.
            ! tt 1...3 for different criteria
            !------------------------------------------------------------------------------
            IF(Domain_status(pntr, tt+1) == 0_ik) THEN
                Domain_status(pntr, tt+1) = 1_ik
            ELSE
                CYCLE
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
            ! Local numerical stiffness
            !------------------------------------------------------------------------------
            local_num_tensor = 0._rk
            mm = (ii-1) * no_lc * no_lc + 1
            DO kk = 1, no_lc
            DO ll = 1, no_lc
                local_num_tensor(kk, ll) = dat_eff_num_stiffnesses(mm)
                mm = mm + 1
            END DO 
            END DO

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
                CALL write_matrix(6, "Effective stiffness of domain "//dmn_char, &
                    local_domain_tensor, fmti='std', unit='MPa')
            END IF 

            !------------------------------------------------------------------------------
            ! Check symmetry of the tensor
            !------------------------------------------------------------------------------
            sym = check_sym (local_domain_tensor)

            !------------------------------------------------------------------------------
            ! Compute additional parameters.
            !------------------------------------------------------------------------------
            tensor(tt)%dmn  = dat_domains(ii)
            tensor(tt)%bvtv = dat_densities(ii)
            tensor(tt)%no_elems = dat_no_elems(ii)
            tensor(tt)%no_nodes = dat_no_nodes(ii)
            tensor(tt)%t_start    = dat_t_start(ii)
            tensor(tt)%t_duration = dat_t_duration(ii)

            tensor(tt)%sym  = sym
            tensor(tt)%mat  = local_domain_tensor

            IF (ma_el_order == 1) THEN
                tensor(tt)%num8  = local_num_tensor

            ELSE IF (ma_el_order == 2) THEN
                tensor(tt)%num20 = local_num_tensor
            END IF 

            tensor(tt)%pos  = local_domain_opt_pos
            tensor(tt)%dmn_size = dmn_size(1)

            grid = get_grid(dmn_size, vdim, spcng)
            tensor(tt)%section = domain_no_to_section(tensor(tt)%dmn, grid)
            
            DO xx=1, 3
                tensor(tt)%phy_dmn_bnds(xx,1) = &
                        tensor(tt)%section(xx) * tensor(tt)%dmn_size

                tensor(tt)%phy_dmn_bnds(xx,2) = &
                    (tensor(tt)%section(xx)+1) * tensor(tt)%dmn_size
            END DO

            tensor(tt)%sym = check_sym(tensor(tt)%mat)
            tensor(tt)%mps = mps(tensor(tt)%mat)
            tensor(tt)%doa_zener = doa_zener(tensor(tt)%mat)
            tensor(tt)%doa_gebert = doa_gebert(tensor(tt)%mat)
            
            !------------------------------------------------------------------------------
            ! Could be assigned globally, but for the sake of clarity, the element type
            ! gets assigned at this place.
            !------------------------------------------------------------------------------
            tensor(tt)%mi_el_type = TRIM(mi_el_type) 
            tensor(tt)%ma_el_order = ma_el_order

            !------------------------------------------------------------------------------
            ! This value is hardcoded in ./f-src/mod_struct_calcmat.f90
            ! The 24 dat_collected entries are not related to 24 load cases/dof !!
            !------------------------------------------------------------------------------
            tensor(tt)%opt_res = 1._rk 

            DO xx=1, 24
                tensor(tt)%collected_logs(xx) = dat_collected_logs(24_ik*(ii-1)+xx)
            END DO
            
            CALL write_tensor_2nd_rank_R66_row(tensor(tt), string)
            WRITE(fh_tens,'(A)') TRIM(string)


            IF (tt==1) THEN
                CALL write_eff_num_stiff_row(tensor(tt), no_lc, string)
                WRITE(fh_covo_num,'(A)') TRIM(string)
            END IF
        END DO ! ii = 1, dat_no_highscore

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
CALL meta_stop_ascii(fh_covo_num, covo_num_suf)

WRITE(6, FMT_TXT_SEP)
WRITE(6, FMT_TXT_xAF0) "Program finished successfully in ", end-start, " seconds."
WRITE(6, FMT_TXT_SEP)

IF(6 /= 6) CALL meta_stop_ascii(fh=6, suf='.6')

END PROGRAM mrd_crawl_tensors

