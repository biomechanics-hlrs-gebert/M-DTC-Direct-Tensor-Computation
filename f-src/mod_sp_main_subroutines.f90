!------------------------------------------------------------------------------
! MODULE: exec_single_domain
!------------------------------------------------------------------------------
!> @author Ralf Schneider  - HLRS - NUM - schneider@hlrs.de
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! DESCRIPTION: 
!> Auxiliary, but central subroputines of the struct-process/DTC
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on: 29.03.2022
!------------------------------------------------------------------------------
MODULE dtc_main_subroutines

USE global_std
USE Operating_System
USE puredat_com
USE chain_routines
USE linfe
USE mpi
USE gen_geometry
USE PETSC
USE calcmat
USE system
USE mpi_system

IMPLICIT NONE 

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: exec_single_domain
!------------------------------------------------------------------------------  
!> @author Ralf Schneider  - HLRS - NUM - schneider@hlrs.de
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Execution chain for the treatment of a single MVE
!
!> @param[inout] root          Root tree of computation
!> @param[in]    lin_domain    Linear, internal, domain number
!> @param[in]    domain        ddc domain number
!> @param[in]    job_dir       Job directory (Rank_...)
!> @param[out]   active        Whether thread is active
!> @param[in]    fh_mpi_worker File handle of sub comm
!> @param[in]    worker_rank_mpi
!> @param[in]    worker_size_mpi
!> @param[in]    comm_mpi
!------------------------------------------------------------------------------
Subroutine exec_single_domain(root, comm_nn, domain, typeraw, &
    job_dir, fh_cluster_log, active, fh_mpi_worker, worker_rank_mpi, worker_size_mpi, comm_mpi)

TYPE(materialcard) :: bone, bone_adjusted
INTEGER(mik), Intent(INOUT), Dimension(no_streams) :: fh_mpi_worker

Character(*) , Intent(in)  :: job_dir
Character(*) , Intent(in)  :: typeraw
INTEGER(mik), Intent(In)  :: worker_rank_mpi, worker_size_mpi, comm_mpi
INTEGER(ik) , intent(in)  :: comm_nn, domain, fh_cluster_log
INTEGER(mik), intent(out) :: active
Type(tBranch)    , Intent(inOut) :: root

REAL(rk), DIMENSION(:), Pointer     :: displ, force
REAL(rk), DIMENSION(:), Allocatable :: glob_displ, glob_force, zeros_R8
REAL(rk) :: factor

INTEGER(mik), Dimension(MPI_STATUS_SIZE) :: status_mpi
INTEGER(mik) :: ierr, petsc_ierr, result_len_mpi_procs

INTEGER(pd_ik), Dimension(:), Allocatable :: serial_pb
INTEGER(pd_ik) :: serial_pb_size

INTEGER(ik) :: domain_elems, ii, jj, kk, id, stat, &
    Istart,Iend, parts, IVstart, IVend, m_size, mem_global, status_global, &
    no_different_hosts, timestamp, HU, lower_limit, upper_limit, &
    macro_order, no_elem_nodes, no_lc

INTEGER(ik), DIMENSION(24) :: collected_logs ! timestamps, memory_usage, pid_returned
INTEGER(ik), Dimension(:), Allocatable :: nodes_in_mesh, gnid_cref
INTEGER(ik), Dimension(:,:), Allocatable :: res_sizes

INTEGER(c_int) :: stat_c_int

CHARACTER(9)   :: domain_char
CHARACTER(40)  :: mssg_fix_len
CHARACTER(mcl) :: timer_name, rank_char, domain_desc, part_desc, &
    desc, mesh_desc, filename, elt_micro, map_sgmnttn=""
CHARACTER(1) :: bin_sgmnttn=""

Character, Dimension(4*mcl) :: c_char_array
Character, Dimension(:), Allocatable :: char_arr

INTEGER(ik), PARAMETER :: host_name_length = 512
CHARACTER(host_name_length), DIMENSION(:), ALLOCATABLE :: host_list, unique_host_list
CHARACTER(host_name_length) :: host_of_part

LOGICAL, PARAMETER :: DEBUG = .TRUE.
logical :: success=.TRUE., host_assumed_unique

Type(tBranch), pointer :: boundary_branch, domain_branch, part_branch
Type(tBranch), pointer :: mesh_branch, meta_para, esd_result_branch

Type(tMat)         :: AA, AA_org
Type(tVec)         :: XX
TYPE(tPETScViewer) :: PetscViewer
Type(tKSP)         :: KSP

INTEGER(ik), Dimension(60)    :: idxm_20, idxn_20
Real(rk),    Dimension(60,60) :: K_loc_20
Type(tVec), Dimension(60) :: FF_20

INTEGER(ik), Dimension(24)    :: idxm_08, idxn_08
Real(rk),    Dimension(24,24) :: K_loc_08
Type(tVec), Dimension(24) :: FF_08

! Init worker_is_active status
Active = 0_mik
collected_logs = 0_ik

write(domain_char,'(I0)') domain

!--------------------------------------------------------------------------
! Get basic infos 
!--------------------------------------------------------------------------
host_of_part = ""
CALL mpi_get_processor_name(host_of_part, result_len_mpi_procs, ierr) 
CALL print_err_stop(std_out, "mpi_get_processor_name failed", ierr)


!------------------------------------------------------------------------------
! Get the mpi communicators total memory usage by the pids of the threads.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

!------------------------------------------------------------------------------
! Tracking (the memory usage is) intended for use during production too.
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN
    collected_logs(1) = INT(time(), ik)
    collected_logs(8) = mem_global
    collected_logs(15) = status_global
END IF


! This function does not work out well.
! CALL get_environment_Variable("HOST", host_of_part)

CALL Search_branch("Input parameters", root, meta_para, success, out_amount)

CALL pd_get(meta_para, "No of mesh parts per subdomain", parts)
CALL pd_get(meta_para, "Physical domain size" , bone%phdsize, 3)
CALL pd_get(meta_para, "Grid spacings" , bone%delta, 3)
CALL pd_get(meta_para, "Young_s modulus" , bone%E)
CALL pd_get(meta_para, "Poisson_s ratio" , bone%nu)

!------------------------------------------------------------------------------
! Rank = 0 -- Local master of comm_mpi
!------------------------------------------------------------------------------
If (worker_rank_mpi == 0) then

    !------------------------------------------------------------------------------
    ! Create job directory in case of non-production run
    !------------------------------------------------------------------------------
    If (out_amount /= "PRODUCTION") then

        c_char_array(1:len(Trim(job_dir)//Char(0))) = str_to_char(Trim(job_dir)//Char(0))

        CALL Stat_Dir(c_char_array, stat_c_int)

        If(stat_c_int /= 0) Then

            CALL exec_cmd_line("mkdir -p "//TRIM(job_dir), stat)

            IF(stat /= 0) CALL print_err_stop(std_out, &
                "Couldn't execute syscall mkdir -p "//TRIM(job_dir), 1)

            CALL Stat_Dir(c_char_array, stat_c_int)

            IF(stat_c_int /= 0) CALL print_err_stop(std_out, &
                "Couldn't create directory "//TRIM(job_dir), 1)

        End If

    End If

    !------------------------------------------------------------------------------
    ! Write log and monitor file
    !------------------------------------------------------------------------------
    Write(un_lf,FMT_MSG_SEP)
    timestamp = time()

    
    WRITE(un_lf, '(A,I0)') 'Start time: ', timestamp

    Write(un_lf, FMT_MSG_xAI0) "Domain No.: ", domain
    Write(un_lf, FMT_MSG)      "Job_dir:    "//Trim(job_dir)
    Write(un_lf,FMT_MSG_SEP)

    !------------------------------------------------------------------------------
    ! Get and log the hosts of all parts
    !------------------------------------------------------------------------------
    IF (.NOT. ALLOCATED(host_list)) ALLOCATE(host_list(parts))
    IF (.NOT. ALLOCATED(unique_host_list)) ALLOCATE(unique_host_list(parts))
    
    host_list = ""
    host_list(1) = host_of_part

    no_different_hosts = 1

    unique_host_list = host_list

    !------------------------------------------------------------------------------
    ! Loop over all ranks of the sub comm
    ! parts refers to ranks an therefore is counted from o, 
    ! but 0 is the master and does not need so send to itself.
    !------------------------------------------------------------------------------
    Do ii = 1, parts-1 

        IF (ii/=0) THEN
            CALL MPI_RECV(host_list(ii), INT(host_name_length, mik), MPI_CHAR, &
                INT(ii,mik), INT(ii,mik), comm_mpi, status_mpi, ierr)

            WRITE(rank_char, "(I20)") ii
            CALL print_err_stop(std_out, &
                "MPI_RECV on host_of_part "//TRIM(ADJUSTL(rank_char))//" didn't succseed", ierr)
        END IF

        !------------------------------------------------------------------------------
        ! Test whether the host is unique in the domain
        !------------------------------------------------------------------------------
        IF (host_list(ii) == "") CYCLE
        
        host_assumed_unique = .TRUE.
        DO jj = 1, parts

            IF((host_list(ii) == unique_host_list(jj)) .AND. (unique_host_list(jj) /= "")) THEN
                host_assumed_unique = .FALSE.
            END IF 
        END DO
        
        IF (host_assumed_unique) THEN
            no_different_hosts = no_different_hosts + 1_ik
            unique_host_list(no_different_hosts) = host_list(ii)
        END IF 

    END Do

    !------------------------------------------------------------------------------
    ! Concatenate the string of the unique host names.
    !------------------------------------------------------------------------------
    mssg_fix_len = "Unique hosts of domain:"
    WRITE(fh_cluster_log, '(A)', ADVANCE="NO") mssg_fix_len//", "//domain_char//", " 

    DO ii = 1, no_different_hosts-1

        WRITE(fh_cluster_log, '(A)', ADVANCE="NO") TRIM(unique_host_list(ii))//","
    END DO
    WRITE(fh_cluster_log, '(A)', ADVANCE="NO") TRIM(unique_host_list(no_different_hosts))

    WRITE(fh_cluster_log, '(A)')
    FLUSH(fh_cluster_log)

ELSE ! (worker_rank_mpi == 0) THEN

    CALL MPI_SEND(host_of_part, INT(host_name_length, mik), MPI_CHAR, 0_mik, worker_rank_mpi, comm_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_SEND of host_of_part didn't succeed", ierr)

END IF 



If (worker_rank_mpi == 0) then

    !------------------------------------------------------------------------------
    ! Generate Geometry
    !------------------------------------------------------------------------------
    Select Case (timer_level)
    Case (1)
        timer_name = "+-- generate_geometry "//trim(domain_char)
    Case default
        timer_name = "generate_geometry"
    End Select

    CALL start_timer(trim(timer_name), .FALSE.)

    CALL generate_geometry(root, domain, job_dir, typeraw, success)

    if (.not. success) write(std_out, FMT_WRN)"generate_geometry() failed."
    
    CALL end_timer(trim(timer_name))
        timestamp = time()


    WRITE(un_lf, '(A,I0)') 'End time: ', timestamp

    !------------------------------------------------------------------------------
    ! Look for the Domain branch
    !------------------------------------------------------------------------------
    domain_desc=''
    Write(domain_desc,'(A,I0)')'Domain ', domain
    
    IF(success) CALL search_branch(trim(domain_desc), root, domain_branch, success, out_amount)

    !------------------------------------------------------------------------------
    ! Get the no of nodes per part
    !------------------------------------------------------------------------------
    mesh_desc = ''
    Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_', domain
    
    IF(success) CALL search_branch(trim(mesh_desc), domain_branch, mesh_branch, success, out_amount)

    !------------------------------------------------------------------------------
    ! Send success variable for gracefully aborting the domain.
    !------------------------------------------------------------------------------
    CALL mpi_bcast(success, 1_mik, MPI_LOGICAL, 0_mik, COMM_MPI, ierr)

    IF (.NOT. success) GOTO 1000

    !------------------------------------------------------------------------------
    ! Only read from the mesh branch if the branch was found.
    ! Otherwise, gracefully abort the domain to proceed with the next one.
    !------------------------------------------------------------------------------
    IF(success) THEN
        CALL pd_get(mesh_branch, 'No of nodes in mesh',  nodes_in_mesh)

        !------------------------------------------------------------------------------
        ! Set the global matrix size
        !------------------------------------------------------------------------------
        m_size = nodes_in_mesh(1) * 3
    ELSE
        m_size = 0
    END IF 

    Do ii = 1, parts-1

        !------------------------------------------------------------------------------
        ! Look for the Part branch
        !------------------------------------------------------------------------------
        domain_desc=''
        Write(part_desc,'(A,I0)')'Part_',ii
        if(success) CALL search_branch(trim(part_desc), domain_branch, part_branch, success)

        If (.NOT. success) Then
            WRITE(mssg,'(A,I0,A,L,A,I0,A)') "Something bad and unexpected happend &
            &in exec_single_domain! Looking for branch of part ",ii," returned ", &
                success, "MPI proc ",worker_rank_mpi," halted."
            CALL print_err_stop(std_out, mssg, 0)
        End If

        !------------------------------------------------------------------------------
        ! Serialize branch with mesh part to send via mpi
        !------------------------------------------------------------------------------
        CALL serialize_branch(part_branch, serial_pb, serial_pb_size, .TRUE.)

        CALL mpi_send(serial_pb_size, 1_mik, MPI_INTEGER8, Int(ii,mik), Int(ii,mik), &
            COMM_MPI, ierr)
        CALL mpi_send(serial_pb, INT(serial_pb_size,mik), MPI_INTEGER8, &
            Int(ii,mik), Int(ii,mik), COMM_MPI, ierr)

        Deallocate(serial_pb)
        
    End Do

    part_desc=''
    Write(part_desc,'(A,I0)')'Part_', parts
    CALL search_branch(trim(part_desc), domain_branch, part_branch, success, out_amount)

    !------------------------------------------------------------------------------
    ! Broadcast matrix size. 
    ! TODO could also be included into part branches.
    !------------------------------------------------------------------------------
    CALL mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)

!------------------------------------------------------------------------------
! Ranks > 0 - Workers
!------------------------------------------------------------------------------
Else

    !------------------------------------------------------------------------------
    ! Check if the mesh branch was read successfully, if not - abort the domain.
    !------------------------------------------------------------------------------
    CALL mpi_bcast(success, 1_mik, MPI_LOGICAL, 0_mik, COMM_MPI, ierr)
    IF(.NOT. success) GOTO 1000

    CALL mpi_recv(serial_pb_size, 1_mik, MPI_INTEGER8, 0_mik, &
        worker_rank_mpi, COMM_MPI, status_mpi, ierr)

    if (allocated(serial_pb)) deallocate(serial_pb)
    
    Allocate(serial_pb(serial_pb_size))

    CALL mpi_recv(serial_pb, INT(serial_pb_size,mik), MPI_INTEGER8, &
        0_mik, worker_rank_mpi, COMM_MPI, status_mpi, ierr)

    !------------------------------------------------------------------------------
    ! Deserialize part branch
    !------------------------------------------------------------------------------
    CALL Start_Timer("Deserialize part branch branch")

    Allocate(part_branch)
    
    CALL deserialize_branch(part_branch, serial_pb, .TRUE.)

    CALL End_Timer("Deserialize part branch branch")

    CALL mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)
            
End If ! (worker_rank_mpi == 0) then

!------------------------------------------------------------------------------
! Abort this domain in case of a fatal error in reading the mesh branch.
!------------------------------------------------------------------------------
IF(.NOT. success) write(*,*) "My process failed, rank:", worker_rank_mpi


IF(.NOT. success) GOTO 1000

!------------------------------------------------------------------------------
! Setup the linear System with a constant system matrix A. Once that
! is done setup the multiple right hand sides and solve the linear
! system multiple times. Save the solutions to calculate effective
! stiffness matirces.
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN   ! Sub Comm Master

    If (out_amount == "DEBUG") THEN 
        Write(un_lf, fmt_dbg_sep)
        Write(un_lf, '(A)') "part branch right after deserialization"
        CALL log_tree(part_branch, un_lf, .TRUE.)
        Write(un_lf, fmt_dbg_sep)
        flush(un_lf)
    END If

    SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- create_Stiffness_matrix "//TRIM(domain_char)
        CASE default
            timer_name = "create_Stiffness_matrix"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 

!------------------------------------------------------------------------------
! Get the mpi communicators total memory usage by the pids of the threads.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

!------------------------------------------------------------------------------
! Tracking (the memory usage is) intended for use during production too.
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN
    collected_logs(2) = INT(time(), ik)
    collected_logs(9) = mem_global
    collected_logs(16) = status_global
END IF
            
!------------------------------------------------------------------------------
! Create Stiffness matrix
! Preallocation avoids dynamic allocations during matassembly.
!------------------------------------------------------------------------------
CALL MatCreate(COMM_MPI, AA    , petsc_ierr)
CALL MatCreate(COMM_MPI, AA_org, petsc_ierr)

CALL MatSetFromOptions(AA,     petsc_ierr)
CALL MatSetFromOptions(AA_org, petsc_ierr)

CALL MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)
CALL MatSetSizes(AA_org,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)

! https://lists.mcs.anl.gov/pipermail/petsc-users/2021-January/042972.html
CALL MatSeqAIJSetPreallocation(AA, 0_ik, PETSC_NULL_INTEGER, petsc_ierr)
CALL MatMPIAIJSetPreallocation(AA, 0_ik, PETSC_NULL_INTEGER, 0_ik, PETSC_NULL_INTEGER, petsc_ierr)

CALL MatSeqAIJSetPreallocation(AA_org, 0_ik, PETSC_NULL_INTEGER, petsc_ierr)
CALL MatMPIAIJSetPreallocation(AA_org, 0_ik, PETSC_NULL_INTEGER, 0_ik, PETSC_NULL_INTEGER, petsc_ierr)

CALL MatSetOption(AA    ,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,petsc_ierr)
CALL MatSetOption(AA_org,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,petsc_ierr)


!------------------------------------------------------------------------------
! Get and write another memory log.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

IF (worker_rank_mpi == 0) THEN
    collected_logs(3) = INT(time(), ik)
    collected_logs(10) = mem_global
    collected_logs(17) = status_global
END IF


!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 
    CALL MatGetOwnershipRange(AA, Istart, Iend, petsc_ierr)
    Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
        "MPI rank : ",worker_rank_mpi,"| matrix ownership for ","A    ", ":",Istart," -- " , Iend

    CALL MatGetOwnershipRange(AA_org, Istart, Iend, petsc_ierr)

    Write(un_lf, "('MM ', A,I4,A,A6,2(A,I9))") &
    "MPI rank : ", worker_rank_mpi, "| matrix ownership for ", "A_org", ":",Istart, " -- " , Iend
End If

!------------------------------------------------------------------------------
! part_branch%leaves(1) : Local  node ids   in part branch
! part_branch%leaves(2) : Coordinate values in part branch
! part_branch%leaves(3) : Global node ids   in part branch
! part_branch%leaves(4) : NOT USED ! (For Element Numbers)
! part_branch%leaves(5) : Topology
! part_branch%leaves(6) : Hounsfield units of the voxels (HU Magnitude)
! No Cross reference global nid to local nid is currently needed since 
! renumbering is deactivated in mod_mesh_partitioning.f90 : part_mesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Segmentation specific information
!------------------------------------------------------------------------------
CALL pd_get(root%branches(1), 'Binary segmentation', char_arr)
bin_sgmnttn = char_to_str(char_arr)

CALL pd_get(root%branches(1), 'Segmentation map', char_arr)
map_sgmnttn = char_to_str(char_arr)

CALL pd_get(root%branches(1), 'Element type  on micro scale', char_arr)
elt_micro = char_to_str(char_arr)

CALL pd_get(root%branches(1), 'Element order on macro scale', macro_order)

!------------------------------------------------------------------------------
! Set macro element specific stuff
!------------------------------------------------------------------------------
IF (macro_order == 1) THEN
    no_elem_nodes = 8
    no_lc = 24
ELSE IF (macro_order == 2) THEN
    no_elem_nodes = 20
    no_lc = 60
ELSE
    CALL print_err_stop(std_out, "Element orders other than 1 or 2 are not supported", 1)
END IF

bone_adjusted%nu = bone%nu

!------------------------------------------------------------------------------
! Assemble matrix
!------------------------------------------------------------------------------
if (TRIM(elt_micro) == "HEX08") then

    domain_elems = part_branch%leaves(5)%dat_no / 8

    Do ii = 1, domain_elems

        kk = 1
        
        ! Translate Topology to global dofs
        Do jj = (ii-1)*8+1 , ii*8
            
            id = part_branch%leaves(5)%p_int8(jj)
            id = (id - 1) * 3

            idxm_08(kk)   = (id    ) !*cr
            idxm_08(kk+1) = (id + 1) !*cr
            idxm_08(kk+2) = (id + 2) !*cr
            
            kk = kk + 3
        End Do
                        
        IF (bin_sgmnttn == "Y") THEN

            HU = part_branch%leaves(6)%p_int8(ii)

            !------------------------------------------------------------------------------
            ! Need to distinguish between > and < 1000 HU. Otherwise, the young modulus
            ! may be negative, which leads to faulty results.
            !------------------------------------------------------------------------------
    
            !------------------------------------------------------------------------------
            ! FH01-4_cl
            !------------------------------------------------------------------------------
            lower_limit = 250_ik
            upper_limit = 14250_ik
    
            !------------------------------------------------------------------------------
            ! FH01-1_mu_121212
            !------------------------------------------------------------------------------
            ! lower_limit = 22000_ik
            ! upper_limit = 35000_ik
    
            !------------------------------------------------------------------------------
            ! FH01-1_mu_121240
            !------------------------------------------------------------------------------
            ! lower_limit = 22000_ik
            ! upper_limit = 35000_ik
    
    
            IF (HU > lower_limit) THEN
                
                IF (HU < 400_ik) THEN
                    factor = 0.2_rk - REAL(400_ik-HU, rk)*0.001_rk
                    ! factor = 0.2_rk - REAL(600_ik-HU, rk)*0.00025_rk
                ELSE
                    factor = 0.2_rk + 0.000057762_rk * REAL(HU-400_ik, rk)
                    ! factor = 0.0_rk 
                END IF 
                
                !------------------------------------------------------------------------------
                ! FH01-4_cl
                !------------------------------------------------------------------------------
                ! factor = 1._rk / REAL((upper_limit-lower_limit),rk) * REAL((HU - lower_limit), rk)
                K_loc_08 = Hexe08_FH01_4_cl(bone) * factor
                !
                !------------------------------------------------------------------------------
                ! FH01-1_mu_121212
                !------------------------------------------------------------------------------
                ! factor = 1._rk / REAL((upper_limit-lower_limit),rk) * REAL((HU - lower_limit), rk) * 0.5
                ! K_loc_08 = Hexe08_5600_03_121212(bone) * factor
                !
                !------------------------------------------------------------------------------
                ! FH01-1_mu_121240
                !------------------------------------------------------------------------------
                ! factor = 1._rk / REAL((upper_limit-lower_limit),rk) * REAL((HU - lower_limit), rk) * 0.5
                ! K_loc_08 = Hexe08_5600_03_121240(bone) * factor
                ! 
                ! K_loc_08 = Hexe08(bone)
                ! K_loc_08 = Hexe08_5600_03_121212(bone)
                ! K_loc_08 = Hexe08_5600_03_121240(bone)
    
                ! K_loc_08 = K_loc_08 * factor
    
            ELSE
                bone_adjusted%E = 0._rk
                ! K_loc_08 = Hexe08(bone_adjusted)
                
                ! K_loc_08 = Hexe08(bone)
                ! K_loc_08 = Hexe08_5600_03_121212(bone)
                ! K_loc_08 = Hexe08_5600_03_121240(bone)
                K_loc_08 = Hexe08_FH01_4_cl(bone_adjusted)
                
            END IF 
        ELSE ! (bin_sgmnttn == "N") THEN
            K_loc_08 = Hexe08(bone_adjusted)
        END IF 

        idxn_08 = idxm_08

        CALL MatSetValues(AA, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)

        CALL MatSetValues(AA_org, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)
        
    End Do

else if (elt_micro == "HEX20") then
    
    K_loc_20 = Hexe20(bone)

    domain_elems = part_branch%leaves(5)%dat_no / 20

    Do ii = 1, domain_elems

        kk = 1
        
        ! Translate Topology to global dofs
        Do jj = (ii-1)*20+1 , ii*20
            
            id = part_branch%leaves(5)%p_int8(jj)
            id = (id - 1) * 3
            
            idxm_20(kk)   = (id    ) !*cr
            idxm_20(kk+1) = (id + 1) !*cr
            idxm_20(kk+2) = (id + 2) !*cr
            
            kk = kk + 3

        End Do

        IF (bin_sgmnttn == "Y") THEN

            HU = part_branch%leaves(6)%p_int8(ii)

            IF (HU > 1000_ik) THEN
                factor = 1._rk / REAL((14250_ik-1000_ik),rk) * REAL((HU - 1000_ik), rk)
                    
                K_loc_20 = Hexe20(bone)
    
                K_loc_20 = K_loc_20 * factor
            ELSE
                bone_adjusted%E = 0._rk
                K_loc_20 = Hexe20(bone_adjusted)
            END IF 
    
        ELSE
            K_loc_20 = Hexe20(bone_adjusted)
        END IF 
                
        idxn_20 = idxm_20

        CALL MatSetValues(AA, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
        CALL MatSetValues(AA_org, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
        
    End Do

end if


IF (worker_rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- MatAssemblyBegin "//TRIM(domain_char)
    CASE default
        timer_name = "MatAssemblyBegin"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 

!------------------------------------------------------------------------------
! Get and write another memory log.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

IF (worker_rank_mpi == 0) THEN
    collected_logs(4) = INT(time(), ik)
    collected_logs(11) = mem_global
    collected_logs(18) = status_global
END IF

CALL MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
CALL MatAssemblyBegin(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)
! Computations can be done while messages are in transition
CALL MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
CALL MatAssemblyEnd(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) CALL end_timer(TRIM(timer_name))

!------------------------------------------------------------------------------
! Crashes on hawk. Out of memory?
!------------------------------------------------------------------------------
!  If (out_amount == "DEBUG") THEN 
!     CALL PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
!     CALL PetscViewerASCIIOpen(COMM_MPI,"AA.output.1",PetscViewer, petsc_ierr);
!     CALL PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
!     CALL MatView(AA, PetscViewer, petsc_ierr)
!     CALL PetscViewerDestroy(PetscViewer, petsc_ierr)
!  End If
    

!------------------------------------------------------------------------------
! At this point the system matrix is assembled. To make it ready to be
! used, the rows and columns of the dofs with prescribed displacements
! have to be eliminated. To do that with MatZeroRowsColumns we need
! right hand side vectors and a solution vector.
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- create_rh_solution_vetors "//TRIM(domain_char)
    CASE default
        timer_name = "create_rh_solution_vetors"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 


!------------------------------------------------------------------------------
! Get and write another memory log.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

IF (worker_rank_mpi == 0) THEN
    collected_logs(5) = INT(time(), ik)
    collected_logs(12) = mem_global
    collected_logs(19) = status_global
END IF

Do ii = 1, no_lc

    !------------------------------------------------------------------------------
    ! Create load vectors
    !------------------------------------------------------------------------------
    IF (macro_order == 1) THEN
        CALL VecCreate(COMM_MPI, FF_08(ii), petsc_ierr)
        CALL VecSetSizes(FF_08(ii), PETSC_DECIDE, m_size, petsc_ierr)
        CALL VecSetFromOptions(FF_08(ii), petsc_ierr)
        CALL VecSet(FF_08(ii), 0._rk,petsc_ierr)
        
        CALL VecGetOwnershipRange(FF_08(ii), IVstart, IVend, petsc_ierr)
    ELSE IF (macro_order == 2) THEN
        CALL VecCreate(COMM_MPI, FF_20(ii), petsc_ierr)
        CALL VecSetSizes(FF_20(ii), PETSC_DECIDE, m_size, petsc_ierr)
        CALL VecSetFromOptions(FF_20(ii), petsc_ierr)
        CALL VecSet(FF_20(ii), 0._rk,petsc_ierr)
        
        CALL VecGetOwnershipRange(FF_20(ii), IVstart, IVend, petsc_ierr)
    END IF

    If (out_amount == "DEBUG") THEN 
        Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
            "MPI rank : ",worker_rank_mpi,"| vector ownership for ","F", ":",IVstart," -- " , IVend
    End If
    
    !------------------------------------------------------------------------------
    IF (macro_order == 1) THEN
        CALL VecAssemblyBegin(FF_08(ii), petsc_ierr)

    ELSE IF (macro_order == 2) THEN
        CALL VecAssemblyBegin(FF_20(ii), petsc_ierr)

    END IF
    
    ! Computations can be done while messages are transferring
End Do

Do ii = 1, no_lc 
    IF (macro_order == 1) THEN
        CALL VecAssemblyEnd(FF_08(ii), petsc_ierr)

    ELSE IF (macro_order == 2) THEN
        CALL VecAssemblyEnd(FF_20(ii), petsc_ierr)

    END IF
    
End Do
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Get Bounds branch of LC 1
! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
!------------------------------------------------------------------------------ 
write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",1
CALL search_branch(trim(desc), part_branch, boundary_branch, success, out_amount)

!------------------------------------------------------------------------------ 
! Setup id reference vector for displacement insertion
! boundary_branch%leaves(2)%dat_no = Number of constraint nodes * 3
!------------------------------------------------------------------------------ 
Allocate(gnid_cref(boundary_branch%leaves(2)%dat_no), zeros_R8(boundary_branch%leaves(2)%dat_no))
zeros_R8 = 0._rk

kk = 1
Do ii = 1, boundary_branch%leaves(1)%dat_no
    id = boundary_branch%leaves(1)%p_int8(ii)
    id = (id - 1)  * 3
    gnid_cref(kk:kk+2) = [id, id+1, id+2]
    kk= kk + 3
End Do

!------------------------------------------------------------------------------ 
! Create solution vector
!------------------------------------------------------------------------------ 
CALL VecCreate(COMM_MPI, XX, petsc_ierr)
CALL VecSetSizes(XX, PETSC_DECIDE, m_size, petsc_ierr)
CALL VecSetFromOptions(XX, petsc_ierr)
CALL VecSet(XX, 0._rk,petsc_ierr)

!------------------------------------------------------------------------------
! Only set prescribed displacements if there are boundary nodes available.
! These are calculated in struct_preprocess subroutine generate_boundaries
!------------------------------------------------------------------------------
IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
    ! Set prescribed displacements of LC1 to solution vector
    CALL VecSetValues(XX, boundary_branch%leaves(2)%dat_no, &
        gnid_cref, -boundary_branch%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
END IF 

CALL VecAssemblyBegin(XX, petsc_ierr)
! Computations can be done while messages are in transition
CALL VecAssemblyEnd(XX, petsc_ierr)

CALL VecGetOwnershipRange(XX, IVstart, IVend, petsc_ierr)

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 
    Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))") &
        "MPI rank : ",worker_rank_mpi,"| vector ownership for ","X", ":",IVstart," -- " , IVend
    CALL PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
    CALL PetscViewerASCIIOpen(COMM_MPI,"FX.output.1",PetscViewer, petsc_ierr);
    CALL PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
    CALL VecView(XX, PetscViewer, petsc_ierr)
    CALL PetscViewerDestroy(PetscViewer, petsc_ierr)
End If

!------------------------------------------------------------------------------
! At this point the right hand side vectors, filled with zeros and
! the solution vector filled with the dirichlet boundary values of
! load case 1 (LC1) are ready to be used.
! Now the right hand side vectors have to be modified with the
! prescribed displacements. This is currently done by multiplying A by
! X filled with the negative displacement values.
! (See above VecSetValues(XX, ... , -boundary_branch%leaves(2)%p_real8, ... )
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- compute_bndry_conditions "//TRIM(domain_char)
    CASE default
        timer_name = "compute_bndry_conditions"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 

!------------------------------------------------------------------------------ 
! Compute dirichlet boundary corrections of first right hand side vector
!------------------------------------------------------------------------------ 
IF (macro_order == 1) THEN
    CALL MatMult(AA,XX,FF_08(1), petsc_ierr);

ELSE IF (macro_order == 2) THEN
    CALL MatMult(AA,XX,FF_20(1), petsc_ierr);
    
END IF

IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
    ! Set zero values for dofs with prescribed displacements

    IF (macro_order == 1) THEN
        CALL VecSetValues(FF_08(1), boundary_branch%leaves(2)%dat_no, &
            gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)        
    ELSE IF (macro_order == 2) THEN
        CALL VecSetValues(FF_20(1), boundary_branch%leaves(2)%dat_no, &
            gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)            
    END IF
    
END IF 

IF (macro_order == 1) THEN
    CALL VecAssemblyBegin(FF_08(1), petsc_ierr)
    CALL VecAssemblyEnd(FF_08(1), petsc_ierr)  
ELSE IF (macro_order == 2) THEN
    CALL VecAssemblyBegin(FF_20(1), petsc_ierr)
    CALL VecAssemblyEnd(FF_20(1), petsc_ierr)         
END IF


!------------------------------------------------------------------------------
! Compute dirichlet boundary corrections of 2nd to 23rd
! right hand side vectors. The first on is already done and
! the 24th will be done afterwards when the columns and rows
! of A set to zero.
!------------------------------------------------------------------------------
Do ii = 2, no_lc-1_ik

    !------------------------------------------------------------------------------
    ! Get Bounds branch of LC ii
    ! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
    ! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
    !------------------------------------------------------------------------------
    write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",ii
    CALL search_branch(trim(desc), part_branch, boundary_branch, success, out_amount)

    IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
        ! Set prescribed displacements of LCii to solution vector
        CALL VecSetValues(XX, boundary_branch%leaves(2)%dat_no, &
            gnid_cref, -boundary_branch%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
    END IF 

    CALL VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition
    CALL VecAssemblyEnd(XX, petsc_ierr)

    !------------------------------------------------------------------------------ 
    ! Compute dirichlet boundary corrections of ii th right hand side vector.
    !------------------------------------------------------------------------------ 
    IF (macro_order == 1) THEN
        CALL MatMult(AA,XX,FF_08(ii), petsc_ierr);

        IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
            ! Set zero values for dofs with prescribed displacements
            CALL VecSetValues(FF_08(ii), boundary_branch%leaves(2)%dat_no, &
                gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
        END IF 
    
        CALL VecAssemblyBegin(FF_08(ii), petsc_ierr)
    ELSE IF (macro_order == 2) THEN
        CALL MatMult(AA,XX,FF_20(ii), petsc_ierr);

        IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
            ! Set zero values for dofs with prescribed displacements
            CALL VecSetValues(FF_20(ii), boundary_branch%leaves(2)%dat_no, &
                gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
        END IF 
    
        CALL VecAssemblyBegin(FF_20(ii), petsc_ierr)     
    END IF

End Do

!------------------------------------------------------------------------------
! Get Bounds branch of the last Load case
! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
!------------------------------------------------------------------------------
write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_", no_lc
CALL search_branch(trim(desc), part_branch, boundary_branch, success, out_amount)

!------------------------------------------------------------------------------
! Only set prescribed displacements if there are boundary nodes available.
! These are calculated in struct_preprocess subroutine generate_boundaries
!------------------------------------------------------------------------------
IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
    ! Set prescribed displacements of LC24 to solution vector
    CALL VecSetValues(XX, boundary_branch%leaves(2)%dat_no, &
        gnid_cref, -boundary_branch%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
END IF 

CALL VecAssemblyBegin(XX, petsc_ierr)
! Computations can be done while messages are in transition
CALL VecAssemblyEnd(XX, petsc_ierr)

!------------------------------------------------------------------------------
! Finalize the open assembleys
!------------------------------------------------------------------------------
IF (macro_order == 1) THEN
    Do ii = 2, no_lc-1
        CALL VecAssemblyEnd(FF_08(ii), petsc_ierr)
    End Do

ELSE IF (macro_order == 2) THEN
    Do ii = 2, no_lc-1
        CALL VecAssemblyEnd(FF_20(ii), petsc_ierr)
    End Do

END IF


!------------------------------------------------------------------------------
! Since we are filling XX with the prescribed displacements
! times -1. , we have to rescale XX before using it in
! MatZeroRowsColumns.
!------------------------------------------------------------------------------
CALL VecScale(XX, -1._rk, petsc_ierr)

!------------------------------------------------------------------------------
! Apply Dirichlet Boundaries to A and the 24th right hand side vector.
!------------------------------------------------------------------------------
IF (macro_order == 1) THEN
    CALL MatZeroRowsColumns(AA, boundary_branch%leaves(2)%dat_no, gnid_cref, &
        0.0_8, XX, FF_08(no_lc), petsc_ierr)
ELSE IF (macro_order == 2) THEN
    CALL MatZeroRowsColumns(AA, boundary_branch%leaves(2)%dat_no, gnid_cref, &
        0.0_8, XX, FF_20(no_lc), petsc_ierr)
END IF

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 

    CALL PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
    CALL PetscViewerASCIIOpen(COMM_MPI,"FF.output.1", PetscViewer, petsc_ierr);
    CALL PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)

    IF (macro_order == 1) THEN
        CALL VecView(FF_08( 1), PetscViewer, petsc_ierr)
    ELSE IF (macro_order == 2) THEN
        CALL VecView(FF_20( 1), PetscViewer, petsc_ierr)
    END IF

    CALL PetscViewerDestroy(PetscViewer, petsc_ierr)
End If

!------------------------------------------------------------------------------
! At this point the right hand side vectors are modified with the      
! prescribed displacements and the rows and columns of A representing  
! dofs with prescribed displacements are filled with zeros.            
! Now a linear solver context can be set up and the linear systems with
! constant operator A and variable right hand side can be solved       
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- solve_system "//TRIM(domain_char)
    CASE default
        timer_name = "solve_system"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 


!------------------------------------------------------------------------------
! Get and write another memory log.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

IF (worker_rank_mpi == 0) THEN
    collected_logs(6) = INT(time(), ik)
    collected_logs(13) = mem_global
    collected_logs(20) = status_global
END IF

!------------------------------------------------------------------------------
! Create linear solver context
!------------------------------------------------------------------------------
CALL KSPCreate(COMM_MPI, ksp, petsc_ierr)

!------------------------------------------------------------------------------
! Set operators. Here the matrix that defines the linear system
! also serves as the preconditioning matrix.
!------------------------------------------------------------------------------
CALL KSPSetOperators(ksp, AA, AA, petsc_ierr)

!------------------------------------------------------------------------------
! Set linear solver defaults for this problem (optional).
! - By extracting the KSP and PC contexts from the KSP context,
!   we can then directly CALL any KSP and PC routines to set
!   various options.
! - The following two statements are optional; all of these
!   parameters could alternatively be specified at runtime via
!   KSPSetFromOptions().  All of these defaults can be
!   overridden at runtime, as indicated below.
!------------------------------------------------------------------------------
CALL KSPSetTolerances(ksp,1.e-2/((m_size+1)*(m_size+1)),1.e-50_rk, &
        PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, petsc_ierr);

!------------------------------------------------------------------------------
! Set runtime options, e.g.,
!     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
! These options will override those specified above as long as
! KSPSetFromOptions() is CALLed _after_ any other customization
! routines.
!------------------------------------------------------------------------------
CALL KSPSetFromOptions(ksp, petsc_ierr)

If ((out_amount == "DEBUG") .AND. (worker_rank_mpi == 0)) THEN 
    filename=''
    write(filename,'(A,I0,A)')trim(job_dir)//trim(project_name)//"_",domain,"_usg.vtk"
    CALL write_data_head(filename, m_size/3)
End If

!------------------------------------------------------------------------------
! Solve the linear system
!------------------------------------------------------------------------------
! Add a results branch to the mesh branch on rank 0
! !! Initial try ... Communicate all results to master and !!
! !! do serial calc_effective_material_parameters          !!
!------------------------------------------------------------------------------
IF (worker_rank_mpi == 0) THEN

    CALL add_branch_to_branch(mesh_branch, esd_result_branch)
    ! Raise branch with 4 child branches
    ! 1. Displacements
    ! 2. Forces
    ! 3. Strains
    ! 4. Stresses
    CALL raise_branch("Results of domain "//domain_char , 4,  0, esd_result_branch)
    CALL raise_branch("Displacements"                   , 0,  0, esd_result_branch%branches(1))
    CALL raise_branch("Forces"                          , 0,  0, esd_result_branch%branches(2))
    CALL raise_branch("Strains"                         , 0,  0, esd_result_branch%branches(3))
    CALL raise_branch("Stresses"                        , 0,  0, esd_result_branch%branches(4))

    CALL log_tree(mesh_branch, un_lf, .FALSE.)
    
    !------------------------------------------------------------------------------
    ! Look again for the Part branch since the part_branch pointer 
    ! gets invalidated by dealloc of the branches array in add_branch_to_branch
    !------------------------------------------------------------------------------
    part_desc=''
    Write(part_desc,'(A,I0)')'Part_', parts
    CALL search_branch(trim(part_desc), mesh_branch, part_branch, success, out_amount)

    ! Allocate global displacement result
    Allocate(glob_displ(0:m_size-1))
    ! Allocate global forces result
    Allocate(glob_force(0:m_size-1))
    ! Allocate local bounds for global result
    Allocate(res_sizes(2,worker_size_mpi-1))
    
End If ! (worker_rank_mpi == 0) THEN

Do jj = 1, no_lc
    
    IF (macro_order == 1) THEN
        CALL KSPSolve(ksp, FF_08(jj), XX, petsc_ierr)
    ELSE IF (macro_order == 2) THEN
        CALL KSPSolve(ksp, FF_20(jj), XX, petsc_ierr)
    END IF 

    !------------------------------------------------------------------------------
    ! Get Bounds branch of LC jj
    ! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
    ! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
    !------------------------------------------------------------------------------
    write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",jj
    CALL search_branch(trim(desc), part_branch, boundary_branch, success, out_amount)
    
    !------------------------------------------------------------------------------
    ! Only set prescribed displacements if there are boundary nodes available.
    ! These are calculated in struct_preprocess subroutine generate_boundaries
    !------------------------------------------------------------------------------
    IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
        ! Extend XX with Boundary displacements
        CALL VecSetValues(XX, boundary_branch%leaves(2)%dat_no, &
        gnid_cref, boundary_branch%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
    END IF 
    
    CALL VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition
    CALL VecAssemblyEnd(XX, petsc_ierr)

    IF (macro_order == 1) THEN
    
        ! Calc reaction forces
        CALL MatMult(AA_org, XX, FF_08(jj), petsc_ierr);
    
        ! Get Pointer to result vector
        CALL VecGetArrayReadF90(XX,displ,petsc_ierr)
        
        ! Get Pointer to force vector
        CALL VecGetArrayReadF90(FF_08(jj),force,petsc_ierr)

    ELSE IF (macro_order == 2) THEN
    
    ! Calc reaction forces
        CALL MatMult(AA_org, XX, FF_20(jj), petsc_ierr);
    
        ! Get Pointer to result vector
        CALL VecGetArrayReadF90(XX,displ,petsc_ierr)
        
        ! Get Pointer to force vector
        CALL VecGetArrayReadF90(FF_20(jj),force,petsc_ierr)
    END IF 

    
    !------------------------------------------------------------------------------
    ! Master/Worker
    !------------------------------------------------------------------------------
    if (worker_rank_mpi > 0) then
        CALL mpi_send([IVstart,IVend], 2_mik, &
            MPI_Integer8,  0_mik, worker_rank_mpi, COMM_MPI, ierr)
        
        CALL mpi_send(displ, Int(IVend-IVstart,mik), &
            MPI_Real8,  0_mik, Int(worker_rank_mpi+worker_size_mpi,mik), COMM_MPI, ierr)

        CALL mpi_send(force, Int(IVend-IVstart,mik), &
            MPI_Real8,  0_mik, Int(worker_rank_mpi+2*worker_size_mpi,mik), COMM_MPI, ierr)

    Else ! Master

        !------------------------------------------------------------------------------
        ! Copy rank 0 local result
        !------------------------------------------------------------------------------
        glob_displ(IVstart:IVend-1) = displ
        glob_force(IVstart:IVend-1) = force
        
        !------------------------------------------------------------------------------
        ! Recv bounds of local results
        !------------------------------------------------------------------------------
        Do ii = 1, worker_size_mpi-1
            CALL mpi_recv(res_sizes(:,ii), 2_mik, MPI_Integer8,  Int(ii, mik), &
                Int(ii, mik), COMM_MPI, status_mpi, ierr)
        End Do
        
        !------------------------------------------------------------------------------
        ! Recv local parts of global result
        !------------------------------------------------------------------------------
        Do ii = 1, worker_size_mpi-1
            CALL mpi_recv(glob_displ(res_sizes(1,ii):res_sizes(2,ii)-1), &
                Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                MPI_Integer8,  Int(ii, mik), Int(ii+worker_size_mpi, mik), &
                COMM_MPI, status_mpi, ierr)

            CALL mpi_recv(glob_force(res_sizes(1,ii):res_sizes(2,ii)-1), &
                Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                MPI_Integer8,  Int(ii, mik), Int(ii+2*worker_size_mpi, mik) , &
                COMM_MPI, status_mpi, ierr)
        End Do

        !------------------------------------------------------------------------------
        ! Add leaf with displacements to the results branch
        !------------------------------------------------------------------------------
        write(desc,'(A)') "Displacements"
        CALL Add_Leaf_to_Branch(esd_result_branch%branches(1), trim(desc), m_size, glob_displ) 

        !------------------------------------------------------------------------------
        ! Add leaf with resultant forces to the results branch
        !------------------------------------------------------------------------------
        write(desc,'(A)') "Reaction Forces"
        CALL Add_Leaf_to_Branch(esd_result_branch%branches(2), trim(desc), m_size, glob_force) 


        If (out_amount == "DEBUG") THEN 
            write(desc,'(A,I2.2)') "DispRes", jj
            CALL write_vtk_data_Real8_vector_1D(matrix = reshape(glob_displ, &
                [3,size(glob_displ)/3]), filename = trim(filename),  &
                desc = trim(desc), head = .FALSE., location="POINT_DATA")

            write(desc,'(A,I2.2)') "ForcRes", jj
            CALL write_vtk_data_Real8_vector_1D(matrix = reshape(glob_force, &
                [3,size(glob_force)/3]), filename = trim(filename),  &
                desc = trim(desc), head = .FALSE., location="POINT_DATA")
        End if
    End If
End Do

!------------------------------------------------------------------------------
! Get and write another memory log.
!------------------------------------------------------------------------------
CALL mpi_system_mem_usage(COMM_MPI, mem_global, status_global) 

IF (worker_rank_mpi == 0) THEN
    collected_logs(7) = INT(time(), ik)
    collected_logs(14) = mem_global
    collected_logs(21) = status_global
END IF

!------------------------------------------------------------------------------
! All 24 linear system solutions are produced. 
! Effective stiffnesses can be calculated.
!------------------------------------------------------------------------------
if (worker_rank_mpi == 0) then

    CALL end_timer(TRIM(timer_name))

    Deallocate(glob_displ, res_sizes, glob_force)
    DEALLOCATE(nodes_in_mesh, zeros_R8)

    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- calc_eff_stiffness "//TRIM(domain_char)
    CASE default
        timer_name = "calc_eff_stiffness"
    End SELECT

    CALL start_timer(TRIM(timer_name), .FALSE.)
    CALL calc_effective_material_parameters(root, comm_nn, domain, &
        fh_mpi_worker, worker_size_mpi, collected_logs)
    CALL end_timer(TRIM(timer_name))
    
ELSE
    DEALLOCATE(part_branch)
End if
   
1000 CONTINUE

!------------------------------------------------------------------------------
! Remove matrices
!------------------------------------------------------------------------------
CALL KSPDestroy(ksp,    petsc_ierr)
CALL MatDestroy(AA,     petsc_ierr)
CALL MatDestroy(AA_org, petsc_ierr)
CALL VecDestroy(XX,     petsc_ierr)

Do ii = 1, no_lc

    IF (macro_order == 1) THEN
    
        CALL VecDestroy(FF_08(ii), petsc_ierr)

    ELSE IF (macro_order == 2) THEN
    
        CALL VecDestroy(FF_20(ii), petsc_ierr)
    END IF 

End Do

End Subroutine exec_single_domain

!------------------------------------------------------------------------------
! SUBROUTINE: Broadcast sequence to stop slave procs correctly
!------------------------------------------------------------------------------
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Stop a program.
!------------------------------------------------------------------------------  
Subroutine stop_slaves()

Integer(mik) :: ierr

CALL mpi_bcast(pro_path, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL mpi_bcast(pro_name, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
!------------------------------------------------------------------------------
! Bcast Serial_root_size = -1 ==> Signal for slave to stop
!------------------------------------------------------------------------------
CALL mpi_bcast(-1_ik, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

End Subroutine stop_slaves

!------------------------------------------------------------------------------
! SUBROUTINE: Broadcast sequence to stop slave procs correctly and print a mssg
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Stop a program. Tailored to the DTC/struct-process. Basically a wrapper.
!------------------------------------------------------------------------------  
SUBROUTINE print_err_stop_slaves(instring, tag)

CHARACTER(*), INTENT(IN), OPTIONAL :: tag

CHARACTER(SCL) :: tag_used, frmt, frmt_stp
CHARACTER(*), INTENT(IN) :: instring

tag_used = ""
IF (PRESENT(tag)) tag_used=tag

SELECT CASE (tag_used)
CASE ("warning")
    frmt     = FMT_WRN
    frmt_stp = FMT_TXT_STOP
CASE ("message")
    frmt     = FMT_MSG
    frmt_stp = FMT_TXT_STOP
CASE DEFAULT
    frmt     = FMT_ERR
    frmt_stp = FMT_ERR_STOP
END SELECT

WRITE(std_out, FMT_TXT_SEP)

! TODO: Repair this routine :-)
! CALL print_trimmed_text(std_out, TRIM(instring), frmt)
WRITE(std_out, frmt_stp)
WRITE(std_out, frmt) TRIM(instring)
WRITE(std_out, FMT_TXT_SEP)

!------------------------------------------------------------------------------  
! \TODO Receive confirmation from slaves, that they are ready to finish by 
! an MPI_BCAST. 
!------------------------------------------------------------------------------  
CALL stop_slaves()

END SUBROUTINE print_err_stop_slaves

End MODULE dtc_main_subroutines