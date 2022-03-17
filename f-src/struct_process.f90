!==============================================================================
!> \mainpage HLRS - Material structure process chain - Code Documentation
!> 
!> <hr>
!> \section desc Description
!>
!> A program collection organised in a tool-chain for the automatic 
!> determination of anisotropic elastic constants of microstructured 
!> materials from micro computer tomography.
!> The Flowchart of how it works can be seen here: \ref current.
!> <hr>
!> \section developers Developers 
!> 
!> Ralf Schneider, Xin Wang, Dmitry Khabi, Johannes Gebert
!>
!> FMPS is developed by Hans Wüstenberg
!>
!> for and with support from
!>
!> <b>High Performance Computing Center Stuttgart &mdash; HLRS</b><BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Johannes Gebert<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Uwe Küster<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dmitry Khabi<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Michael M. Resch
!>
!> <b>LASSO Ingenieurgesellschaft mbH</b><BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ulrich Hindenlang<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Artur Kurz<BR>
!>
!> <hr>
!> \section version Version
!>
!> 4.0 alpha
!>
!> <hr>
!> \section svn SVN 
!>
!> Versioning control is provided by FusionForge
!> https://gforge.hlrs.de/projects/struct-process
!>
!> <hr>
!> \section install Installation
!>
!> Installation is done by
!> <ol>
!>   <li> source setenv.sh <system>
!>   <li> make
!>
!> </ol>
!>
!> Where <system> can currently take the values "zeus, "julius" or "hawk".
!> Please not that in order to install the material structure process chain working copies of metis and petsc have to be installed on the executing system. Installation instructions for the dependencies can be found at this place: 
!> \ref prerequisites.
!>
!> <hr>
!> \section testing Testing
!>
!> Description of the provided regression and correctness checks can be found at this place: 
!> \ref regression.
!>
!> <hr>
!>
!==============================================================================

!==============================================================================
!> Module containing the execution of a single domain
!--
!------------------------------------------------------------------------------
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on: 03.01.2022
!------------------------------------------------------------------------------
Module sp_aux_routines

USE global_std
Use Operating_System
use puredat_com
use chain_routines
use linfe
use mpi
use gen_geometry
USE PETSC
Use calcmat

Implicit None
  
Contains

!============================================================================
!> Execution chain for the treatment of a single MVE
!>
!> \todo FMPS read .epp and .err from file is not that a good idea
Subroutine exec_single_domain(root, lin_domain, domain, job_dir, Active, fh_mpi_worker, &
    rank_mpi, size_mpi, comm_mpi)

TYPE(materialcard) :: bone
LOGICAL, PARAMETER :: DEBUG=.TRUE.
Integer(kind=mik), Intent(INOUT), Dimension(no_streams) :: fh_mpi_worker
Type(tBranch) :: domain_tree

Integer(kind=mik), intent(out) :: Active

Type(tBranch), Intent(inOut) :: root
Integer(kind=ik), Intent(in) :: domain, lin_domain
Character(LEN=*), Intent(in) :: job_dir
Integer(kind=mik), Intent(In) :: rank_mpi, size_mpi, comm_mpi


!----------------------------------------------------------------
Integer(kind=mik), Dimension(MPI_STATUS_SIZE) :: status_mpi

Type(tBranch), pointer :: boundary_branch, domain_branch, part_branch
Type(tBranch), pointer :: mesh_branch, meta_para, esd_result_branch

Integer(kind=mik)  :: ierr

Character(Len=mcl) :: desc, mesh_desc, filename
Character(Len=mcl) :: elt_micro

Character, Dimension(4*mcl) :: c_char_array

Integer(kind=c_int) :: stat_c_int
Integer             :: stat
Integer             :: umon

Character(len=mcl)  :: env_var

Integer(kind=ik)    :: m_size

Character(len=9)    :: domain_char

logical             :: success

Character(len=mcl)  :: timer_name, domain_desc, part_desc

Integer(kind=mik)         :: petsc_ierr
Type(tMat)                :: AA, AA_org
Type(tVec)                :: XX
Type(tVec), Dimension(24) :: FF
TYPE(tPETScViewer)        :: PetscViewer
Type(tKSP)                :: KSP
Integer(Kind=ik)          :: Istart,Iend, parts, IVstart, IVend
Integer(Kind=ik), Dimension(:), Allocatable :: nodes_in_mesh

Integer(kind=pd_ik), Dimension(:), Allocatable :: serial_pb
Integer(kind=pd_ik)                            :: serial_pb_size

Integer(Kind=pd_ik)  , Dimension(no_streams)  :: no_data, removed_data
Integer(kind=ik)                              :: domain_elems, ii, jj, kk, id
Integer(kind=ik), Dimension(:), Allocatable   :: gnid_cref
Integer(kind=ik), Dimension(:,:), Allocatable :: res_sizes

Real(KIND=rk), DIMENSION(:), Pointer       :: displ, force
Real(KIND=rk), DIMENSION(:), Allocatable   :: glob_displ, glob_force
Real(KIND=rk), DIMENSION(:), Allocatable   :: zeros_R8

Character, Dimension(:), Allocatable :: char_arr

CHARACTER(LEN=8) :: date
CHARACTER(LEN=10) :: time
CHARACTER(LEN=5)  :: timezone
CHARACTER(len=mcl) :: str

!--------------------------------------------------------------------------

Integer(Kind=ik), Dimension(60)    :: idxm_20, idxn_20
Real(kind=rk),    Dimension(60,60) :: K_loc_20

Integer(Kind=ik), Dimension(24)    :: idxm_08, idxn_08
Real(kind=rk),    Dimension(24,24) :: K_loc_08

integer(Kind=ik) :: preallo
        
!--------------------------------------------------------------------------

! Init worker_is_active status
Active = 0_mik

write(domain_char,'(I0)') domain

!--------------------------------------------------------------------------
! Get basic infos 
!--------------------------------------------------------------------------
CALL Search_branch("Input parameters", root, meta_para, success, DEBUG)

CALL pd_get(meta_para, "No of mesh parts per subdomain", parts)
CALL pd_get(meta_para, "Physical domain size" , bone%phdsize, 3)
CALL pd_get(meta_para, "Grid spacings" , bone%delta, 3)
CALL pd_get(meta_para, "Young_s modulus" , bone%E)
CALL pd_get(meta_para, "Poisson_s ratio" , bone%nu)
    
!------------------------------------------------------------------------------
! Rank = 0 -- Local master of comm_mpi
!------------------------------------------------------------------------------
If (rank_mpi == 0) then

    !------------------------------------------------------------------------------
    ! Create local tree
    !------------------------------------------------------------------------------
    ! CALL raise_tree(Trim(project_name), domain_tree)
    ! CALL include_branch_into_branch(s_b=meta_para, t_b=root, blind=.TRUE.)

    !------------------------------------------------------------------------------
    ! Add part branch to domain_tree
    !------------------------------------------------------------------------------
    ! part_desc=''
    ! Write(part_desc,'(A,I0)')'Part_', parts

    ! CALL add_branch_to_branch(domain_tree, domain_branch)

    ! CALL raise_branch(trim(part_desc), 0_pd_ik, 0_pd_ik, domain_branch)
    
    !------------------------------------------------------------------------------
    ! Create job directory in case of non-production run
    !------------------------------------------------------------------------------
    If (out_amount /= "PRODUCTION") then

        c_char_array(1:len(Trim(job_dir)//Char(0))) = str_to_char(Trim(job_dir)//Char(0))

        CALL Stat_Dir(c_char_array, stat_c_int)

        If(stat_c_int /= 0) Then

            CALL execute_command_line("mkdir -p "//trim(job_dir), CMDSTAT=stat)

            IF(stat /= 0) CALL print_err_stop(std_out, "Couldn't execute syscall mkpir -p "//TRIM(job_dir), 1)

            CALL Stat_Dir(c_char_array, stat_c_int)

            IF(stat_c_int /= 0) CALL print_err_stop(std_out, "Couldn't create directory "//TRIM(job_dir), 1)

        End If

        Write(un_lf,FMT_MSG_SEP)

    End If

    !------------------------------------------------------------------------------
    ! Write log and monitor file
    !------------------------------------------------------------------------------
    Write(un_lf,FMT_MSG_xAI0) "Domain No. : ", domain
    Write(un_lf,FMT_MSG)      "Job_dir    : "//Trim(job_dir)

    CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

    str = ''
    str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
    str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
    str = TRIM(str)//' '//timezone
    
    WRITE(un_lf, '(2A)') 'Start time: ', TRIM(str)

    CALL get_environment_Variable("HOSTNAME", env_var)
    Write(un_lf,FMT_MSG) "Host       : "//Trim(env_var)

    umon = give_new_unit()

    Open(unit=umon, file=Trim(job_dir)//Trim(project_name)//".mon",action="write", &
            status="replace")

    Write(umon,FMT_MSG_SEP)
    Write(umon,FMT_MSG_xAI0) "Domain No. : ", domain

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

    CALL generate_geometry(root, lin_domain, domain, job_dir, fh_mpi_worker, success)

    if (.not.success) then
        write(*,FMT_WRN)"generate_geometry() returned .FALSE."
    End if

    CALL end_timer(trim(timer_name))
    
    CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

    str = ''
    str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
    str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
    str = TRIM(str)//' '//timezone
    
    WRITE(un_lf, '(2A)') 'End time: ', TRIM(str)

    close(umon)

    !------------------------------------------------------------------------------
    ! Look for the Domain branch
    !------------------------------------------------------------------------------
    domain_desc=''
    Write(domain_desc,'(A,I0)')'Domain ', domain
    
    CALL search_branch(trim(domain_desc), root, domain_branch, success, DEBUG)

    !------------------------------------------------------------------------------
    ! Get the no of nodes per part
    !------------------------------------------------------------------------------
    mesh_desc = ''
    Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_', domain
    
    CALL search_branch(trim(mesh_desc), domain_branch, mesh_branch, success, DEBUG)
    CALL pd_get(mesh_branch, 'No of nodes in mesh',  nodes_in_mesh)

    !------------------------------------------------------------------------------
    ! Set the global matrix size
    !------------------------------------------------------------------------------
    m_size = nodes_in_mesh(1) * 3

    Do ii = 1, parts-1

        !------------------------------------------------------------------------------
        ! Look for the Part branch
        !------------------------------------------------------------------------------
        domain_desc=''
        Write(part_desc,'(A,I0)')'Part_',ii
        CALL search_branch(trim(part_desc), domain_branch, part_branch, success)

        If (.NOT. success) Then
            WRITE(mssg,'(A,I0,A,L,A,I0,A)') "Something bad and unexpected happend &
            &in exec_single_domain! Looking for branch of part ",ii," returned ", &
                success, "MPI proc ",rank_mpi," halted."
            CALL print_err_stop(std_out, mssg, 1)
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
        
    CALL search_branch(trim(part_desc), domain_branch, part_branch, success, DEBUG)

    !------------------------------------------------------------------------------
    ! Broadcast matrix size. 
    ! TODO could also be included into part branches.
    !------------------------------------------------------------------------------
    CALL mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)

!------------------------------------------------------------------------------
! Ranks > 0 - Workers
!------------------------------------------------------------------------------
Else

    CALL mpi_recv(serial_pb_size, 1_mik, MPI_INTEGER8, 0_mik, &
        rank_mpi, COMM_MPI, status_mpi, ierr)

    if (allocated(serial_pb)) deallocate(serial_pb)
    
    Allocate(serial_pb(serial_pb_size))

    CALL mpi_recv(serial_pb, INT(serial_pb_size,mik), MPI_INTEGER8, &
        0_mik, rank_mpi, COMM_MPI, status_mpi, ierr)

    ! Deserialize part branch
    CALL Start_Timer("Deserialize part branch branch")

    Allocate(part_branch)
    
    CALL deserialize_branch(part_branch, serial_pb, .TRUE.)

    CALL End_Timer("Deserialize part branch branch")

    CALL mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)
            
End If ! (rank_mpi == 0) then

!------------------------------------------------------------------------------
! Setup the linear System with a constant system matrix A. Once that
! is done setup the multiple right hand sides and solve the linear
! system multiple times. Save the solutions to calculate effective
! stiffness matirces.
!------------------------------------------------------------------------------
IF (rank_mpi == 0) THEN   ! Sub Comm Master

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
! Calculate amount of memory to allocate.
!------------------------------------------------------------------------------
preallo = (part_branch%leaves(5)%dat_no * 3) / parts + 1

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

CALL MatSeqAIJSetPreallocation(AA, preallo, PETSC_NULL_INTEGER, petsc_ierr)
CALL MatMPIAIJSetPreallocation(AA, preallo, PETSC_NULL_INTEGER, preallo, PETSC_NULL_INTEGER, petsc_ierr)

CALL MatSeqAIJSetPreallocation(AA_org, preallo, PETSC_NULL_INTEGER, petsc_ierr)
CALL MatMPIAIJSetPreallocation(AA_org, preallo, PETSC_NULL_INTEGER, preallo, PETSC_NULL_INTEGER, petsc_ierr)

CALL MatSetOption(AA    ,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,petsc_ierr)
CALL MatSetOption(AA_org,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,petsc_ierr)

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 
    CALL MatGetOwnershipRange(AA, Istart, Iend, petsc_ierr)
    Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
        "MPI rank : ",rank_mpi,"| matrix ownership for ","A    ", ":",Istart," -- " , Iend

    CALL MatGetOwnershipRange(AA_org, Istart, Iend, petsc_ierr)

    Write(un_lf, "('MM ', A,I4,A,A6,2(A,I9))") &
    "MPI rank : ", rank_mpi, "| matrix ownership for ", "A_org", ":",Istart, " -- " , Iend
End If

!------------------------------------------------------------------------------
! part_branch%leaves(1) : Local  node ids   in part branch
! part_branch%leaves(2) : Coordinate values in part branch
! part_branch%leaves(3) : Global node ids   in part branch
! part_branch%leaves(4) : NOT USED ! (For Element Numbers)
! part_branch%leaves(5) : Topology
! No Cross reference global nid to local nid is currently needed since 
! renumbering is deactivated in mod_mesh_partitioning.f90 : part_mesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Assemble matrix
!------------------------------------------------------------------------------
CALL pd_get(root%branches(1), 'Element type  on micro scale', char_arr)
elt_micro = char_to_str(char_arr)

if (TRIM(elt_micro) == "HEX08") then

    K_loc_08 = Hexe08(bone)

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
        
        idxn_08 = idxm_08

        CALL MatSetValues(AA, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)

        CALL MatSetValues(AA_org, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)
        
    End Do

else if (elt_micro == "HEX20") then
    
    K_loc_20 = Hexe20()

    domain_elems = part_branch%leaves(5)%dat_no / 20

    Do ii = 1, domain_elems

        kk = 1
        
        ! Translate Topology to global dofs *
        Do jj = (ii-1)*20+1 , ii*20
            
            id = part_branch%leaves(5)%p_int8(jj)
            id = (id - 1) * 3
            
            idxm_20(kk)   = (id    ) !*cr
            idxm_20(kk+1) = (id + 1) !*cr
            idxm_20(kk+2) = (id + 2) !*cr
            
            kk = kk + 3
            
        End Do
        
        idxn_20 = idxm_20

        CALL MatSetValues(AA, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
        CALL MatSetValues(AA_org, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
        
    End Do

end if


IF (rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- MatAssemblyBegin "//TRIM(domain_char)
    CASE default
        timer_name = "MatAssemblyBegin"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 

CALL MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
CALL MatAssemblyBegin(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)
! Computations can be done while messages are in transition
CALL MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
CALL MatAssemblyEnd(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

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
IF (rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- create_rh_solution_vetors "//TRIM(domain_char)
    CASE default
        timer_name = "create_rh_solution_vetors"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
END IF 

Do ii = 1, 24

    !------------------------------------------------------------------------------
    ! Create load vectors
    !------------------------------------------------------------------------------
    CALL VecCreate(COMM_MPI, FF(ii), petsc_ierr)
    CALL VecSetSizes(FF(ii), PETSC_DECIDE, m_size, petsc_ierr)
    CALL VecSetFromOptions(FF(ii), petsc_ierr)
    CALL VecSet(FF(ii), 0._rk,petsc_ierr)
    
    CALL VecGetOwnershipRange(FF(ii), IVstart, IVend, petsc_ierr)

    If (out_amount == "DEBUG") THEN 
        Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
            "MPI rank : ",rank_mpi,"| vector ownership for ","F", ":",IVstart," -- " , IVend
    End If
    
    !------------------------------------------------------------------------------
    CALL VecAssemblyBegin(FF(ii), petsc_ierr)
    ! Computations can be done while messages are transferring
End Do
Do ii = 1, 24 
    CALL VecAssemblyEnd(FF(ii), petsc_ierr)
End Do
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Get Bounds branch of LC 1
! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
!------------------------------------------------------------------------------ 
write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",1
CALL search_branch(trim(desc), part_branch, boundary_branch, success, DEBUG)

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
IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 
    Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))") &
        "MPI rank : ",rank_mpi,"| vector ownership for ","X", ":",IVstart," -- " , IVend
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
IF (rank_mpi == 0) THEN   ! Sub Comm Master
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
CALL MatMult(AA,XX,FF(1), petsc_ierr);

IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
    ! Set zero values for dofs with prescribed displacements
    CALL VecSetValues(FF(1), boundary_branch%leaves(2)%dat_no, &
        gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
END IF 

CALL VecAssemblyBegin(FF(1), petsc_ierr)
CALL VecAssemblyEnd(FF(1), petsc_ierr)

!------------------------------------------------------------------------------
! Compute dirichlet boundary corrections of 2nd to 23rd
! right hand side vectors. The first on is already done and
! the 24th will be done afterwards when the columns and rows
! of A set to zero.
!------------------------------------------------------------------------------
Do ii = 2, 23

    !------------------------------------------------------------------------------
    ! Get Bounds branch of LC ii
    ! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
    ! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
    !------------------------------------------------------------------------------
    write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",ii
    CALL search_branch(trim(desc), part_branch, boundary_branch, success, DEBUG)

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
    CALL MatMult(AA,XX,FF(ii), petsc_ierr);

    IF(boundary_branch%leaves(2)%dat_no /= 0_ik) THEN
        ! Set zero values for dofs with prescribed displacements
        CALL VecSetValues(FF(ii), boundary_branch%leaves(2)%dat_no, &
            gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
    END IF 

    CALL VecAssemblyBegin(FF(ii), petsc_ierr)
    
End Do

!------------------------------------------------------------------------------
! Get Bounds branch of LC 24
! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
!------------------------------------------------------------------------------
write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",24
CALL search_branch(trim(desc), part_branch, boundary_branch, success, DEBUG)

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
Do ii = 2, 23
    CALL VecAssemblyEnd(FF(ii), petsc_ierr)
End Do

!------------------------------------------------------------------------------
! Since we are filling XX with the prescribed displacements
! times -1. , we have to rescale XX before using it in
! MatZeroRowsColumns.
!------------------------------------------------------------------------------
CALL VecScale(XX, -1._rk, petsc_ierr)

!------------------------------------------------------------------------------
! Apply Dirichlet Boundaries to A and the 24th right hand side vector.
!------------------------------------------------------------------------------
CALL MatZeroRowsColumns(AA, boundary_branch%leaves(2)%dat_no, gnid_cref, &
        0.0_8, XX, FF(24), petsc_ierr)

!------------------------------------------------------------------------------
! End timer
!------------------------------------------------------------------------------
IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

If (out_amount == "DEBUG") THEN 
    CALL PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
    CALL PetscViewerASCIIOpen(COMM_MPI,"FF.output.1", PetscViewer, petsc_ierr);
    CALL PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
    CALL VecView(FF( 1), PetscViewer, petsc_ierr)
    CALL PetscViewerDestroy(PetscViewer, petsc_ierr)
End If

!------------------------------------------------------------------------------
! At this point the right hand side vectors are modified with the      
! prescribed displacements and the rows and columns of A representing  
! dofs with prescribed displacements are filled with zeros.            
! Now a linear solver context can be set up and the linear systems with
! constant operator A and variable right hand side can be solved       
!------------------------------------------------------------------------------
IF (rank_mpi == 0) THEN   ! Sub Comm Master
    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- solve_system "//TRIM(domain_char)
    CASE default
        timer_name = "solve_system"
    End SELECT
    
    CALL start_timer(TRIM(timer_name), .FALSE.)
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

If ((out_amount == "DEBUG") .AND. (rank_mpi == 0)) THEN 
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
IF (rank_mpi == 0) THEN

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
    CALL search_branch(trim(part_desc), mesh_branch, part_branch, success, DEBUG)

    ! Allocate global displacement result
    Allocate(glob_displ(0:m_size-1))
    ! Allocate global forces result
    Allocate(glob_force(0:m_size-1))
    ! Allocate local bounds for global result
    Allocate(res_sizes(2,size_mpi-1))
    
End If

Do jj = 1,24
    
    CALL KSPSolve(ksp, FF(jj), XX, petsc_ierr)
    
    !------------------------------------------------------------------------------
    ! Get Bounds branch of LC jj
    ! boundary_branch%leaves(1)%p_int8  : Boundary displacement node global ids
    ! boundary_branch%leaves(2)%p_real8 : Boundary displacement values
    !------------------------------------------------------------------------------
    write(desc,'(A,I0)') "Boundaries_"//trim(domain_char)//"_",jj
    CALL search_branch(trim(desc), part_branch, boundary_branch, success, DEBUG)
    
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

    ! Calc reaction forces
    CALL MatMult(AA_org, XX, FF(jj), petsc_ierr);
    
    ! Get Pointer to result vector
    CALL VecGetArrayReadF90(XX,displ,petsc_ierr)
    
    ! Get Pointer to force vector
    CALL VecGetArrayReadF90(FF(jj),force,petsc_ierr)
    
    !------------------------------------------------------------------------------
    ! Master/Worker
    !------------------------------------------------------------------------------
    if (rank_mpi > 0) then
        CALL mpi_send([IVstart,IVend], 2_mik, &
            MPI_Integer8,  0_mik, rank_mpi, COMM_MPI, ierr)
        
        CALL mpi_send(displ, Int(IVend-IVstart,mik), &
            MPI_Real8,  0_mik, Int(rank_mpi+size_mpi,mik), COMM_MPI, ierr)

        CALL mpi_send(force, Int(IVend-IVstart,mik), &
            MPI_Real8,  0_mik, Int(rank_mpi+2*size_mpi,mik), COMM_MPI, ierr)

    Else ! Master

        !------------------------------------------------------------------------------
        ! Copy rank 0 local result
        !------------------------------------------------------------------------------
        glob_displ(IVstart:IVend-1) = displ
        glob_force(IVstart:IVend-1) = force
        
        !------------------------------------------------------------------------------
        ! Recv bounds of local results
        !------------------------------------------------------------------------------
        Do ii = 1, size_mpi-1
            CALL mpi_recv(res_sizes(:,ii), 2_mik, MPI_Integer8,  Int(ii, mik), &
                Int(ii, mik), COMM_MPI, status_mpi, ierr)
        End Do
        
        !------------------------------------------------------------------------------
        ! Recv local parts of global result
        !------------------------------------------------------------------------------
        Do ii = 1, size_mpi-1
            CALL mpi_recv(glob_displ(res_sizes(1,ii):res_sizes(2,ii)-1), &
                Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                MPI_Integer8,  Int(ii, mik), Int(ii+size_mpi, mik), &
                COMM_MPI, status_mpi, ierr)

            CALL mpi_recv(glob_force(res_sizes(1,ii):res_sizes(2,ii)-1), &
                Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                MPI_Integer8,  Int(ii, mik), Int(ii+2*size_mpi, mik) , &
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
! All 24 linear system solutions are produced. 
! Effective stiffnesses can be calculated.
!------------------------------------------------------------------------------
if (rank_mpi == 0) then

    CALL end_timer(TRIM(timer_name))

    Deallocate(glob_displ, res_sizes, glob_force)

    SELECT CASE (timer_level)
    CASE (1)
        timer_name = "+-- calc_eff_stiffness "//TRIM(domain_char)
    CASE default
        timer_name = "calc_eff_stiffness"
    End SELECT

    CALL start_timer(trim(timer_name), .FALSE.)
    CALL calc_effective_material_parameters(root, lin_domain, domain, fh_mpi_worker)
    CALL end_timer(trim(timer_name))

End if

!------------------------------------------------------------------------------
! Remove matrices
!------------------------------------------------------------------------------
CALL MatDestroy(AA,     petsc_ierr)
CALL MatDestroy(AA_org, petsc_ierr)
CALL VecDestroy(XX,     petsc_ierr)

Do ii = 1, 24
    CALL VecDestroy(FF(ii), petsc_ierr)
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
SUBROUTINE print_err_stop_slaves(text, tag)

CHARACTER(LEN=*) , INTENT(IN) :: text
CHARACTER(LEN=*) , INTENT(IN), OPTIONAL :: tag

CHARACTER(LEN=SCL) :: tag_used, frmt, frmt_stp
Integer(mik) :: ierr

tag_used = ""
IF (PRESENT(tag)) tag_used=tag

SELECT CASE (tag_used)
CASE ("warning")
    frmt     = FMT_WRN
    frmt_stp = FMT_WRN_STOP
CASE ("message")
    frmt     = FMT_MSG
    frmt_stp = FMT_MSG_STOP
CASE DEFAULT
    frmt     = FMT_ERR
    frmt_stp = FMT_ERR_STOP
END SELECT

WRITE(std_out, frmt) TRIM(text)
WRITE(std_out, frmt_stp)

CALL stop_slaves()

END SUBROUTINE print_err_stop_slaves

End Module sp_aux_routines

!------------------------------------------------------------------------------
!> Struct Process main programm
!------------------------------------------------------------------------------
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on: 03.01.2022
!>
!> \todo At >> Recieve Job_Dir << calculate the Job_Dir from the Domain number
!>       and the Domain Decomposition parameters.
!> \todo During geometry generation do not allocate the Topology field for the
!>       max possible number of elements. Instead check with count for actual
!>       number of elements
!------------------------------------------------------------------------------
Program main_struct_process

USE iso_c_binding
USE global_std
USE user_interaction
USE formatted_plain
USE puredat 
USE meta
USE meta_puredat_interface
USE auxiliaries
USE chain_routines
USE MPI
USE decomp 
USE sp_aux_routines
USE PETSC
USE petsc_opt

Implicit None
! Parameter
INTEGER(KIND=ik), PARAMETER :: debug = 2   ! Choose an even integer!!

! Always provide in/out for meta driven environments
TYPE(materialcard) :: bone

!-- MPI Variables
INTEGER(KIND=mik) :: ierr, rank_mpi, size_mpi
INTEGER(KIND=mik) :: petsc_ierr
INTEGER(KIND=mik) :: worker_rank_mpi, worker_size_mpi
INTEGER(KIND=mik) :: Active, request, finished, worker_comm, WRITE_ROOT_COMM

INTEGER(KIND=mik), Dimension(no_streams)       :: fh_mpi_root, fh_mpi_worker
INTEGER(KIND=mik), Dimension(MPI_STATUS_SIZE)  :: status_mpi
INTEGER(KIND=mik), Dimension(:,:), Allocatable :: statuses_mpi
INTEGER(KIND=mik), Dimension(:)  , Allocatable :: worker_is_active, req_list
INTEGER(KIND=mik) :: mii, mjj

CHARACTER(Len=mcl)  :: link_name = 'struct process'

INTEGER(kind=c_int) :: stat_c_int

Type(tBranch)       :: groot, root, phi_tree
Type(tBranch), pointer :: ddc, meta_para, result_branch, mesh_branch

CHARACTER           , DIMENSION(4*mcl)          :: c_char_array
CHARACTER           , DIMENSION(:), ALLOCATABLE :: char_arr
CHARACTER(LEN=4*mcl), DIMENSION(:), ALLOCATABLE :: domain_path
CHARACTER(LEN=mcl)  , DIMENSION(:), ALLOCATABLE :: m_rry      

CHARACTER(LEN=4*mcl) :: job_dir
CHARACTER(LEN=mcl)   :: cmd_arg_history='', tmp_fn=''
CHARACTER(LEN=mcl)   :: muCT_pd_path, muCT_pd_name, domain_desc, mesh_desc, binary
CHARACTER(LEN=8)     :: elt_micro, output
CHARACTER(LEN=1)     :: restart='N', restart_cmd_arg='U' ! U = 'undefined'

REAL(KIND=rk), DIMENSION(3) :: delta
REAL(KIND=rk) :: strain

INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE :: Domains, workers_assigned_domains, nn_D
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE :: Domain_stats

INTEGER(KIND=ik), DIMENSION(3) :: xa_d, xe_d, vdim

INTEGER(KIND=ik) :: nn, ii, jj, kk, dc, stat, computed_domains = 0
INTEGER(KIND=ik) :: No_of_domains, path_count, tmp_un
INTEGER(KIND=ik) :: alloc_stat, aun, free_file_handle
INTEGER(KIND=ik) :: Domain, llimit, parts, elo_macro

INTEGER(KIND=pd_ik), DIMENSION(:), ALLOCATABLE :: serial_root
INTEGER(KIND=pd_ik), DIMENSION(no_streams) :: dsize, removed_data

INTEGER(KIND=pd_ik) :: serial_root_size

LOGICAL :: success, stat_exists, heaxist, stp = .FALSE., create_new_header = .FALSE.

!----------------------------------------------------------------------------
CALL mpi_init(ierr)
CALL print_err_stop(std_out, "MPI_INIT didn't succeed", ierr)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't be retrieved", ierr)

!------------------------------------------------------------------------------
! Rank 0 -- Init (Master) Process and broadcast init parameters 
!------------------------------------------------------------------------------
If (rank_mpi == 0) THEN

    If (size_mpi < 2) THEN
        mssg = "We need at least 2 MPI processes to execute this program."
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    CALL Start_Timer("Init Process")

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, stp, restart_cmd_arg, cmd_arg_history)
    IF(stp) GOTO 1000
    
    IF (in%full=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        mssg = "No input file given"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'dtc' 
    global_meta_program_keyword = 'TENSOR_COMPUTATION'
    CALL meta_append(m_rry)
        
    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title(["Dr.-Ing. Ralf Schneider (HLRS, NUM)", "Johannes Gebert, M.Sc. (HLRS, NUM) "])

    IF(debug >=0) WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/temporary.std_out"
    WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          

    !------------------------------------------------------------------------------
    ! Set input paths
    !------------------------------------------------------------------------------
    ! It is strongly recommended not to play around with these paths carelessly.
    ! Some of the dependencies are easily overlooked.
    !------------------------------------------------------------------------------
    muCT_pd_path = TRIM(in%path)
    muCT_pd_name = TRIM(in%bsnm)

    !------------------------------------------------------------------------------
    ! Output directory and 
    ! Implicitly creates a subdirectory.
    !------------------------------------------------------------------------------    
    outpath = TRIM(out%path)//TRIM(out%bsnm)//"/"
    project_name = "results" ! TRIM(out%bsnm)//

    pro_path = outpath
    pro_name = project_name

    CALL meta_read('MICRO_ELMNT_TYPE' , m_rry, elt_micro)
    CALL meta_read('OUT_FMT'          , m_rry, output)
    CALL meta_read('RESTART'          , m_rry, restart)
write(*,*) "Restart: ", restart
    CALL meta_read('SIZE_DOMAIN'      , m_rry, bone%phdsize)
    CALL meta_read('SPACING'          , m_rry, bone%delta)
    CALL meta_read('DIMENSIONS'       , m_rry, vdim)
    CALL meta_read('LO_BNDS_DMN_RANGE', m_rry, xa_d)
    CALL meta_read('UP_BNDS_DMN_RANGE', m_rry, xe_d)
    CALL meta_read('BINARIZE_LO'      , m_rry, llimit)
    CALL meta_read('MESH_PER_SUB_DMN' , m_rry, parts)
    CALL meta_read('RVE_STRAIN'       , m_rry, strain)
    CALL meta_read('YOUNG_MODULUS'    , m_rry, bone%E)
    CALL meta_read('POISSON_RATIO'    , m_rry, bone%nu)
    CALL meta_read('MACRO_ELMNT_ORDER', m_rry, elo_macro)

    !------------------------------------------------------------------------------
    ! Restart handling
    !------------------------------------------------------------------------------
    ! Standard lock file of MeRaDat not used, because it would interfere with the 
    ! status file of the struct process. 
    !------------------------------------------------------------------------------
    ! CALL meta_handle_lock_file(restart, restart_cmd_arg)
    CALL meta_compare_restart(restart, restart_cmd_arg)

    !------------------------------------------------------------------------------
    ! Spawn a results file
    !------------------------------------------------------------------------------
    ! This log file may collide with the original log file (!)
    ! The regular struct_process log file contains still has the "old" basename!
    !------------------------------------------------------------------------------
    CALL meta_start_ascii(fh_mon, mon_suf)

    If (out_amount == "DEBUG") CALL meta_start_ascii(fh_log, log_suf)
    
    IF (std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL meta_write('DBG_LVL', out_amount )

    !------------------------------------------------------------------------------
    ! Warning / Error handling
    !------------------------------------------------------------------------------
    IF ( (bone%phdsize(1) /= bone%phdsize(2)) .OR. (bone%phdsize(1) /= bone%phdsize(3)) ) THEN
        mssg = 'Currently, all 3 dimensions of the physical domain size must be equal!'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF
    
    IF ( (delta(1) /= delta(2)) .OR. (delta(1) /= delta(3)) ) THEN
        mssg = 'Currently, the spacings of all 3 dimensions must be equal!'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    IF ( (xa_d(1) > xe_d(1)) .OR. (xa_d(2) > xe_d(2)) .or. (xa_d(3) > xe_d(3)) ) THEN
        mssg = 'Input parameter error: Start value of domain range larger than end value.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    !------------------------------------------------------------------------------
    ! Program breaks if the phdsize is not taking the boundary nodes into account (!).
    ! Therefore, the boundaries are calculated with + 2 Voxels
    !------------------------------------------------------------------------------
    IF ((bone%phdsize(1) > (vdim(1) + 2_ik)*bone%delta(1)) .OR. & 
        (bone%phdsize(2) > (vdim(2) + 2_ik)*bone%delta(2)) .OR. & 
        (bone%phdsize(3) > (vdim(3) + 2_ik)*bone%delta(3))) THEN
        mssg = 'The domains are larger than the field of view.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF
    
    !------------------------------------------------------------------------------
    ! Each subdomain gets computed by a user defined amount of processors. This 
    ! amount of processors equals to a specific amount of mesh parts.
    ! Ideally, all processors are used. Therefore, MOD(size_mpi-1, parts) shall 
    ! resolve without a remainder. "-1" to take the master process into account.
    !------------------------------------------------------------------------------
    ! This dependency is checked after the restart handling again.
    !------------------------------------------------------------------------------
    IF (MOD(size_mpi-1, parts) /= 0) THEN
        mssg = 'mod(size_mpi-1,parts) /= 0 ! This case is not supported.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 
    
    !------------------------------------------------------------------------------
    ! Raise and build meta_para tree
    ! Hardcoded, implicitly given order of the leafs. 
    ! DO NOT CHANGE ORDER WITHOUT MODIFYING ALL OTHER INDICES REGARDING »meta_para«
    !------------------------------------------------------------------------------
    Allocate(meta_para)
    CALL raise_tree("Input parameters", meta_para)

    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_path"                , mcl , str_to_char(muCT_pd_path)) 
    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_name"                , mcl , str_to_char(muCT_pd_name)) 
    CALL add_leaf_to_branch(meta_para, "Physical domain size"                 , 3_ik, bone%phdsize) 
    CALL add_leaf_to_branch(meta_para, "Lower bounds of selected domain range", 3_ik, xa_d) 
    CALL add_leaf_to_branch(meta_para, "Upper bounds of selected domain range", 3_ik, xe_d)      
    
    CALL add_leaf_to_branch(meta_para, "Grid spacings"                        , 3_rk, bone%delta) 
    CALL add_leaf_to_branch(meta_para, "Lower limit of iso value"      , 1_ik, [llimit])      
    CALL add_leaf_to_branch(meta_para, "Element type  on micro scale"  , len(elt_micro) , str_to_char(elt_micro))      
    CALL add_leaf_to_branch(meta_para, "No of mesh parts per subdomain", 1_ik           , [parts]) 
    CALL add_leaf_to_branch(meta_para, "Output Format"                 , len(output)    , str_to_char(output)) 
    
    CALL add_leaf_to_branch(meta_para, "Average strain on RVE"         , 1_ik           , [strain])    
    CALL add_leaf_to_branch(meta_para, "Young_s modulus"               , 1_ik           , [bone%E]) 
    CALL add_leaf_to_branch(meta_para, "Poisson_s ratio"               , 1_ik           , [bone%nu]) 
    CALL add_leaf_to_branch(meta_para, "Element order on macro scale"  , 1_ik           , [elo_macro]) 
    CALL add_leaf_to_branch(meta_para, "Output amount"                 , len(out_amount), str_to_char(out_amount)) 
    
    CALL add_leaf_to_branch(meta_para, "Restart"                       , 1_ik           , str_to_char(restart)) 
    CALL add_leaf_to_branch(meta_para, "Number of voxels per direction", 3_ik           , vdim) 

    !------------------------------------------------------------------------------
    ! Prepare output directory via CALLing the c function.
    ! Required, because INQUIRE only acts on files, not on directories.
    ! File exists if stat_c_int = 0 
    !------------------------------------------------------------------------------
    c_char_array(1:LEN(TRIM(outpath)//CHAR(0))) = str_to_char(TRIM(outpath)//CHAR(0))
    CALL Stat_Dir(c_char_array, stat_c_int)

    IF(stat_c_int /= 0) THEN

        CALL execute_command_line("mkdir -p "//TRIM(outpath),CMDSTAT=stat)

        IF(stat /= 0) THEN
            mssg = 'Could not execute syscall »mkdir -p '//trim(outpath)//'«.'
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF
        
        CALL Stat_Dir(c_char_array, stat_c_int)

        IF(stat_c_int /= 0) THEN
            mssg = 'Could not create the output directory »'//TRIM(outpath)//'«.'
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF
    ELSE 
        WRITE(fh_mon, FMT_MSG) "Reusing the output directory"
        WRITE(fh_mon, FMT_MSG) TRIM(outpath)
    END IF

    CALL link_start(link_name, .TRUE., .FALSE., success)
    IF (.NOT. success) THEN
        mssg = "Something went wrong during link_start"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF
        
    !------------------------------------------------------------------------------
    ! Allocate and init field for selected domain range
    ! Required for the worker_is_active tracking
    ! Required for steering the allocation of memory for writing data into 
    ! the stream files (raise leaves etc.)
    !------------------------------------------------------------------------------
    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

    Allocate(Domains(No_of_domains), stat=alloc_stat)
    CALL alloc_err("Domains", alloc_stat)
write(*,*) "No of domains", No_of_domains
    Allocate(Domain_stats(No_of_domains), stat=alloc_stat)
    CALL alloc_err("Domain_stats", alloc_stat)
    Domain_stats = -9999999_ik

    Allocate(domain_path(0:No_of_domains))
    domain_path = ''

    !------------------------------------------------------------------------------
    ! Check whether there already is a project header and an worker_is_active tracker
    ! Header-files are PureDat, not MeRaDat. Therefore via pro_path/pro_name
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM( in%p_n_bsnm)//".head", EXIST=heaxist)
    INQUIRE(FILE = TRIM(out%p_n_bsnm)//".status", EXIST=stat_exists)
    aun = give_new_unit()

    !------------------------------------------------------------------------------
    ! Check if input header exists. If not --> create one.
    !------------------------------------------------------------------------------
    IF(.NOT. heaxist) THEN
        free_file_handle = give_new_unit()
        CALL convert_meta_to_puredat(free_file_handle, m_rry)
    END IF

    !------------------------------------------------------------------------------
    ! Check if output header exists. If not --> Internally reset to restart='N'.
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM(pro_path)//TRIM(pro_name)//".head", EXIST=heaxist)

    IF((restart == 'Y') .AND. (.NOT. heaxist)) create_new_header=.TRUE.



write(*,*) "pro_path: ", TRIM(pro_path)
write(*,*) "pro_name: ", TRIM(pro_name)

write(*,*) "Restart: ", restart
    !------------------------------------------------------------------------------
    ! The Output name normally is different than the input name.
    ! Therefore, an existing header implies a restart.
    !------------------------------------------------------------------------------
    IF ((restart == 'N') .OR. (create_new_header)) THEN ! BasiCALLy a completely new computation

        !------------------------------------------------------------------------------
        ! project_name --> out%p_n_bsnm/bsnm --> subdirectory with file name = bsnm.suf
        !------------------------------------------------------------------------------
        IF ((stat_exists) .AND. (.NOT. create_new_header)) THEN
            mssg = "No restart requested, but the file '"&
                //TRIM(out%p_n_bsnm)//".status' already exists."
            CALL print_err_stop_slaves(mssg, "warning"); GOTO 1000
        END IF 

write(*,*) "Restart: dafuq!"
        !------------------------------------------------------------------------------
        ! project_name --> out%p_n_bsnm/bsnm --> subdirectory with file name = bsnm.suf
        !------------------------------------------------------------------------------
        ! Tree which is fed back = root. A collection of stream paths and Null pointers
        !------------------------------------------------------------------------------
        CALL raise_tree(Trim(project_name), root)

        !------------------------------------------------------------------------------
        ! Source branch / target branch
        !------------------------------------------------------------------------------
write (*,*) "IM INCLUDING THE meta_param"
        CALL include_branch_into_branch(s_b=meta_para, t_b=root, blind=.TRUE.)

        !------------------------------------------------------------------------------
        ! Read an existing input tree (with microfocus ct data).
        !------------------------------------------------------------------------------
        ! Load puredat tree of micro-CT data and calculate the global
        ! parameters of the domain decomposition
        !------------------------------------------------------------------------------
        pro_path = muCT_pd_path
        pro_name = muCT_pd_name

        phi_tree = read_tree()

        !------------------------------------------------------------------------------
        ! Set project name and path of global domain decomposition     
        !------------------------------------------------------------------------------
        pro_path = outpath
        pro_name = project_name

        allocate(ddc)
        ddc = calc_general_ddc_params(bone%phdsize, phi_tree)
        
        CALL include_branch_into_branch(s_b=ddc, t_b=root, blind=.TRUE.)
        
    ELSE ! restart = 'Y'         

        !------------------------------------------------------------------------------
        ! Read an existing output tree (with microfocus ct data).
        ! Characteristics of PureDat and MeRaDat require this chicken dance...
        !------------------------------------------------------------------------------
        ! pro_path = TRIM(in%path)
        ! pro_name = TRIM(in%bsnm)
write(*,*) "root path: "//TRIM(pro_path)//TRIM(pro_name)

        root = read_tree()

        ! pro_path = outpath
        ! pro_name = project_name
write(*,*) "root path: "//TRIM(pro_path)//TRIM(pro_name)
        !------------------------------------------------------------------------------
        
        If (out_amount == "DEBUG") THEN 
            WRITE(un_lf, fmt_dbg_sep)
            Write(un_lf,'(A)') "root right after restart read"
            CALL log_tree(root, un_lf,.FALSE.)
            WRITE(un_lf, fmt_dbg_sep)
            flush(un_lf)
        END If
       
        !------------------------------------------------------------------------------
        ! This CALLing sequence is only valid since "Averaged Material 
        ! Properties" only contains r8 data added at the end of the 
        ! r8-stream. More correct would be a routine that ensures data
        ! integrity and compresses potentially missing stream data in
        ! an efficient way.
        !------------------------------------------------------------------------------
        CALL delete_branch_from_branch("Averaged Material Properties", root, dsize)

        CALL get_stream_size(root, dsize)
        root%streams%dim_st = dsize
        root%streams%ii_st  = dsize + 1

        CALL read_streams(root)

        CALL connect_pointers(root%streams, root)

        CALL search_branch("Global domain decomposition", root, ddc, success)

        IF (.NOT. success) THEN
            mssg = "No branch named 'Global domain decomposition', however a restart was requested."
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF

        CALL search_branch("Input parameters", root, meta_para, success)

        IF (.NOT. success) then
            mssg = "No branch named 'Input parameters', however a restart was requested."
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF

        !------------------------------------------------------------------------------
        ! Reset Output amount and Restart in loaded param branch
        ! Hardcoded, implicitly given order of the leaves. 
        ! DO NOT CHANGE INDICES WITHOUT MODYFING THE »add_leaf_to_branch« SEQUENCES.
        !------------------------------------------------------------------------------
        meta_para%leaves(15)%p_char = str_to_char(out_amount)
        meta_para%leaves(16)%p_char = "Y"

        If (out_amount == "DEBUG") THEN 
            Write(un_lf,fmt_dbg_sep)
            Write(un_lf,'(A)') "root right after restart and deletion of avg mat props branch"
            CALL log_tree(root, un_lf,.FALSE.)
            Write(un_lf, fmt_dbg_sep)
            flush(un_lf)
        END If

    END IF ! restart == Yes/No
  
    !------------------------------------------------------------------------------
    ! Initial setup of the worker_is_active tracker
    !------------------------------------------------------------------------------
    ! Tracker writes -9.999.999 if: Space within tracking is not occupied by a domain
    ! Tracker writes    +Number if: Domain is computed successfully
    !
    ! If the position within the tracking file is implicitely connected to the 
    ! domain number, changes to the domain range before restart will crash the 
    ! computations (!)
    !------------------------------------------------------------------------------
    IF (stat_exists) THEN
        OPEN(aun, FILE=TRIM(out%p_n_bsnm)//".status", &
            ACTION="READWRITE", STATUS="OLD", ACCESS="STREAM")
        
        READ(aun) Domain_stats
        REWIND(aun)
 
        !------------------------------------------------------------------------------
        ! Check whether computation will use resources properly.
        !------------------------------------------------------------------------------
        DO ii=1, SIZE(domain_stats)
            IF (domain_stats(ii) >= 0) computed_domains = computed_domains + 1 
        END DO

        IF (No_of_domains == computed_domains) THEN 
            mssg = "Job is already finished. No restart required."
            CALL print_err_stop_slaves(mssg, "message"); GOTO 1000
        END IF 
        
    ELSE
        OPEN(aun, FILE=TRIM(out%p_n_bsnm)//".status", &
            ACTION="READWRITE", STATUS="NEW", ACCESS="STREAM")
        WRITE(aun) Domain_stats
        FLUSH(aun)
    END IF ! (stat_exists)

    !------------------------------------------------------------------------------
    ! Check whether computation will use resources properly.
    ! Check outside Yes/No if/else because the amount of remaining domains
    ! determines the (re-)start of the job.
    !------------------------------------------------------------------------------
    IF ((No_of_domains - computed_domains) * parts /= size_mpi - 1) THEN
        mssg = "Remaining domains * parts /= /= size_mpi - 1"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF
    
    CALL add_branch_to_branch(root, result_branch)
    CALL raise_branch("Averaged Material Properties", 0_pd_ik, 18_pd_ik, result_branch)
    
    !------------------------------------------------------------------------------
    ! To use pd_store, the memory must be allocated by raising the leaves.
    !------------------------------------------------------------------------------
    CALL raise_leaves(no_leaves = 18_pd_ik, &
        desc = [ &
        "Domain forces                                     ", &
        "Effective numerical stiffness                     ", &
        "Symmetry deviation - effective numerical stiffness", &
        "Averaged stresses                                 ", &
        "Averaged strains                                  ", &
        "Effective stiffness                               ", &
        "Symmetry deviation - effective stiffness          ", &
        "Averaged Effective stiffness                      ", &
        "Symmetry deviation - Averaged effective stiffness ", &
        "Rotation Angle CR_1                               ", &
        "Rotation Vector CR_1                              ", &
        "Final coordinate system CR_1                      ", &
        "Optimized Effective stiffness CR_1                ", &
        "Rotation Angle CR_2                               ", &
        "Rotation Vector CR_2                              ", &
        "Final coordinate system CR_2                      ", &
        "Optimized Effective stiffness CR_2                ", &
        "Effective density                                 "], &
        dat_ty = [(5_1, ii=1,18)], &
        dat_no = [ &
        No_of_domains * 24*24, No_of_domains * 24*24, No_of_domains        , &
        No_of_domains *  6*24, No_of_domains *  6*24, No_of_domains *  6* 6, &
        No_of_domains        , No_of_domains *  6* 6, No_of_domains        , &
        No_of_domains        , No_of_domains *     3, No_of_domains *     9, &
        No_of_domains *  6* 6, No_of_domains        , No_of_domains *     3, &
        No_of_domains *     9, No_of_domains *  6* 6, No_of_domains       ], &
        branch = result_branch)
    
    result_branch%leaves(:)%pstat = -1
    
    CALL set_bounds_in_branch(root, root%streams)
    
    If (out_amount == "DEBUG") THEN 
        Write(un_lf, fmt_dbg_sep)
        Write(un_lf, '(A)') "root right before serialisation"
        CALL log_tree(root,un_lf,.True.)
        Write(un_lf, fmt_dbg_sep)
        flush(un_lf)
    END If

    !------------------------------------------------------------------------------
    ! serialize root branch
    !------------------------------------------------------------------------------
    CALL serialize_branch(root, serial_root, serial_root_size, .TRUE.)

    !------------------------------------------------------------------------------
    ! Get number of domains on each axis
    !------------------------------------------------------------------------------
    CALL pd_get(ddc,"nn_D", nn_D)

    !------------------------------------------------------------------------------
    ! Init Domain Cross Reference and domain paths
    !------------------------------------------------------------------------------
    dc = 0
    domain_path(0) = Trim(pro_path)//Trim(pro_name)//"_domains" ! "_results"

    nn = 1
    Do kk = xa_d(3), xe_d(3)
    Do jj = xa_d(2), xe_d(2)
    Do ii = xa_d(1), xe_d(1)
        dc         = dc + 1_mik
        path_count = dc / 2_mik

        Write(domain_path(dc),'(A,"/",I0)')Trim(domain_path(path_count)),dc

        Domains(nn) = ii + jj * nn_D(1) + kk * nn_D(1)*nn_D(2)

        nn = nn + 1_mik
    End Do
    End Do
    End Do

    !------------------------------------------------------------------------------
    ! Generate worker_is_active_List
    !------------------------------------------------------------------------------
    Allocate(worker_is_active(size_mpi-1), stat=alloc_stat)
    CALL alloc_err("worker_is_active_List", alloc_stat)

    worker_is_active=0

    Allocate(workers_assigned_domains(size_mpi-1), stat=alloc_stat)
    CALL alloc_err("workers_assigned_domains",alloc_stat)

    workers_assigned_domains = 0
    CALL End_Timer("Init Process")

    !------------------------------------------------------------------------------
    ! Start Workers
    !------------------------------------------------------------------------------
    CALL Start_Timer("Broadcast Init meta_para")

    CALL mpi_bcast(outpath,      INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(project_name, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

    ! write(un_lf, FMT_MSG_xAI0) "Broadcasting serialized root of size [Byte] ", serial_root_size*8

    CALL mpi_bcast(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

    CALL mpi_bcast(serial_root, INT(serial_root_size,mik), MPI_INTEGER8, 0_mik,&
        MPI_COMM_WORLD, ierr)

    CALL End_Timer("Broadcast Init meta_para")

    !------------------------------------------------------------------------------
    ! Execute collective mpi_comm_split. Since mpi_comm_world rank 0 is
    ! the head master worker_comm is not needed and it should not be in
    ! any worker group and communicator. With MPI_UNDEFINED passed as
    ! color worker_comm gets the value MPI_COMM_NULL
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank_mpi, worker_comm, ierr)
    IF (ierr /=0) THEN
        mssg = "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

!------------------------------------------------------------------------------
! Ranks > 0 -- Worker slaves
!------------------------------------------------------------------------------
Else
    CALL Start_Timer("Broadcast Init meta_para")
     
    !------------------------------------------------------------------------------
    ! Broadcast recieve init parameters
    !------------------------------------------------------------------------------
    CALL mpi_bcast(outpath      , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(project_name , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
    
    pro_path = outpath
    pro_name = project_name

    !------------------------------------------------------------------------------
    ! Serial_root_size == -1 ==> Signal that No_of_domains < size_mpi-1
    !------------------------------------------------------------------------------
    If (serial_root_size == -1) then
        CALL End_Timer("Broadcast Init meta_para")
        Goto 1000
    End If
    
    Allocate(serial_root(serial_root_size))
    
    CALL mpi_bcast(serial_root, INT(serial_root_size, mik), MPI_INTEGER8, 0_mik, &
        MPI_COMM_WORLD, ierr)

    CALL End_Timer("Broadcast Init meta_para")

    !------------------------------------------------------------------------------
    ! Deserialize root branch
    !------------------------------------------------------------------------------
    CALL Start_Timer("Deserialize root branch")

    CALL deserialize_branch(root, serial_root, .TRUE.)
    CALL assign_pd_root(root)
    CALL set_bounds_in_branch(root,root%streams)

    !------------------------------------------------------------------------------
    ! Commented out since out_amount is a parametrized global variable 
    ! CALL pd_get(root%branches(1),"Output amount", char_arr)
    ! out_amount = char_to_str(char_arr)
    ! deallocate(char_arr)
    !------------------------------------------------------------------------------
    CALL pd_get(root%branches(1),"Restart", char_arr)
    restart = char_to_str(char_arr)
    deallocate(char_arr)

    CALL End_Timer("Deserialize root branch")

    !------------------------------------------------------------------------------
    ! Init Domain Cross Reference
    !------------------------------------------------------------------------------
    write(*,*) "XA_D", xa_d
    CALL pd_get(root%branches(1), "Lower bounds of selected domain range", xa_d, 3)
    CALL pd_get(root%branches(1), "Upper bounds of selected domain range", xe_d, 3)
    
    write(*,*) "XA_D", xa_d
    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)
    
    Allocate(Domains(No_of_domains),stat=alloc_stat)
    CALL alloc_err("Domains",alloc_stat)

    CALL pd_get(root%branches(2),"nn_D", nn_D)
    
    !------------------------------------------------------------------------------
    ! Linear vs ddc (3D) domain numbers. nn is a purely internal variable.
    !------------------------------------------------------------------------------
    nn = 1
    Do kk = xa_d(3), xe_d(3)
    Do jj = xa_d(2), xe_d(2)
    Do ii = xa_d(1), xe_d(1)
        Domains(nn) = ii + jj * nn_D(1) + kk * nn_D(1)*nn_D(2)
        nn = nn + 1_mik
    End Do
    End Do
    End Do

    CALL pd_get(root%branches(1), "No of mesh parts per subdomain", parts)

    !------------------------------------------------------------------------------
    ! All Worker Ranks -- Init worker Communicators
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, Int((rank_mpi-1)/parts, mik), rank_mpi, worker_comm, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", ierr)
    
    CALL MPI_COMM_RANK(WORKER_COMM, worker_rank_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't retrieve worker_rank_mpi", ierr)

    CALL MPI_COMM_SIZE(WORKER_COMM, worker_size_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't retrieve worker_size_mpi", ierr)

    !------------------------------------------------------------------------------
    ! This sets the options for PETSc in-core. To alter the options
    ! add them in Set_PETSc_Options in Module pets_opt in file
    ! f-src/mod_parameters.f90
    !------------------------------------------------------------------------------
    CALL Set_PETSc_Options()

    PETSC_COMM_WORLD = worker_comm
    CALL PetscInitialize(PETSC_NULL_CHARACTER, petsc_ierr)
End If

WRITE_ROOT_COMM = MPI_COMM_WORLD

!------------------------------------------------------------------------------
! Open stream files
!------------------------------------------------------------------------------
Call Open_Stream_Files(root%streams, "write", "replace", fh_mpi_root, WRITE_ROOT_COMM)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
Allocate(req_list(size_mpi-1), stat=alloc_stat)
CALL alloc_err("req_list", alloc_stat)
req_list=1 

Allocate(statuses_mpi(MPI_STATUS_SIZE, size_mpi-1), stat=alloc_stat)
CALL alloc_err("statuses_mpi", alloc_stat)

!------------------------------------------------------------------------------
! Global Rank 0 -- Process master Start working process
!------------------------------------------------------------------------------
If (rank_mpi==0) Then
    !------------------------------------------------------------------------------
    ! Only global rank 0 -- Init MPI IO-System
    !------------------------------------------------------------------------------
    CALL get_stream_size(root, dsize)

    If (out_amount == "DEBUG") &
        write(fh_mon,FMT_MSG_AxI0) "On rank zero, stream sizes: ", dsize
    
    nn = 1; mii = 1

    !  Call Start_Timer("Write Root Branch")
    !  call store_parallel_branch(root, fh_mpi_root)
    !  Call Write_Tree(root)
    !  Call End_Timer("Write Root Branch")

    !------------------------------------------------------------------------------
    ! mii is incremented by mii = mii + parts
    ! mii --> Worker master
    !------------------------------------------------------------------------------
    DO WHILE (nn <= No_of_domains) 
    
    ! (mii <= (size_mpi-1_mik))

        !------------------------------------------------------------------------------
        ! Skip this DO WHILE cycle if Domain_stat of the domain (nn) > 0 and therefore
        ! marked as "already computed".
        !------------------------------------------------------------------------------
        If (Domain_stats(nn) >= 0) then
            nn = nn + 1_ik
            cycle
        End If

        !------------------------------------------------------------------------------
        ! Send information for all parts to the corresponding ranks.
        ! From global "first worker rank of specific domain" to "global last worker 
        ! rank of specific domain"
        !------------------------------------------------------------------------------
        ! jj --> Worker
        !------------------------------------------------------------------------------
        Do mjj = mii, mii + parts-1
            
            req_list(mjj) = 1_mik
            worker_is_active(mjj) = 1_mik
            workers_assigned_domains(mjj) = nn

            !------------------------------------------------------------------------------
            ! Start worker
            !------------------------------------------------------------------------------
            CALL mpi_send(worker_is_active(mjj), 1_mik, MPI_INTEGER4, mjj, mjj, MPI_COMM_WORLD, ierr)
            CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)

            !------------------------------------------------------------------------------
            ! Send Domain number
            !------------------------------------------------------------------------------
            CALL mpi_send(nn, 1_mik, MPI_INTEGER8, mjj, mjj, MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", ierr)
            
            if (out_amount /= "PRODUCTION") then
               CALL mpi_send(domain_path(nn), Int(4*mcl,mik), MPI_CHARACTER, mjj, mjj, MPI_COMM_WORLD,ierr)
               CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", ierr)
            End if

            !------------------------------------------------------------------------------
            ! Worker has finished
            !------------------------------------------------------------------------------
            CALL MPI_IRECV(worker_is_active(mjj), 1_mik, MPI_INTEGER4, mjj, mjj, MPI_COMM_WORLD, req_list(mjj), ierr)
            CALL print_err_stop(std_out, "MPI_IRECV of worker_is_active(mjj) didn't succeed", ierr)
        End Do
         
        !------------------------------------------------------------------------------
        ! Log to monitor file. Only for master worker!
        !------------------------------------------------------------------------------
        IF(out_amount == "DEBUG") THEN
            WRITE(fh_mon, FMT_MSG_xAI0) "Domain ", nn, " at Ranks ", mii, " to ", mii + parts-1
            flush(fh_mon)
        END IF


        !------------------------------------------------------------------------------
        ! After all ranks have their first work package
        !------------------------------------------------------------------------------
        IF ((size_mpi-1_mik)/parts >= nn) THEN
            !------------------------------------------------------------------------------
            ! Only the master worker is allowed to send a 'finished' flag. 
            ! Otherwise, index mii will screw up.
            !------------------------------------------------------------------------------
            CALL MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
write(*,*) "Request list during send: ", req_list
            CALL print_err_stop(std_out, &
            "MPI_WAITANY on req_list for IRECV of worker_is_active(mii) didn't succeed", ierr)

            IF(finished /= MPI_UNDEFINED) THEN
                mii = finished

                Domain_stats(nn) = Domains(nn)

                !------------------------------------------------------------------------------
                ! Job is finished.Now, write domain number to position of linear domain number.
                !------------------------------------------------------------------------------


write(*,*) "This nn: ", nn   
write(*,*) "This domain number: ", Domains(nn)   
                write(aun, pos=((nn-1) * ik + 1)) INT(Domains(nn), KIND=ik)
                flush(aun)
            END IF
        END IF

        !------------------------------------------------------------------------------
        ! iterate with step with of parts since the DO WHILE acts on all parts
        ! in one iteration with a dedicated loop (jj = mii, mii+parts-1)
        !------------------------------------------------------------------------------
        mii = mii + Int(parts,mik)

        !------------------------------------------------------------------------------
        ! Iterate over domain
        !------------------------------------------------------------------------------
        nn = nn + 1_mik
    End Do
    
    !------------------------------------------------------------------------------
    ! Add leaf with analyzed cube numbers to root
    !------------------------------------------------------------------------------
    If ( restart == "N" ) then
        CALL add_leaf_to_branch   (root, "Number of domains", No_of_domains, Domains)
        CALL set_bounds_in_branch (root, root%streams)
    End If

    write(*,*) "SOLLTE SCHREIBEN"

    !------------------------------------------------------------------------------
    ! store_parallel_branch(root, fh_mpi_root) will not print anything?!
    !------------------------------------------------------------------------------
    CALL Start_Timer("Write Root Branch")

    CALL store_parallel_branch(root%branches(1), fh_mpi_root)
    CALL store_parallel_branch(root%branches(2), fh_mpi_root)
    CALL Write_Tree(root)

    CALL End_Timer("Write Root Branch")

    
    !------------------------------------------------------------------------------
    ! Write Root header and input parameters
    ! Stream size small, data remain in rank-directories
    !------------------------------------------------------------------------------
    ! CALL Start_Timer("Write Root Branch")
    
    ! CALL open_stream_files(root%streams, "write", "replace")
    ! Save from workers?
    ! CALL store_branch(root, root%branches%streams, .TRUE.)
    
    ! CALL Write_Tree(root)
    ! CALL close_stream_files(root)

    ! CALL End_Timer("Write Root Branch")
    
    ! FHS??
    
        ! CALL store_parallel_branch(root, FH_MPI)
    
! write(*,*) "Rank: ", rank_mpi, " root%streams%dim_st = ", root%streams%dim_st
! write(*,*) "Rank: ", rank_mpi, " root%streams%no_branches = ", root%streams%no_branches
! write(*,*) "Rank: ", rank_mpi, " root%streams%no_leaves = ", root%streams%no_leaves
! write(*,*) "Rank: ", rank_mpi, " pro_path:", pro_path
! write(*,*) "Rank: ", rank_mpi, " pro_name:", pro_name



    ! CALL store_parallel_branch(root%branches(3), fh_mpi_worker) ! Branch with 'Averaged Material Properties'

        
    !------------------------------------------------------------------------------
    ! A CALL to MPI_Waitany can be used to wait for 
    ! the completion of one out of several requests.
    ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Waitany.3.php
    !------------------------------------------------------------------------------
    ! CALL MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
    ! CALL print_err_stop(std_out, &
    ! "MPI_WAITANY on req_list for IRECV of worker_is_active(ii) didn't succeed", ierr)

    ! IF(finished /= MPI_UNDEFINED) THEN
    !     ii = finished

    !     Domain_stats(nn) = Domains(nn)

    !     !------------------------------------------------------------------------------
    !     ! Job is finished(?)
    !     !------------------------------------------------------------------------------
    !     write(aun, pos=(Domain_stats(nn)-1) * ik + 1) INT(Domains(nn), KIND=ik)
    !     flush(aun)
    ! END IF

    ! !------------------------------------------------------------------------------
    ! ! After all worker masters have a job, the succeeding iterations start.
    ! !------------------------------------------------------------------------------
    ! Do While (nn <= No_of_domains)

    !     !------------------------------------------------------------------------------
    !     ! Skip this DO WHILE cycle if Domain_stat of the domain (nn) > 0 and therefore
    !     ! marked as "already computed".
    !     !------------------------------------------------------------------------------
    !     If (Domain_stats(nn) > 0) then
    !         nn = nn + 1_ik
    !         cycle
    !     End If

    !     Do jj = ii, ii + parts- 1

    !         worker_is_active(jj) = 1_mik
    !         workers_assigned_domains(jj) = nn

    !         CALL mpi_send(worker_is_active(jj), 1_mik, MPI_INTEGER4 , Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
    !         CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)

    !         CALL mpi_send(nn, 1_mik, MPI_INTEGER8, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
    !         CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", ierr)

    !         if (out_amount /= "PRODUCTION") then
    !             CALL mpi_send(domain_path(nn), Int(4_mik*mcl,mik), &
    !                     MPI_CHARACTER, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
    !             CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", ierr)
    !         End if
    !     End Do
         
    !     !------------------------------------------------------------------------------
    !     ! Log to monitor file. Only for master worker!
    !     !------------------------------------------------------------------------------
    !     IF(out_amount == "DEBUG") THEN
    !         WRITE(fh_mon, FMT_MSG_xAI0) "Domain ", nn, " at Rank ", ii
    !         flush(fh_mon)
    !     END IF
        
    !     !------------------------------------------------------------------------------
    !     ! Worker confirmes start
    !     !------------------------------------------------------------------------------
    !     CALL MPI_IRECV(worker_is_active(ii), 1_mik, MPI_INTEGER4, Int(ii,mik), &
    !                 Int(ii,mik), MPI_COMM_WORLD, REQ_LIST(ii), IERR)
    !     CALL print_err_stop(std_out, "MPI_IRECV of worker_is_active(ii) didn't succeed", ierr)

    !     CALL MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
    !     CALL print_err_stop(std_out, "MPI_WAITANY for IRECV of worker_is_active(ii) didn't succeed", ierr)

    !     IF(finished /= MPI_UNDEFINED) THEN
    !         ii = finished

    !     END IF

    !     If (out_amount == "DEBUG") WRITE(fh_mon, FMT_MSG_xAI0) "Domain ", nn, " finished"
        
    !     Domain_stats(nn) = Domains(nn)
    !     ! Domain_stats(workers_assigned_domains(ii)) = worker_is_active(ii)

    !     !------------------------------------------------------------------------------
    !     ! Job is finished(?)
    !     !------------------------------------------------------------------------------
    !     write(aun, pos=(Domain_stats(nn)-1) * ik + 1) INT(Domains(nn), KIND=ik)
    !     flush(aun)
        
    !     nn = nn + 1_mik

    ! End Do
write(*,*) "Request list: ", req_list
    !------------------------------------------------------------------------------
    ! Wait for all workers to return
    !------------------------------------------------------------------------------
    CALL MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
    CALL print_err_stop(std_out, &
        "MPI_WAITALL on req_list for IRECV of worker_is_active(ii) didn't succeed", ierr)
Write(*,*) "Shut the fuck up"
    !------------------------------------------------------------------------------
    ! TODO refactor domain_cross reference from size size_mpi-1 to
    ! (size_mpi-1)/parts
    !------------------------------------------------------------------------------
! HOPE IT IS NOT REQUIRED?

    ! Do ii = 1, size_mpi-1, parts
    !     write(aun, pos=(workers_assigned_domains(ii)-1)*ik+1) INT(worker_is_active(ii), KIND=ik)
    ! End Do

    ! flush(aun)

    !------------------------------------------------------------------------------
    ! Stop workers
    !------------------------------------------------------------------------------
    worker_is_active = -1_mik
    
    Do mii = 1_mik, size_mpi-1_mik
        CALL mpi_send(worker_is_active(mii), 1_mik, MPI_INTEGER4, mii, mii, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)
    End Do
Write(*,*) "oh no"

!------------------------------------------------------------------------------
! Global ranks > 0 -- Workers
!------------------------------------------------------------------------------
Else

    IF (worker_rank_mpi == 0) THEN
        !------------------------------------------------------------------------------
        ! Extend project_name and outpath with rank - only for domain specific master!
        !------------------------------------------------------------------------------
        WRITE(outpath,'(A,A,I7.7,A)')    TRIM(outpath), "Rank_", rank_mpi, "/"
        WRITE(project_name,'(A,A,I7.7)') TRIM(project_name), "_", rank_mpi

        !------------------------------------------------------------------------------
        ! Prepare output directory via CALLing the c function.
        ! File exists if stat_c_int = 0 
        !------------------------------------------------------------------------------
        c_char_array(1:LEN(TRIM(outpath)//CHAR(0))) = str_to_char(TRIM(outpath)//CHAR(0))
        CALL Stat_Dir(c_char_array, stat_c_int)

        IF(stat_c_int /= 0) THEN

            CALL execute_command_line("mkdir -p "//TRIM(outpath),CMDSTAT=stat)

            IF(stat /= 0) CALL print_err_stop(std_out, &
                'Could not execute syscall »mkdir -p '//trim(outpath)//'«.', 1)

            CALL Stat_Dir(c_char_array, stat_c_int)

            IF(stat_c_int /= 0) THEN
                mssg = 'Could not create the output directory »'//TRIM(outpath)//'«.'
                CALL print_err_stop(std_out, mssg, 1)
            END IF
        END IF

        CALL link_start(link_name, .True., .True.)
    END IF 

    !------------------------------------------------------------------------------
    ! Send new directory/filename
    ! It is the Rank_ and not the root directory! Do not delete this sequence!
    ! Root dir was send via broadcast, because the worker masters are not 
    ! determined at the beginning of the program.
    !------------------------------------------------------------------------------
    CALL MPI_BCAST(outpath,      INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)
    CALL MPI_BCAST(project_name, INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)
write(*,*)"Test0"

    !------------------------------------------------------------------------------
    ! Open Stream files to write data to master-worker rank directories
    !------------------------------------------------------------------------------
    pro_path = outpath
    pro_name = project_name     

    CALL set_stream_filenames(root%streams)
    CALL open_stream_files(root%streams, "write", "replace", fh_mpi_worker, WORKER_COMM)

    !------------------------------------------------------------------------------
    ! Worker Loop
    !------------------------------------------------------------------------------
    Do
        !------------------------------------------------------------------------------
        ! Start workers
        !------------------------------------------------------------------------------
        CALL MPI_RECV(active, 1_mik, MPI_INTEGER4, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Active didn't succseed", ierr)

        !------------------------------------------------------------------------------
        ! Stop workers
        !------------------------------------------------------------------------------
        IF (active == -1) EXIT

        !------------------------------------------------------------------------------
        ! Receive domain numbers
        !------------------------------------------------------------------------------
        CALL mpi_recv(nn, 1_mik, MPI_INTEGER8, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Domain didn't succeed", ierr)
write(*,*)"Test_nn"

        Domain = Domains(nn)

        IF (out_amount /= "PRODUCTION") THEN
           !------------------------------------------------------------------------------
           ! Receive Job_Dir
           !------------------------------------------------------------------------------
           CALL mpi_recv(job_dir, 4_mik*int(mcl,mik), mpi_character, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
           CALL print_err_stop(std_out, "MPI_RECV on Domain path didn't succeed", ierr)
        ELSE
           job_dir = outpath
        End if
        
        IF (job_dir(len(job_dir):len(job_dir)) /= "/") job_dir = trim(job_dir)//"/"

write(*,*)"Test1"
        !======================================================================
        CALL exec_single_domain(root, nn, Domain, job_dir, Active, fh_mpi_worker, &
             worker_rank_mpi, worker_size_mpi, worker_comm)
        !======================================================================
write(*,*)"Test2"

        !------------------------------------------------------------------------------
        ! Organize Results
        !------------------------------------------------------------------------------
        IF (worker_rank_mpi==0) THEN
            CALL Start_Timer("Write Worker Root Branch")

            ! root%streams%ii_st  = 1
            ! root%streams%dim_st = 0
            ! CALL reset_bounds_in_branch(root, root%streams)
            
            ! root%streams%ii_st  = 1
            ! root%streams%dim_st = 0
            ! CALL homogenize_branch(root, root%streams)

            !------------------------------------------------------------------------------
            ! Set these variables again...
            !------------------------------------------------------------------------------
            pro_path = outpath
            pro_name = project_name     
                                   
            !------------------------------------------------------------------------------
            ! Store header file
            !------------------------------------------------------------------------------
            ! CALL Write_Tree(root)
            CALL Write_Tree(root%branches(3)) ! Branch with 'Averaged Material Properties'
            ! ../../../bin/pd_dump_leaf_x86_64 $PWD/ results_0000001 7
        
        
        
        
        DO ii = 1, SIZE(root%branches)
            write(*, '(A, I0, 2A, T80, I0, T84, A)') &
                "Branch(", ii, ") of the tree: ", &
                TRIM(root%branches(ii)%desc), &
                SIZE(root%branches(ii)%leaves), " leaves."
        END DO



            CALL End_Timer("Write Worker Root Branch")
        END IF

        !------------------------------------------------------------------------------
        ! Store stream files. MPI-Parallel. All must work together.
        !------------------------------------------------------------------------------
        ! CALL store_parallel_branch(root, fh_mpi_worker)
        ! ../../../bin/pd_dump_Eleaf_x86_64 $PWD/ results_0000001 35
        CALL store_parallel_branch(root%branches(3), fh_mpi_worker) ! Branch with 'Averaged Material Properties'
        CALL close_stream_files(root%streams, fh_mpi_worker)

        active = 0_mik

        CALL MPI_ISEND(active, 1_mik, MPI_INTEGER4, 0_mik, rank_mpi, MPI_COMM_WORLD, request, ierr)
        CALL print_err_stop(std_out, "MPI_ISEND on Active didn't succeed", ierr)

        CALL MPI_WAIT(request, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_WAIT on request for ISEND Active didn't succeed", ierr)

    End Do

    CALL PetscFinalize(petsc_ierr) 
    
End If

write(*,*) "my_rank: ", rank_mpi, " after goto 1000"

CALL mpi_bcast(pro_path, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL mpi_bcast(pro_name, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

write(*,*) "Rank: ", rank_mpi, " pro_path:", pro_path
write(*,*) "Rank: ", rank_mpi, " pro_name:", pro_name

IF(rank_mpi == 0) THEN

    If (out_amount == "DEBUG") THEN
        Write(un_lf, fmt_dbg_sep)
        Write(un_lf, fmt_MSG_xAI0) "Final Root pointer proc", rank_mpi
        CALL log_tree(root, un_lf, .True.)
        Write(un_lf, fmt_dbg_sep)

        !------------------------------------------------------------------------------
        ! Give information about the tree structure
        !------------------------------------------------------------------------------
        DO ii = 1, SIZE(root%branches)
            write(un_lf, '(A, I0, 2A, T80, I0, T84, A)') &
                "Branch(", ii, ") of the tree: ", &
                TRIM(root%branches(ii)%desc), &
                SIZE(root%branches(ii)%leaves), " leaves."
        END DO
    END If

    CALL link_end(link_name,.True.)

    CALL meta_signing(binary)
    CALL meta_close(size_mpi)

    If (out_amount == "DEBUG") CALL meta_stop_ascii(fh_log, log_suf)

    CALL meta_stop_ascii(fh_mon, mon_suf)
    CALL meta_stop_ascii(aun   , ".status")

    IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')
END IF ! (rank_mpi == 0)

1000 Continue

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", ierr)

End Program main_struct_process