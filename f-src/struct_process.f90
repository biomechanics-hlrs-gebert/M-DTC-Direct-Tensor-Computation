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
  Subroutine exec_single_domain(root, lin_nn, nn, job_dir, Active, fh_mpi, &
       rank_mpi, size_mpi, comm_mpi)

    TYPE(materialcard) :: bone

    LOGICAL, PARAMETER :: DEBUG=.TRUE.
    
    Integer(kind=mik), intent(out) :: Active

    Type(tBranch), Intent(inOut) :: root
    Integer(kind=ik), Intent(in) :: nn, lin_nn
    Character(LEN=*), Intent(in) :: job_dir
    Integer(kind=mik), Intent(In) :: rank_mpi, size_mpi, comm_mpi
    Integer(kind=mik), Intent(In), Dimension(no_streams) :: fh_mpi
    
    !----------------------------------------------------------------
    Integer(kind=mik), Dimension(MPI_STATUS_SIZE)   :: status_mpi
    Type(tBranch), pointer :: bb, db, pb, mb, meta_para,resb
    Integer(kind=mik) :: ierr
    
    Character(Len=mcl) :: desc, mesh_desc, filename
    Character(Len=mcl) :: elt_micro
    
    Character, Dimension(4*mcl) :: c_char_array

    Integer(kind=c_int) :: stat_c_int
    Integer             :: stat
    Integer             :: umon
 
    Character(len=mcl)  :: env_var
 
    Integer(kind=ik)    :: m_size
 
    Character(len=9)    :: nn_char
 
    logical             :: success
 
    Character(len=mcl)  :: timer_name, domain_desc, part_desc

    Integer(kind=mik)         :: petsc_ierr
    Type(tMat)                :: AA, AA_org
    Type(tVec)                :: XX
    Type(tVec), Dimension(24) :: FF
    TYPE(tPETScViewer)        :: PetscViewer
    Type(tKSP)                :: KSP
    Integer(Kind=ik)          :: Istart,Iend, parts_per_subdomain, IVstart, IVend
    Integer(Kind=ik), Dimension(:), Allocatable :: nodes_in_mesh

    Integer(kind=pd_ik), Dimension(:), Allocatable :: serial_pb
    Integer(kind=pd_ik)                            :: serial_pb_size

    Integer(Kind=pd_ik)  , Dimension(no_streams)  :: no_data
    Integer(kind=ik)                              :: nn_elems, ii, jj, kk, id
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
           
    !--------------------------------------------------------------------------

    !** Init Activity status ****
    Active = 0_mik

    write(nn_char,'(I0)')nn

    !** Get basic infos ------------------------------------------
    CALL Search_branch("Input parameters", root, meta_para, success, DEBUG)
    CALL pd_get(meta_para, "No of mesh parts per subdomain", parts_per_subdomain)
    CALL pd_get(meta_para, "Physical domain size" , bone%phdsize, 3)
    CALL pd_get(meta_para, "Grid spacings" , bone%delta, 3)
    CALL pd_get(meta_para, "Young_s modulus" , bone%E)
    CALL pd_get(meta_para, "Poisson_s ratio" , bone%nu)
     
    !****************************************************************************
    !** Rank = 0 -- Local master of comm_mpi ************************************
    !****************************************************************************
    If ( rank_mpi == 0) then

       !**************************************************************************
       !** Create job directory in case of non-production run
       !**************************************************************************
       If (out_amount /= "PRODUCTION") then

          c_char_array(1:len(Trim(job_dir)//Char(0))) = str_to_char(Trim(job_dir)//Char(0))

          Call Stat_Dir(c_char_array, stat_c_int)

          If ( stat_c_int /= 0 ) Then

             CALL execute_command_line("mkdir -p "//trim(job_dir), CMDSTAT=stat)

             IF(stat /= 0) CALL print_err_stop(std_out, "Couldn't execute syscall mkpir -p "//TRIM(job_dir), 1)

             CALL Stat_Dir(c_char_array, stat_c_int)

             IF(stat_c_int /= 0) CALL print_err_stop(std_out, "Couldn't create directory "//TRIM(job_dir), 1)

          Else

             !***********************************************************************
             !usta = give_new_unit()
             !.... Check sta file

          End If

          Write(un_lf,FMT_MSG_SEP)

       End If

       !**************************************************************************
       !** Write log and monitor file
       !**************************************************************************
       Write(un_lf,FMT_MSG_xAI0) "Domain No. : ",nn
       Write(un_lf,FMT_MSG)      "Job_dir    : "//Trim(job_dir)

       CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)
 
       str = ''
       str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
       str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
       str = TRIM(str)//' '//timezone
       
       WRITE(un_lf, '(2A)') 'Start time: ', TRIM(str)

       Call get_environment_Variable("HOSTNAME", env_var)
       Write(un_lf,FMT_MSG) "Host       : "//Trim(env_var)

       umon = give_new_unit()

       Open(unit=umon, file=Trim(job_dir)//Trim(project_name)//".mon",action="write", &
            status="replace")

       Write(umon,FMT_MSG_SEP)
       Write(umon,FMT_MSG_xAI0) "Domain No. : ",nn

       !**************************************************************************
       !** Generate Geometry
       !**************************************************************************
       Select Case (timer_level)
       Case (1)
          timer_name = "+-- generate_geometry "//trim(nn_char)
       Case default
          timer_name = "generate_geometry"
       End Select

       Call start_timer(trim(timer_name), .FALSE.)
       Call generate_geometry(root, lin_nn, nn, job_dir, fh_mpi, success)
       if (.not.success) then
          write(*,FMT_WRN)"generate_geometry() returned .FALSE."
       End if
       Call end_timer(trim(timer_name))
     
       CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)
 
       str = ''
       str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
       str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
       str = TRIM(str)//' '//timezone
       
       WRITE(un_lf, '(2A)') 'End time: ', TRIM(str)

       close(umon)

       !** Look for the Domain branch ****************************************
       domain_desc=''
       Write(domain_desc,'(A,I0)')'Domain ',nn
       
       Call search_branch(trim(domain_desc), root, db, success, DEBUG)

       !** Get the no of nodes per part **************************************
       mesh_desc = ''
       Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_',nn
       
       Call search_branch(trim(mesh_desc), db, mb, success, DEBUG)
       call pd_get(mb, 'No of nodes in mesh',  nodes_in_mesh)

       !** Set the global matrix size ****************************************
       m_size = nodes_in_mesh(1) * 3

       Do ii = 1, parts_per_subdomain-1

          ! Look for the Part branch *************************
          domain_desc=''
          Write(part_desc,'(A,I0)')'Part_',ii
          
          !** serialize branch with mesh part *********************************
          Call search_branch(trim(part_desc), db, pb, success)
          If (.NOT. success) Then
            WRITE(mssg,'(2(A,I9))') "Exec_single_domain, rank ", rank_mpi ,": Branch of part ",ii
            CALL print_err_stop(std_out, mssg, 1)
          End If
          
          Call serialize_branch(pb,serial_pb,serial_pb_size,.TRUE.)

          Call mpi_send(serial_pb_size, 1_mik, MPI_INTEGER8, Int(ii,mik), Int(ii,mik), &
               COMM_MPI, ierr)
          Call mpi_send(serial_pb, INT(serial_pb_size,mik), MPI_INTEGER8, &
               Int(ii,mik), Int(ii,mik), COMM_MPI, ierr)

          Deallocate(serial_pb)
          
       End Do

       part_desc=''
       Write(part_desc,'(A,I0)')'Part_',parts_per_subdomain

       Call search_branch(trim(part_desc), db, pb, success, DEBUG)

       !** Broadcast matrix size. TODO could also be included into part branches.
       Call mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)

    !****************************************************************************
    !** Ranks > 0 -- Workers ****************************************************
    !****************************************************************************
    Else

       Call mpi_recv(serial_pb_size, 1_mik, mpi_integer8, 0_mik, &
            rank_mpi, COMM_MPI, status_mpi, ierr)

       if (allocated(serial_pb)) deallocate(serial_pb)
       
       Allocate(serial_pb(serial_pb_size))

       Call mpi_recv(serial_pb, INT(serial_pb_size,mik), mpi_integer8, &
            0_mik, rank_mpi, COMM_MPI, status_mpi, ierr)

       !** Deserialize part branch ****************************
       Call Start_Timer("Deserialize part branch branch")

       Allocate(pb)
       
       Call deserialize_branch(pb, serial_pb, .TRUE.)

       Call End_Timer("Deserialize part branch branch")

       Call mpi_bcast(m_size, 1_mik, MPI_INTEGER8, 0_mik, COMM_MPI, ierr)
              
    End If

    ! DEBUG INFORMATION
    If (out_amount == "DEBUG") THEN 
       Write(un_lf,fmt_dbg_sep)
       Write(un_lf,'(A)')"part branch right after deserialization"
       Call log_tree(pb,un_lf,.TRUE.)
       Write(un_lf,fmt_dbg_sep)
       flush(un_lf)
    END If
    ! DEBUG INFORMATION       


    !**************************************************************************
    !** Setup the linear System with a constant system matrix A. Once that  ***
    !** is done setup the multiple right hand sides and solve the linear    ***
    !** system multiple times. Save the solutions to calculate effective    ***
    !** stiffness matirces.                                                 ***
    !**************************************************************************
    IF (rank_mpi == 0) THEN   ! Sub Comm Master
        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- create_Stiffness_matrix "//TRIM(nn_char)
        CASE default
            timer_name = "create_Stiffness_matrix"
        End SELECT
        
        CALL start_timer(TRIM(timer_name), .FALSE.)
    END IF 


    !** Create Stiffness matrix **************************************
    call MatCreate(COMM_MPI, AA    , petsc_ierr)
    call MatCreate(COMM_MPI, AA_org, petsc_ierr)

    call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)
    call MatSetFromOptions(AA,petsc_ierr)
    call MatSetUp(AA,petsc_ierr)

    call MatSetSizes(AA_org,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)
    call MatSetFromOptions(AA_org,petsc_ierr)
    call MatSetUp(AA_org,petsc_ierr)

    !------------------------------------------------------------------------------
    ! End timer
    !------------------------------------------------------------------------------
    IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

    ! DEBUG INFORMATION
    If (out_amount == "DEBUG") THEN 
       call MatGetOwnershipRange(AA, Istart , Iend,  petsc_ierr)
       Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
            "MPI rank : ",rank_mpi,"| matrix ownership for ","A    ", ":",Istart," -- " , Iend

       call MatGetOwnershipRange(AA_org, Istart , Iend,  petsc_ierr)

       Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
       "MPI rank : ",rank_mpi,"| matrix ownership for ","A_org", ":",Istart," -- " , Iend
    End If
    ! DEBUG INFORMATION       

    !** pb%leaves(1) : Local  node ids   in part branch
    !** pb%leaves(2) : Coordinate values in part branch
    !** pb%leaves(3) : Global node ids   in part branch
    !** pb%leaves(4) : NOT USED ! (For Element Numbers)
    !** pb%leaves(5) : Topology
    !** No Cross reference global nid to local nid is currently
    !** needed since renumbering is
    !** deactivated in mod_mesh_partitioning.f90 : part_mesh

    !** Assemble matrix **********************************************

    Call pd_get(root%branches(1),'Element type  on micro scale',char_arr)
    elt_micro = char_to_str(char_arr)

    if (TRIM(elt_micro) == "HEX08") then

       K_loc_08 = Hexe08(bone)

       nn_elems = pb%leaves(5)%dat_no / 8
    
       Do ii = 1, nn_elems

          kk = 1
          
          !** Translate Topology to global dofs *************************
          Do jj = (ii-1)*8+1 , ii*8
             
             id = pb%leaves(5)%p_int8(jj)
             id = (id - 1) * 3
             
             idxm_08(kk)   = (id    ) !*cr
             idxm_08(kk+1) = (id + 1) !*cr
             idxm_08(kk+2) = (id + 2) !*cr
             
             kk = kk + 3
             
          End Do
          
          idxn_08 = idxm_08

          Call MatSetValues(AA, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)

          Call MatSetValues(AA_org, 24_8, idxm_08, 24_8 ,idxn_08, K_loc_08, ADD_VALUES, petsc_ierr)
          
       End Do

    else if (elt_micro == "HEX20") then
      
       K_loc_20 = Hexe20()

       nn_elems = pb%leaves(5)%dat_no / 20
    
       Do ii = 1, nn_elems

          kk = 1
          
          !** Translate Topology to global dofs *************************
          Do jj = (ii-1)*20+1 , ii*20
             
             id = pb%leaves(5)%p_int8(jj)
             id = (id - 1) * 3
             
             idxm_20(kk)   = (id    ) !*cr
             idxm_20(kk+1) = (id + 1) !*cr
             idxm_20(kk+2) = (id + 2) !*cr
             
             kk = kk + 3
             
          End Do
          
          idxn_20 = idxm_20

          Call MatSetValues(AA, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
          Call MatSetValues(AA_org, 60_8, idxm_20, 60_8 ,idxn_20, K_loc_20, ADD_VALUES, petsc_ierr)
          
       End Do
    
    end if


    IF (rank_mpi == 0) THEN   ! Sub Comm Master
        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- MatAssemblyBegin "//TRIM(nn_char)
        CASE default
            timer_name = "MatAssemblyBegin"
        End SELECT
        
        CALL start_timer(TRIM(timer_name), .FALSE.)
    END IF 

    Call MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    Call MatAssemblyBegin(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    ! Computations can be done while messages are in transition
    Call MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    Call MatAssemblyEnd(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)

    !------------------------------------------------------------------------------
    ! End timer
    !------------------------------------------------------------------------------
    IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

    !------------------------------------------------------------------------------
    ! Crashes on hawk. Out of memory?
    !------------------------------------------------------------------------------
    !  If (out_amount == "DEBUG") THEN 
    !     Call PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
    !     Call PetscViewerASCIIOpen(COMM_MPI,"AA.output.1",PetscViewer, petsc_ierr);
    !     Call PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
    !     Call MatView(AA, PetscViewer, petsc_ierr)
    !     Call PetscViewerDestroy(PetscViewer, petsc_ierr)
    !  End If
        

    !***************************************************************************
    !** At this point the system matrix is assembled. To make it ready to be ***
    !** used, the rows and columns of the dofs with prescribed displacements ***
    !** have to be eliminated. To do that with MatZeroRowsColumns we need    ***
    !** right hand side vectors and a solution vector.                       ***
    !***************************************************************************
    IF (rank_mpi == 0) THEN   ! Sub Comm Master
        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- create_rh_solution_vetors "//TRIM(nn_char)
        CASE default
            timer_name = "create_rh_solution_vetors"
        End SELECT
        
        CALL start_timer(TRIM(timer_name), .FALSE.)
    END IF 

    Do ii = 1, 24
    
       !** Create load vectors ***************************************
       call VecCreate(COMM_MPI, FF(ii), petsc_ierr)
       call VecSetSizes(FF(ii), PETSC_DECIDE, m_size, petsc_ierr)
       call VecSetFromOptions(FF(ii), petsc_ierr)
       Call VecSet(FF(ii), 0._rk,petsc_ierr)
       
       call VecGetOwnershipRange(FF(ii), IVstart, IVend, petsc_ierr)
 
       ! DEBUG INFORMATION
       If (out_amount == "DEBUG") THEN 
          Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
               "MPI rank : ",rank_mpi,"| vector ownership for ","F", ":",IVstart," -- " , IVend
       End If
       ! DEBUG INFORMATION       
       
       Call VecAssemblyBegin(FF(ii), petsc_ierr)
       !** Computations can be done while messages are in transition

    End Do

    Do ii = 1, 24 
       Call VecAssemblyEnd(FF(ii), petsc_ierr)
    End Do

    !** Get Bounds branch of LC 1 ****************************
    !** bb%leaves(1)%p_int8  : Boundary displacement node global ids
    !** bb%leaves(2)%p_real8 : Boundary displacement values
    write(desc,'(A,I0)')"Boundaries_"//trim(nn_char)//"_",1
    Call search_branch(trim(desc), pb, bb, success, DEBUG)

    !** Setup id reference vector for displacement insertion *********
    !** bb%leaves(2)%dat_no = Number of constraint nodes * 3
    Allocate(gnid_cref(bb%leaves(2)%dat_no), zeros_R8(bb%leaves(2)%dat_no))
    zeros_R8 = 0._rk
    
    kk = 1
    Do ii = 1, bb%leaves(1)%dat_no
       id = bb%leaves(1)%p_int8(ii)
       id = (id - 1)  * 3
       gnid_cref(kk:kk+2) = [id, id+1, id+2]
       kk= kk + 3
    End Do
    
    !** Create solution vector ***************************************
    call VecCreate(COMM_MPI, XX, petsc_ierr)
    call VecSetSizes(XX, PETSC_DECIDE, m_size, petsc_ierr)
    call VecSetFromOptions(XX, petsc_ierr)
    Call VecSet(XX, 0._rk,petsc_ierr)

    !------------------------------------------------------------------------------
    ! Only set prescribed displacements if there are boundary nodes available.
    ! These are calculated in struct_preprocess subroutine generate_boundaries
    !------------------------------------------------------------------------------
    IF(bb%leaves(2)%dat_no /= 0_ik) THEN
        !** Set prescribed displacements of LC1 to solution vector *******
        Call VecSetValues(XX, bb%leaves(2)%dat_no, &
            gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
    END IF 

    Call VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition
    Call VecAssemblyEnd(XX, petsc_ierr)


    call VecGetOwnershipRange(XX, IVstart, IVend, petsc_ierr)
 
    !------------------------------------------------------------------------------
    ! End timer
    !------------------------------------------------------------------------------
    IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))


    If (out_amount == "DEBUG") THEN 
        Write(un_lf,"('MM ', A,I4,A,A6,2(A,I9))")&
            "MPI rank : ",rank_mpi,"| vector ownership for ","X", ":",IVstart," -- " , IVend
        Call PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
        Call PetscViewerASCIIOpen(COMM_MPI,"FX.output.1",PetscViewer, petsc_ierr);
        Call PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
        Call VecView(XX, PetscViewer, petsc_ierr)
        Call PetscViewerDestroy(PetscViewer, petsc_ierr)
    End If
 
    !***************************************************************************
    !** At this point the right hand side vectors, filled with zeros and     ***
    !** the solution vector filled with the dirichlet boundary values of     ***
    !** load case 1 (LC1) are ready to be used.                              ***
    !** Now the right hand side vectors have to be modified with the         ***
    !** prescribed displacements. This is currently done by multiplying A by ***
    !** X filled with the negative displacement values.                      ***
    !** (See above VecSetValues(XX, ... , -bb%leaves(2)%p_real8, ... )       ***
    !***************************************************************************
    IF (rank_mpi == 0) THEN   ! Sub Comm Master
        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- compute_bndry_conditions "//TRIM(nn_char)
        CASE default
            timer_name = "compute_bndry_conditions"
        End SELECT
        
        CALL start_timer(TRIM(timer_name), .FALSE.)
    END IF 

    !** Compute dirichlet boundary corrections of first right  *******
    !** hand side vector.                                      *******
    Call MatMult(AA,XX,FF(1), petsc_ierr);

    IF(bb%leaves(2)%dat_no /= 0_ik) THEN
        !** Set zero values for dofs with prescribed displacements *******
        Call VecSetValues(FF(1), bb%leaves(2)%dat_no, &
            gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
    END IF 

    Call VecAssemblyBegin(FF(1), petsc_ierr)
    Call VecAssemblyEnd(FF(1), petsc_ierr)

    !** Compute dirichlet boundary corrections of 2nd to 23rd      ***
    !** right hand side vectors. The first on is already done and  ***
    !** the 24th will be done afterwards when the columns and rows ***
    !** of A set to zero.                                          ***
    Do ii = 2, 23

       !** Get Bounds branch of LC ii **********************************
       !** bb%leaves(1)%p_int8  : Boundary displacement node global ids
       !** bb%leaves(2)%p_real8 : Boundary displacement values
       write(desc,'(A,I0)')"Boundaries_"//trim(nn_char)//"_",ii
       Call search_branch(trim(desc), pb, bb, success, DEBUG)

      IF(bb%leaves(2)%dat_no /= 0_ik) THEN
         !** Set prescribed displacements of LCii to solution vector *****
         Call VecSetValues(XX, bb%leaves(2)%dat_no, &
               gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
      END IF 

       Call VecAssemblyBegin(XX, petsc_ierr)
       ! Computations can be done while messages are in transition
       Call VecAssemblyEnd(XX, petsc_ierr)

       !** Compute dirichlet boundary corrections of ii th right *******
       !** hand side vector.                                     *******
       Call MatMult(AA,XX,FF(ii), petsc_ierr);

      IF(bb%leaves(2)%dat_no /= 0_ik) THEN
         !** Set zero values for dofs with prescribed displacements ******
         Call VecSetValues(FF(ii), bb%leaves(2)%dat_no, &
               gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
      END IF 

       Call VecAssemblyBegin(FF(ii), petsc_ierr)
       
    End Do

    !** Get Bounds branch of LC 24 **********************************
    !** bb%leaves(1)%p_int8  : Boundary displacement node global ids
    !** bb%leaves(2)%p_real8 : Boundary displacement values
    write(desc,'(A,I0)')"Boundaries_"//trim(nn_char)//"_",24
    Call search_branch(trim(desc), pb, bb, success, DEBUG)
    
    !------------------------------------------------------------------------------
    ! Only set prescribed displacements if there are boundary nodes available.
    ! These are calculated in struct_preprocess subroutine generate_boundaries
    !------------------------------------------------------------------------------
    IF(bb%leaves(2)%dat_no /= 0_ik) THEN
       !** Set prescribed displacements of LC24 to solution vector *****
       Call VecSetValues(XX, bb%leaves(2)%dat_no, &
          gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
    END IF 
 
    Call VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition
 
    Call VecAssemblyEnd(XX, petsc_ierr)

   
    !** Finalize the open assembleys *************
    Do ii = 2, 23
       Call VecAssemblyEnd(FF(ii), petsc_ierr)
    End Do


    !** Since we are filling XX with the prescribed displacements  ***
    !** times -1. , we have to rescale XX before using it in       ***
    !** MatZeroRowsColumns.                                        ***
    Call VecScale(XX, -1._rk, petsc_ierr)
    
    !** Apply Dirichlet Boundaries to A and the 24th right hand    ***
    !** side vector                                                ***
    call MatZeroRowsColumns(AA, bb%leaves(2)%dat_no, gnid_cref, &
         0.0_8, XX, FF(24), petsc_ierr)
    
    !------------------------------------------------------------------------------
    ! End timer
    !------------------------------------------------------------------------------
    IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))


    ! DEBUG INFORMATION
    If (out_amount == "DEBUG") THEN 
       Call PetscViewerCreate(COMM_MPI, PetscViewer, petsc_ierr)
       Call PetscViewerASCIIOpen(COMM_MPI,"FF.output.1", PetscViewer, petsc_ierr);
       Call PetscViewerSetFormat(PetscViewer, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
       Call VecView(FF( 1), PetscViewer, petsc_ierr)
       Call PetscViewerDestroy(PetscViewer, petsc_ierr)
    End If
   
    !***************************************************************************
    !** At this point the right hand side vectors are modified with the       **
    !** prescribed displacements and the rows and columns of A representing   **
    !** dofs with prescribed displacements are filled with zeros.             **
    !** Now a linear solver context can be set up and the linear systems with **
    !** constant operator A and variable right hand side can be solved        **
    !***************************************************************************
    IF (rank_mpi == 0) THEN   ! Sub Comm Master
        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- solve_system "//TRIM(nn_char)
        CASE default
            timer_name = "solve_system"
        End SELECT
        
        CALL start_timer(TRIM(timer_name), .FALSE.)
    END IF 

    !** Create linear solver context *********************************
    CALL KSPCreate(COMM_MPI, ksp, petsc_ierr)
    
    ! Set operators. Here the matrix that defines the linear system
    ! also serves as the preconditioning matrix.
    CALL KSPSetOperators(ksp, AA, AA, petsc_ierr)
    
    ! Set linear solver defaults for this problem (optional).
    ! - By extracting the KSP and PC contexts from the KSP context,
    !   we can then directly call any KSP and PC routines to set
    !   various options.
    ! - The following two statements are optional; all of these
    !   parameters could alternatively be specified at runtime via
    !   KSPSetFromOptions().  All of these defaults can be
    !   overridden at runtime, as indicated below.
    Call KSPSetTolerances(ksp,1.e-2/((m_size+1)*(m_size+1)),1.e-50_rk, &
         PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, petsc_ierr);
    
    ! Set runtime options, e.g.,
    !     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    ! These options will override those specified above as long as
    ! KSPSetFromOptions() is called _after_ any other customization
    ! routines.
    CALL KSPSetFromOptions(ksp, petsc_ierr)
    
    ! DEBUG INFORMATION
    If (out_amount == "DEBUG") THEN 
       
       if ( rank_mpi == 0 ) then
          filename=''
          write(filename,'(A,I0,A)')trim(job_dir)//trim(project_name)//"_",nn,"_usg.vtk"
          Call write_data_head(filename, m_size/3)
       End if

    End If
    ! DEBUG INFORMATION
    
    !*****************************************************************
    !**                  Solve the linear system                   ***
    !*****************************************************************

!!!!!!!!!!!!!>> Development >>>>>>>>>>>>>>>>>>>>
!!!!!!!!!!!!!>> Development >>>>>>>>>>>>>>>>>>>>    
!!!!!!!!!!!!!>> Development >>>>>>>>>>>>>>>>>>>>

    !** Add a results branch to the mesh branch on rank 0 ************
    !** !! Initial try ... Communicate all results to master and !! **
    !** !! do serial calc_effective_material_parameters          !! **
    If (rank_mpi == 0) then

       call add_branch_to_branch(mb,resb)
       !** Raise branch with 4 child branches ************************
       !** 1. Displacements
       !** 2. Forces
       !** 3. Strains
       !** 4. Stresses
       Call raise_branch("Results of domain "//nn_char, 4, 0, resb)
       Call raise_branch("Displacements", 0, 0, resb%branches(1))
       Call raise_branch("Forces"       , 0, 0, resb%branches(2))
       Call raise_branch("Strains"      , 0, 0, resb%branches(3))
       Call raise_branch("Stresses"     , 0, 0, resb%branches(4))

       call log_tree(mb,un_lf,.FALSE.)
       
       !** Look again for the Part branch since the pb pointer  ******
       !** gets invalidated by dealloc of the branches array in ******
       !** add_branch_to_branch                                 ******
       domain_desc=''
       Write(part_desc,'(A,I0)')'Part_',parts_per_subdomain
       Call search_branch(trim(part_desc), mb, pb, success, DEBUG)

       ! Allocate global displacement result **********************
       Allocate(glob_displ(0:m_size-1))
       ! Allocate global forces result ****************************
       Allocate(glob_force(0:m_size-1))
       ! Allocate local bounds for global result ******************
       Allocate(res_sizes(2,size_mpi-1))
       
    End If

    Do jj = 1,24
       
       Call KSPSolve(ksp, FF(jj), XX, petsc_ierr)
       
       !** Get Bounds branch of LC jj **********************************
       !** bb%leaves(1)%p_int8  : Boundary displacement node global ids
       !** bb%leaves(2)%p_real8 : Boundary displacement values
       write(desc,'(A,I0)')"Boundaries_"//trim(nn_char)//"_",jj
       Call search_branch(trim(desc), pb, bb, success, DEBUG)
       
      !------------------------------------------------------------------------------
      ! Only set prescribed displacements if there are boundary nodes available.
      ! These are calculated in struct_preprocess subroutine generate_boundaries
      !------------------------------------------------------------------------------
      IF(bb%leaves(2)%dat_no /= 0_ik) THEN
         !** Extend XX with Boundary displacements *********************
         Call VecSetValues(XX, bb%leaves(2)%dat_no, &
            gnid_cref, bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
      END IF 
       
       Call VecAssemblyBegin(XX, petsc_ierr)
       ! Computations can be done while messages are in transition
       Call VecAssemblyEnd(XX, petsc_ierr)

       !** Calc reaction forces **************************************
       Call MatMult(AA_org, XX, FF(jj), petsc_ierr);
       
       !** Get Pointer to result vector **********
       Call VecGetArrayReadF90(XX,displ,petsc_ierr)
       
       !** Get Pointer to force vector ***********
       Call VecGetArrayReadF90(FF(jj),force,petsc_ierr)
       
       if ( rank_mpi > 0 ) then
          Call mpi_send([IVstart,IVend], 2_mik, &
               MPI_Integer8,  0_mik, rank_mpi, COMM_MPI, ierr)
          
          Call mpi_send(displ, Int(IVend-IVstart,mik), &
               MPI_Real8,  0_mik, Int(rank_mpi+size_mpi,mik), &
               COMM_MPI, ierr)

          Call mpi_send(force, Int(IVend-IVstart,mik), &
               MPI_Real8,  0_mik, Int(rank_mpi+2*size_mpi,mik), &
               COMM_MPI, ierr)
       Else
          
          ! Copy rank 0 local result *********************************
          glob_displ(IVstart:IVend-1) = displ
          glob_force(IVstart:IVend-1) = force
          
          ! Recv bounds of local results *****************************
          Do ii = 1, size_mpi-1
             Call mpi_recv(res_sizes(:,ii), 2_mik, &
                  MPI_Integer8,  Int(ii, mik), Int(ii, mik), &
                  COMM_MPI, status_mpi, ierr)
          End Do
          
          ! Recv local parts_per_subdomain of global result ************************
          Do ii = 1, size_mpi-1
             Call mpi_recv(glob_displ(res_sizes(1,ii):res_sizes(2,ii)-1), &
                  Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                  MPI_Integer8,  Int(ii, mik), Int(ii+size_mpi, mik), &
                  COMM_MPI, status_mpi, ierr)

             Call mpi_recv(glob_force(res_sizes(1,ii):res_sizes(2,ii)-1), &
                  Int(res_sizes(2,ii)-res_sizes(1,ii),mik), &
                  MPI_Integer8,  Int(ii, mik), Int(ii+2*size_mpi, mik) , &
                  COMM_MPI, status_mpi, ierr)
          End Do

          !** Add leaf with displacements to the results branch ******
          write(desc,'(A)')"Displacements"
          call Add_Leaf_to_Branch(resb%branches(1), trim(desc), &
                                  m_size,           glob_displ) 

          !** Add leaf with resultant forces to the results branch ***
          write(desc,'(A)')"Reaction Forces"
          call Add_Leaf_to_Branch(resb%branches(2), trim(desc), &
                                  m_size,           glob_force) 

          If (out_amount == "DEBUG") THEN 
          
             write(desc,'(A,I2.2)')"DispRes",jj
             Call write_vtk_data_Real8_vector_1D(&
                  matrix = reshape(glob_displ,[3,size(glob_displ)/3]), &
                  filename = trim(filename),  &
                  desc = trim(desc), head = .FALSE.,location="POINT_DATA")

             write(desc,'(A,I2.2)')"ForcRes",jj
             Call write_vtk_data_Real8_vector_1D(&
                  matrix = reshape(glob_force,[3,size(glob_force)/3]), &
                  filename = trim(filename),  &
                  desc = trim(desc), head = .FALSE.,location="POINT_DATA")
             
          End if
                 
       End If
    End Do

    !------------------------------------------------------------------------------
    ! End timer
    !------------------------------------------------------------------------------
    IF (rank_mpi == 0) CALL end_timer(TRIM(timer_name))

    !*****************************************************************
    !** All 24 linear system solutions are produced. Effective    ****
    !** stiffnesses can be calculated                             ****
    !*****************************************************************
    if ( rank_mpi == 0 ) then

        Deallocate(glob_displ, res_sizes, glob_force)

        SELECT CASE (timer_level)
        CASE (1)
            timer_name = "+-- calc_eff_stiffness "//TRIM(nn_char)
        CASE default
            timer_name = "calc_eff_stiffness"
        End SELECT

        Call start_timer(trim(timer_name), .FALSE.)
        call calc_effective_material_parameters(root, nn)
        Call end_timer(trim(timer_name))

    End if
    
!!!!!!!!!!!!!<< Development <<<<<<<<<<<<<<<<<<<<
    
    call MatDestroy(AA, petsc_ierr)
    Call VecDestroy(XX, petsc_ierr)
    Do ii = 1, 24
       Call VecDestroy(FF(ii) , petsc_ierr)
    End Do
    
    if ( rank_mpi > 0 ) then
       no_data = 0
       call destroy_tree(pb,no_data)
       deallocate(pb)
    End if

  End Subroutine exec_single_domain

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
 
  !-- MPI Variables -------------------------------------------------------------------------------
  INTEGER(KIND=mik) :: ierr, rank_mpi, size_mpi
  INTEGER(KIND=mik) :: petsc_ierr
  INTEGER(KIND=mik) :: worker_rank_mpi, worker_size_mpi
  INTEGER(KIND=mik) :: Active, request, finished, worker_comm
  INTEGER(KIND=mik), Dimension(no_streams) :: fh_mpi

  INTEGER(KIND=mik), Dimension(MPI_STATUS_SIZE)  :: status_mpi
  INTEGER(KIND=mik), Dimension(:,:), Allocatable :: statuses_mpi
  INTEGER(KIND=mik), Dimension(:)  , Allocatable :: Activity, req_list

  !------------------------------------------------------------------------------------------------
  CHARACTER(Len=mcl)  :: link_name = 'struct process'
  INTEGER(kind=c_int) :: stat_c_int
  INTEGER             :: stat
  Type(tBranch)       :: root, phi_tree
  Type(tBranch), pointer :: ddc, meta_para, res
  
  CHARACTER, DIMENSION(4*mcl) :: c_char_array
  CHARACTER, DIMENSION(:), ALLOCATABLE :: char_arr
  CHARACTER(LEN=4*mcl), DIMENSION(:), ALLOCATABLE :: domain_path
  CHARACTER(LEN=mcl)  , DIMENSION(:), ALLOCATABLE :: m_rry      
  CHARACTER(LEN=4*mcl) :: job_dir
  CHARACTER(LEN=1) :: restart='N', restart_cmd_arg='U' ! U = 'undefined'
  CHARACTER(LEN=mcl) :: cmd_arg_history=''
  CHARACTER(LEN=mcl) :: muCT_pd_path, muCT_pd_name, domain_desc, binary
  CHARACTER(LEN=8) :: elt_micro, output
  
  
  REAL(KIND=rk), DIMENSION(3) :: delta
  REAL(KIND=rk) :: strain
  
  INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE   :: Domains, Domain_stats, act_domains, nn_D
  INTEGER(KIND=ik), DIMENSION(3) :: xa_d, xe_d, vdim
  
  INTEGER(KIND=ik) :: nn, ii, jj, kk, dc
  INTEGER(KIND=ik) :: amount_domains, path_count
  INTEGER(KIND=ik) :: alloc_stat, aun, free_file_handle
  INTEGER(KIND=ik) :: Domain, llimit, parts_per_subdomain, elo_macro
  
  INTEGER(KIND=pd_ik), DIMENSION(:), ALLOCATABLE :: serial_root
  INTEGER(KIND=pd_ik), DIMENSION(no_streams) :: dsize
  INTEGER(KIND=pd_ik) :: serial_root_size
  LOGICAL :: success, fexist, heaxist, stp = .FALSE.

  !----------------------------------------------------------------------------
 
  CALL mpi_init(ierr)
  CALL print_err_stop(std_out, "MPI_INIT didn't succeed", INT(ierr, KIND=ik))

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)
  CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", INT(ierr, KIND=ik))
 
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
  CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't be retrieved", INT(ierr, KIND=ik))
 
  If (size_mpi < 2) CALL print_err_stop(std_out, "We need at least 2 MPI processes to execute this program.", 1)
  
  !------------------------------------------------------------------------------
  ! Rank 0 -- Init (Master) Process and broadcast init parameters 
  !------------------------------------------------------------------------------
  If (rank_mpi==0) Then
 
      Call Start_Timer("Init Process")

      !------------------------------------------------------------------------------
      ! Parse the command arguments
      !------------------------------------------------------------------------------
      CALL get_cmd_args(binary, in%full, stp, restart_cmd_arg, cmd_arg_history)
      IF(stp) GOTO 1001
      
      IF (in%full=='') THEN
         CALL usage(binary)    

         !------------------------------------------------------------------------------
         ! On std_out since file of std_out is not spawned
         !------------------------------------------------------------------------------
         CALL print_err_stop(6, "No input file given", 1)
      END IF
    
      !------------------------------------------------------------------------------
      ! Check and open the input file; Modify the Meta-Filename / Basename
      ! Define the new application name first
      !------------------------------------------------------------------------------
      global_meta_prgrm_mstr_app = 'ddtc' 
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

      CALL show_title()
   
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
      CALL meta_read('SIZE_DOMAIN'      , m_rry, bone%phdsize)
      CALL meta_read('SPACING'          , m_rry, bone%delta)
      CALL meta_read('DIMENSIONS'       , m_rry, vdim)
      CALL meta_read('LO_BNDS_DMN_RANGE', m_rry, xa_d)
      CALL meta_read('UP_BNDS_DMN_RANGE', m_rry, xe_d)
      CALL meta_read('BINARIZE_LO'      , m_rry, llimit)
      CALL meta_read('MESH_PER_SUB_DMN' , m_rry, parts_per_subdomain)
      CALL meta_read('RVE_STRAIN'       , m_rry, strain)
      CALL meta_read('YOUNG_MODULUS'    , m_rry, bone%E)
      CALL meta_read('POISSON_RATIO'    , m_rry, bone%nu)
      CALL meta_read('MACRO_ELMNT_ORDER', m_rry, elo_macro)

      !------------------------------------------------------------------------------
      ! Restart handling
      ! Done after meta_io to decide based on keywords
      !------------------------------------------------------------------------------
      CALL meta_handle_lock_file(restart, restart_cmd_arg)

      !------------------------------------------------------------------------------
      ! Spawn a log file and a results file
      !------------------------------------------------------------------------------
      ! This log file may collide with the original log file (!)
      ! The regular struct_process log file contains still has the "old" basename!
      !------------------------------------------------------------------------------
      ! CALL meta_start_ascii(fhl, log_suf)
      CALL meta_start_ascii(fhmon, mon_suf)
      
      IF (std_out/=6) THEN
         CALL meta_start_ascii(std_out, '.std_out')

         CALL show_title() 
      END IF

      CALL meta_write('DBG_LVL', out_amount )

      !------------------------------------------------------------------------------
      ! Warning / Error handling
      !------------------------------------------------------------------------------
      IF ( (bone%phdsize(1) /= bone%phdsize(2)) .OR. (bone%phdsize(1) /= bone%phdsize(3)) ) THEN
         CALL print_err_stop(std_out, 'Currently, all 3 dimensions of the physical domain size must be equal!', 1)
      END IF
      
      IF ( (delta(1) /= delta(2)) .OR. (delta(1) /= delta(3)) ) THEN
         CALL print_err_stop(std_out, 'Currently, the spacings of all 3 dimensions must be equal!', 1)
      END IF

      IF ( (xa_d(1) > xe_d(1)) .OR. (xa_d(2) > xe_d(2)) .or. (xa_d(3) > xe_d(3)) ) THEN
         CALL print_err_stop(std_out, 'Input parameter error: Start value of domain range larger than end value.', 1)
      END IF

      ! Program breaks if the phdsize is not taking the boundary nodes into account (!).
      ! Therefore, the boundaries are calculated with + 2 Voxels
      IF((bone%phdsize(1) > (vdim(1) + 2_ik)*bone%delta(1)) .OR. & 
         (bone%phdsize(2) > (vdim(2) + 2_ik)*bone%delta(2)) .OR. & 
         (bone%phdsize(3) > (vdim(3) + 2_ik)*bone%delta(3))) THEN
         CALL print_err_stop(std_out, &
            'The domains are larger than the field of view.', 1)
      END IF
      
      !------------------------------------------------------------------------------
      ! Each subdomain gets computed by a user defined amount of processors. This 
      ! amount of processors equals to a specific amount of mesh parts_per_subdomain.
      ! Ideally, all processors are used. Therefore, MOD(size_mpi-1, parts_per_subdomain) shall 
      ! resolve without a remainder. "-1" to take the master process into account.
      !------------------------------------------------------------------------------
      IF (MOD(size_mpi-1, parts_per_subdomain) /= 0) THEN
         CALL print_err_stop(std_out, 'mod(size_mpi-1,parts_per_subdomain) /= 0 ! This case is not supported.', 1)
      END IF

      !------------------------------------------------------------------------------
      ! Raise and build meta_para tree
      ! Hardcoded, implicitly given order of the leafs. 
      ! DO NOT CHANGE ORDER WITHOUT MODIFYING ALL OTHER INDICES REGARDING »meta_para«
      !------------------------------------------------------------------------------
      Allocate(meta_para)
      Call raise_tree("Input parameters", meta_para)

      CALL add_leaf_to_branch(meta_para, "muCT puredat pro_path"                , mcl , str_to_char(muCT_pd_path))
      CALL add_leaf_to_branch(meta_para, "muCT puredat pro_name"                , mcl , str_to_char(muCT_pd_name))
      CALL add_leaf_to_branch(meta_para, "Physical domain size"                 , 3_ik, bone%phdsize)
      CALL add_leaf_to_branch(meta_para, "Lower bounds of selected domain range", 3_ik, xa_d)
      CALL add_leaf_to_branch(meta_para, "Upper bounds of selected domain range", 3_ik, xe_d)     
      CALL add_leaf_to_branch(meta_para, "Grid spacings"                        , 3_rk, bone%delta)
      CALL add_leaf_to_branch(meta_para, "Lower limit of iso value"      , 1_ik, [llimit])     
      CALL add_leaf_to_branch(meta_para, "Element type  on micro scale"  , len(elt_micro) , str_to_char(elt_micro))     
      CALL add_leaf_to_branch(meta_para, "No of mesh parts per subdomain", 1_ik           , [parts_per_subdomain])
      CALL add_leaf_to_branch(meta_para, "Output Format"                 , len(output)    , str_to_char(output))
      CALL add_leaf_to_branch(meta_para, "Average strain on RVE"         , 1_ik           , [strain])   
      CALL add_leaf_to_branch(meta_para, "Young_s modulus"               , 1_ik           , [bone%E])
      CALL add_leaf_to_branch(meta_para, "Poisson_s ratio"               , 1_ik           , [bone%nu])
      CALL add_leaf_to_branch(meta_para, "Element order on macro scale"  , 1_ik           , [elo_macro])
      CALL add_leaf_to_branch(meta_para, "Output amount"                 , len(out_amount), str_to_char(out_amount))
      CALL add_leaf_to_branch(meta_para, "Restart"                       , 1_ik           , str_to_char(restart))
      CALL add_leaf_to_branch(meta_para, "Number of voxels per direction", 3_ik           , vdim)

      !------------------------------------------------------------------------------
      ! Prepare output directory via calling the c function.
      ! Required, because INQUIRE only acts on files, not on directories.
      ! File exists if stat_c_int = 0 
      !------------------------------------------------------------------------------
      c_char_array(1:LEN(TRIM(outpath)//CHAR(0))) = str_to_char(TRIM(outpath)//CHAR(0))
      CALL Stat_Dir(c_char_array, stat_c_int)

      IF(stat_c_int /= 0) THEN

         CALL execute_command_line("mkdir -p "//TRIM(outpath),CMDSTAT=stat)

         IF(stat /= 0) THEN
            CALL print_err_stop(std_out, 'Could not execute syscall »mkdir -p '//trim(outpath)//'«.', 1)
         END IF 

         CALL Stat_Dir(c_char_array, stat_c_int)

         IF(stat_c_int /= 0) THEN
            CALL print_err_stop(std_out, 'Could not create the output directory »'//TRIM(outpath)//'«.', 1)
         END IF
      ELSE 
         WRITE(un_mon, FMT_MSG) "Reusing the output directory"
         WRITE(un_mon, FMT_MSG) TRIM(outpath)
      END IF

      CALL link_start(link_name, .TRUE., .FALSE., success)
      IF (.NOT. success) CALL print_err_stop(std_out, "Something went wrong during link_start", 1)
   

      !------------------------------------------------------------------------------
      ! Allocate and init field for selected domain range
      !------------------------------------------------------------------------------
      amount_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

      Allocate(Domains(amount_domains),stat=alloc_stat)
      Call alloc_err("Domains",alloc_stat)

      Allocate(Domain_stats(amount_domains),stat=alloc_stat)
      Call alloc_err("Domain_stats",alloc_stat)

      Allocate(domain_path(0:amount_domains))
      domain_path = ''


      !------------------------------------------------------------------------------
      ! New activity tracker unit
      !------------------------------------------------------------------------------
      aun = give_new_unit()

      !------------------------------------------------------------------------------
      ! Check whether there already is a project header
      !------------------------------------------------------------------------------
      INQUIRE(FILE=TRIM(pro_path)//TRIM(pro_name)//'.head', EXIST=heaxist)

      !------------------------------------------------------------------------------
      ! The Output name normally is different than the input name.
      ! Therefore, an existing header implies a restart.
      !------------------------------------------------------------------------------
      IF (restart == 'N') THEN ! Project header not available

         IF(heaxist) CALL print_err_stop(std_out, &
            "Restart requested, but header already exists", 1)

         !------------------------------------------------------------------------------
         ! Create a PureDat header file, based on the meta file and the raw bin blob
         !------------------------------------------------------------------------------
         free_file_handle = give_new_unit()
         CALL convert_meta_to_puredat(free_file_handle, m_rry)

         !------------------------------------------------------------------------------
         ! project_name --> out%p_n_bsnm/bsnm --> subdirectory with file name = bsnm.suf
         !------------------------------------------------------------------------------
         ! Tree which is fed back = root. A collection of stream paths and Null pointers
         !------------------------------------------------------------------------------
         Call raise_tree(Trim(project_name),root)

         !------------------------------------------------------------------------------
         ! Source branch / target branch
         !------------------------------------------------------------------------------
         Call include_branch_into_branch(s_b=meta_para, t_b=root, blind=.TRUE.)

         !------------------------------------------------------------------------------
         ! Read an existing input tree (with microfocus ct data).
         !------------------------------------------------------------------------------
         !** Load puredat tree of micro-CT data and calculate the global
         !** parameters of the domain decomposition
         pro_path = muCT_pd_path
         pro_name = muCT_pd_name

         phi_tree = read_tree()

         ! !** Set project name and path of global domain decomposition     
         pro_path = outpath
         pro_name = project_name


         allocate(ddc)
         ddc = calc_general_ddc_params(bone%phdsize, phi_tree)
         
         call include_branch_into_branch(s_b=ddc, t_b=root, blind=.TRUE.)

         !------------------------------------------------------------------------------
         ! Initialize the activity tracker.
         !------------------------------------------------------------------------------
         OPEN(aun, FILE=TRIM(outpath)//"/"//trim(project_name)//".status", &
            ACTION="WRITE", STATUS="REPLACE", ACCESS="STREAM")

         Domain_stats = 1_ik

         WRITE(aun) Domain_stats
         FLUSH(aun)

      ELSE ! restart = 'Y'         

         IF(.NOT. heaxist) THEN  
            mssg="Restart requested, but header does not exist. &
            &Please specify 'RESTART = N' if there was no previous computation. &
            &Currently, the restart procedure does not support an automatic switch &
            &to 'RESTART = N'."
            CALL print_err_stop(std_out, mssg, 1)
         END IF

        !------------------------------------------------------------------------------
        ! Read an existing output tree (with microfocus ct data).
        !------------------------------------------------------------------------------
        root = read_tree()

        ! DEBUG INFORMATION
        If (out_amount == "DEBUG") THEN 
           WRITE(un_lf,FMT_DBG_SEP)
           Write(un_lf,'(A)')"root right after restart read"
           Call log_tree(root,un_lf,.FALSE.)
           WRITE(un_lf,FMT_DBG_SEP)
           flush(un_lf)
        END If
        ! DEBUG INFORMATION
       
        !------------------------------------------------------------------------------
        ! This calling sequence is only valid since "Averaged Material 
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
            CALL print_err_stop(std_out, mssg, 1)
         END IF

         CALL search_branch("Input parameters", root, meta_para, success)

         IF (.NOT. success) then
            mssg = "No branch named 'Input parameters', however a restart was requested."
            CALL print_err_stop(std_out, mssg, 1)
         END IF

         !------------------------------------------------------------------------------
         ! Reset Output amount and Restart in loaded param branch
         ! Hardcoded, implicitly given order of the leafs. 
         ! DO NOT CHANGE INDICES WITHOUT MODYFING THE »add_leaf_to_branch« SEQUENCES.
         !------------------------------------------------------------------------------
         meta_para%leaves(14)%p_char = str_to_char(out_amount)
         meta_para%leaves(15)%p_char = "Y"

         ! DEBUG INFORMATION
         If (out_amount == "DEBUG") THEN 
            Write(un_lf,fmt_dbg_sep)
            Write(un_lf,'(A)')"root right after restart and deletion of avg mat props branch"
            Call log_tree(root,un_lf,.FALSE.)
            Write(un_lf,fmt_dbg_sep)
            flush(un_lf)
         END If

         !------------------------------------------------------------------------------
         ! Read the activity tracker.
         !------------------------------------------------------------------------------
         INQUIRE(aun, EXIST=fexist)

         IF (.NOT. fexist) THEN
            mssg='The file '//TRIM(outpath)//"/"//trim(project_name)//"does not exist."//'.'
            CALL print_err_stop(std_out, TRIM(ADJUSTL(mssg)), 1_ik)
         END IF

         ! Open to read
         OPEN(aun, FILE=TRIM(outpath)//"/"//trim(project_name)//".status", &
            ACTION="READ", STATUS="OLD", ACCESS="STREAM")

         READ(aun) Domain_stats

         CLOSE(aun)

         ! Open to write
         OPEN(aun, FILE=TRIM(outpath)//"/"//trim(project_name)//".status", &
            ACTION="WRITE", STATUS="OLD", ACCESS="STREAM")

      END IF ! restart == Yes/No
         
      Call pd_get(ddc,"nn_D",nn_D)

     !** Init Result branch ***************************************************
     !if ( Restart == "N" ) Then
     
     Call add_branch_to_branch(root,res)
     Call raise_branch("Averaged Material Properties", 0_pd_ik, 18_pd_ik, res)
     
     Call raise_leaves(no_leaves = 18_pd_ik, &
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
          "Effective density                                 "] , &
          dat_ty = [(5_1,ii=1,18)], &
          dat_no = [ &
          amount_domains * 24*24, amount_domains * 24*24, amount_domains        , &
          amount_domains *  6*24, amount_domains *  6*24, amount_domains *  6* 6, &
          amount_domains        , amount_domains *  6* 6, amount_domains        , &
          amount_domains        , amount_domains *     3, amount_domains *     9, &
          amount_domains *  6* 6, amount_domains        , amount_domains *     3, &
          amount_domains *     9, amount_domains *  6* 6, amount_domains           ], &
          branch = res)
     
     res%leaves(:)%pstat = -1
     
     Call set_bounds_in_branch(root, root%streams)
     
     !End If
     
     ! DEBUG INFORMATION
     If (out_amount == "DEBUG") THEN 
        Write(un_lf,fmt_dbg_sep)
        Write(un_lf,'(A)')"root right before serialisation"
        Call log_tree(root,un_lf,.True.)
        Write(un_lf,fmt_dbg_sep)
        flush(un_lf)
     END If
     ! DEBUG INFORMATION

     !** serialize root branch ************************************************
     Call serialize_branch(root,serial_root,serial_root_size,.TRUE.)

     !** Init Domain Cross Reference and domain paths *************************
     dc = 0
     domain_path(0) = Trim(pro_path)//Trim(pro_name)//"_domains" ! "_results"

     nn = 1
     Do kk = xa_d(3), xe_d(3)
        Do jj = xa_d(2), xe_d(2)
           Do ii = xa_d(1), xe_d(1)
              
              dc         = dc + 1_mik
              path_count = dc / 2_mik
              Write(domain_path(dc),'(A,"/",I0)')Trim(domain_path(path_count)),dc

              Domains(nn) = ii + jj * nn_D(1)         + &
                                 kk * nn_D(1)*nn_D(2)
              nn = nn + 1_mik

           End Do
        End Do
     End Do

     !** Generate Activity_List ***********************************************
     Allocate(Activity(size_mpi-1), stat=alloc_stat)
     Call alloc_err("Activity_List", alloc_stat)

     Activity=1

     Allocate(act_domains(size_mpi-1),stat=alloc_stat)
     Call alloc_err("act_domains",alloc_stat)

     act_domains = 0
     
     If ( (amount_domains*parts_per_subdomain < size_mpi-1) .OR. &
          ( count( Domain_stats < 10 ) == 0 ) .OR. &
          ((count( Domain_stats < 10 )*parts_per_subdomain) < size_mpi-1) ) Then

        If (amount_domains*parts_per_subdomain < size_mpi-1) then
           Write(un_mon,FMT_ERR_AxI0)"amount_domains: ", amount_domains
           Write(un_mon,FMT_ERR_AxI0)"parts_per_subdomain: ", parts_per_subdomain
           Write(un_mon,FMT_ERR_AxI0)"size_mpi-1: ", size_mpi-1
           Write(un_mon,FMT_ERR)"amount_domains*parts_per_subdomain < size_mpi-1"

        Else If ((count( Domain_stats < 10 )*parts_per_subdomain) < size_mpi-1) then
           Write(un_mon, FMT_ERR)   "Remaining amount_domains < Number of Solution Master"
           Write(un_mon, FMT_ERR_xAI0) "Remaining amount_domains:   ", count( Domain_stats < 10 )
           Write(un_mon, FMT_ERR_xAI0) "Number of solution masters: ", (size_mpi-1)/parts_per_subdomain
        Else

           Write(un_mon, FMT_ERR)"Restart on fully finished job."
        End If
        Write(un_mon, FMT_ERR)"This case is not supported."
        
        Call mpi_bcast(pro_path, INT(mcl,mik), MPI_CHAR, 0_mik,&
                       MPI_COMM_WORLD, ierr)

        Call mpi_bcast(pro_name, INT(mcl,mik), MPI_CHAR, 0_mik,&
             MPI_COMM_WORLD, ierr)
        
        !** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
        Call mpi_bcast(-1_ik, 1_mik, MPI_INTEGER8, 0_mik,&
                       MPI_COMM_WORLD, ierr)

        Call End_Timer("Init Process")
 
        Goto 1000

     End If

     Call End_Timer("Init Process")

     !*************************************************************************
     !** Start Workers ********************************************************

     Call Start_Timer("Broadcast Init meta_para")
     
     Call mpi_bcast(pro_path, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

     Call mpi_bcast(pro_name, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

     write(un_lf, FMT_MSG_xAI0)"Broadcasting serialized root of size [Byte] ", serial_root_size*8
     
     Call mpi_bcast(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
     
     Call mpi_bcast(serial_root, INT(serial_root_size,mik), MPI_INTEGER8, 0_mik,&
          MPI_COMM_WORLD, ierr)

     Call End_Timer("Broadcast Init meta_para")

     !** Execute collective mpi_comm_split. Since mpi_comm_world rank 0 is
     !** the head master worker_comm is not needed and it should not be in
     !** any worker group and communicator. With MPI_UNDEFINED passed as
     !** color worker_comm gets the value MPI_COMM_NULL
     Call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, &
          rank_mpi, worker_comm, ierr)
     CALL print_err_stop(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", INT(ierr, KIND=ik))


  !****************************************************************************
  !** Ranks > 0 -- Worker slaves **********************************************
  !****************************************************************************
  Else
     
     !** Broadcast recieve init parameters ************************************
     Call Start_Timer("Broadcast Init meta_para")
     
     Call mpi_bcast(outpath         , INT(mcl,mik), MPI_CHAR    , 0_mik, MPI_COMM_WORLD, ierr)
     Call mpi_bcast(project_name    , INT(mcl,mik), MPI_CHAR    , 0_mik, MPI_COMM_WORLD, ierr)
     Call mpi_bcast(serial_root_size, 1_mik       , MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

     !** Serial_root_size == -1 ==> Signal that amount_domains < size_mpi-1 ****
     If ( serial_root_size == -1 ) then

        Call End_Timer("Broadcast Init meta_para")
        
        Goto 1001

     End If
     
     Allocate(serial_root(serial_root_size))
     
     Call mpi_bcast(serial_root, INT(serial_root_size,mik), MPI_INTEGER8, 0_mik,&
          MPI_COMM_WORLD, ierr)

     Call End_Timer("Broadcast Init meta_para")

     !** Deserialize root branch **********************************************
     Call Start_Timer("Deserialize root branch")

     Call deserialize_branch(root, serial_root, .TRUE.)
     Call assign_pd_root (root)
     Call set_bounds_in_branch(root,root%streams)

     ! Commented out since out_amount is a parametrized global variable 
     ! Call pd_get(root%branches(1),"Output amount", char_arr)
     ! out_amount = char_to_str(char_arr)
     ! deallocate(char_arr)

     Call pd_get(root%branches(1),"Restart", char_arr)
     restart = char_to_str(char_arr)
     deallocate(char_arr)

     Call End_Timer("Deserialize root branch")

     !** Init Domain Cross Reference ******************************************
     Call pd_get(root%branches(1), "Lower bounds of selected domain range", xa_d, 3)
     Call pd_get(root%branches(1), "Upper bounds of selected domain range", xe_d, 3)
     
     amount_domains= (xe_d(1)-xa_d(1)+1) * &
                     (xe_d(2)-xa_d(2)+1) * &
                     (xe_d(3)-xa_d(3)+1)
     
     Allocate(Domains(amount_domains),stat=alloc_stat)
     Call alloc_err("Domains",alloc_stat)

     Call pd_get(root%branches(2),"nn_D",nn_D)
     
     nn = 1
     Do kk = xa_d(3), xe_d(3)
        Do jj = xa_d(2), xe_d(2)
           Do ii = xa_d(1), xe_d(1)
              
              Domains(nn) = ii + jj * nn_D(1)         + &
                                 kk * nn_D(1)*nn_D(2)
              nn = nn + 1_mik

           End Do
        End Do
     End Do

     Call pd_get(root%branches(1),"No of mesh parts per subdomain",parts_per_subdomain)

     !*************************************************************************
     !** All Worker Ranks -- Init worker Communicators ************************
     !*************************************************************************
     Call MPI_Comm_split(MPI_COMM_WORLD, Int((rank_mpi-1)/parts_per_subdomain,mik), &
                         rank_mpi, worker_comm, ierr)
     CALL print_err_stop(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", INT(ierr, KIND=ik))
     
     Call MPI_COMM_RANK(WORKER_COMM, worker_rank_mpi, ierr)
     CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't retrieve worker_rank_mpi", INT(ierr, KIND=ik))

 
     Call MPI_COMM_SIZE(WORKER_COMM, worker_size_mpi, ierr)
     CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't retrieve worker_size_mpi", INT(ierr, KIND=ik))

     !** This sets the options for PETSc in-core. To alter the options ***
     !** add them in Set_PETSc_Options in Module pets_opt in file      ***
     !** f-src/mod_parameters.f90                                      ***
     Call Set_PETSc_Options()

     PETSC_COMM_WORLD = worker_comm
     CALL PetscInitialize(PETSC_NULL_CHARACTER, petsc_ierr)
  End If

  !****************************************************************************
  !** All Ranks -- Init MPI request and status lists **************************
  !****************************************************************************
  Allocate(req_list(size_mpi-1),stat=alloc_stat)
  Call alloc_err("req_list",alloc_stat)
  req_list=0
  
  Allocate(statuses_mpi(MPI_STATUS_SIZE,size_mpi-1),stat=alloc_stat)
  Call alloc_err("statuses_mpi",alloc_stat)
  
  !****************************************************************************
  !** All Ranks -- Init MPI IO-System *****************************************
  !****************************************************************************
  Call Open_Stream_Files(root%streams, "write", "replace", fh_mpi)

  !****************************************************************************
  !** Rank 0 -- Process master Start working process **************************
  !****************************************************************************
  If (rank_mpi==0) Then

     call get_stream_size(root, dsize)

     ! DEBUG INFORMATION
     If (out_amount == "DEBUG") THEN
        write(un_mon,FMT_MSG_AxI0)"On rank zero, stream sizes: ",dsize
     End If
     ! DEBUG INFORMATION
     
     nn = 1
     ii = 1

      !------------------------------------------------------------------------------
      ! Supply all worker masters  with their first work package
      ! ii is incremented by ii = ii + parts_per_subdomain
      !------------------------------------------------------------------------------
      Do While (ii <= (size_mpi-1_mik))

         if (nn > amount_domains) exit

         If ( Domain_stats(nn) /= 1 ) then
            nn = nn + 1_mik
            cycle
         End If

         act_domains(ii) = nn

         Do jj = ii, ii + parts_per_subdomain-1
            
            !** Activity = 1 (Set above during init) ***
            CALL mpi_send(Activity(jj), 1_mik, mpi_integer4, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

            CALL mpi_send(nn, 1_mik, mpi_integer8, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", INT(ierr, KIND=ik))
            
            if (out_amount /= "PRODUCTION") then
               Call mpi_send(domain_path(nn), Int(4_mik*mcl,mik), &
                     MPI_CHARACTER, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
               CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", INT(ierr, KIND=ik))
            End if
         End Do
         
         !** Log to global stdout **********************************************
         Write(un_mon, FMT_MSG_xAI0)"MPI rank: ",ii, "      Domain number: ",Domains(nn)
         flush(un_mon)

         nn = nn + 1_mik
         
         Call MPI_IRECV(Activity(ii), 1_mik, MPI_INTEGER4, Int(ii,mik), Int(ii,mik), &
               MPI_COMM_WORLD, REQ_LIST(ii), IERR)
         CALL print_err_stop(std_out, "MPI_IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

         ii = ii + Int(parts_per_subdomain,mik)

      End Do

      If ( restart == "N" ) then
         !** Add leaf with analyzed cube numbers to root ***********************
         Call add_leaf_to_branch(root, "Domain Numbers", amount_domains, Domains)
         Call set_bounds_in_branch(root, root%streams)
      End If

      !** Write Root header and input parameters *******************************
      Call Start_Timer("Write Root Branch")
      call store_parallel_branch(root, FH_MPI)
      Call Write_Tree(root)
      Call End_Timer("Write Root Branch")
         
      Call MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
      CALL print_err_stop(std_out, &
      "MPI_WAITANY on req_list for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

      ii = finished

      Domain_stats(act_domains(ii)) = Activity(ii)

      write(aun,pos=(act_domains(ii)-1)*ik+1) INT(Activity(ii), KIND=ik)
      flush(aun)
      act_domains(ii) = nn

      Do While (nn <= amount_domains)

         If ( Domain_stats(nn) /= 1 ) then
            if (nn > amount_domains) exit
            nn = nn + 1_mik
            cycle
         End If

         Do jj = ii, ii + parts_per_subdomain- 1

            Activity(jj) = 1_mik
            
            Call mpi_send(Activity(jj), 1_mik, mpi_integer4 , Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

            Call mpi_send(nn          , 1_mik, mpi_integer8, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", INT(ierr, KIND=ik))

            if (out_amount /= "PRODUCTION") then
               Call mpi_send(domain_path(nn), Int(4_mik*mcl,mik), &
                     MPI_CHARACTER, Int(jj,mik), Int(jj,mik), MPI_COMM_WORLD,ierr)
               CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", INT(ierr, KIND=ik))
            End if
         End Do
         
         !** Log to global stdout **********************************************
         Write(un_mon,'(2(A,I10))')"MPI rank : ",ii, " ; Domain number: ",Domains(nn)
         flush(un_mon)
         
         nn = nn + 1_mik
         
         Call MPI_IRECV(Activity(ii), 1_mik, MPI_INTEGER4, Int(ii,mik), &
                        Int(ii,mik), MPI_COMM_WORLD, REQ_LIST(ii), IERR)
         CALL print_err_stop(std_out, "MPI_IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

         Call MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
         CALL print_err_stop(std_out, "MPI_WAITANY on   for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

         ii = finished

         ! DEBUG INFORMATION
         If (out_amount == "DEBUG") THEN
            Write(un_mon,*)"Domain ",ii, Domains(act_domains(ii))," finished"
         End If
         ! DEBUG INFORMATION
         
         Domain_stats(act_domains(ii)) = Activity(ii)

         write(aun,pos=(act_domains(ii)-1)*ik+1) INT(Activity(ii), KIND=ik)
         flush(aun)
         act_domains(ii) = nn

      End Do

      !** Write last data element to ensure correct file size ******************
      Call MPI_FILE_WRITE_AT(FH_MPI(5), &
            Int(res%leaves(18)%lbound-1+amount_domains, MPI_OFFSET_KIND), &
            0_pd_rk, Int(1,pd_mik), MPI_Real8, status_mpi, ierr)

         Call MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
         CALL print_err_stop(std_out, &
         "MPI_WAITANY on req_list for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))


      !** TODO refactor domain_cross reference from size size_mpi-1 to *********
      !** (size_mpi-1)/parts_per_subdomain ***************************************************
      Do ii = 1, size_mpi-1, parts_per_subdomain
         write(aun,pos=(act_domains(ii)-1)*ik+1) INT(Activity(ii), KIND=ik)
      End Do

      flush(aun)


   Activity = -1
   
   Do ii = 1_mik, size_mpi-1_mik
      Call mpi_send(Activity(ii), 1_mik, mpi_integer4, Int(ii,mik), &
                     Int(ii,mik), MPI_COMM_WORLD,ierr)
      CALL print_err_stop(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))
   End Do

  !****************************************************************************
  !** Ranks > 0 -- Workers ****************************************************
  !****************************************************************************
  Else 

      !------------------------------------------------------------------------------
      ! Extend project_name and outpath with rank bone%number
      !------------------------------------------------------------------------------
      WRITE(outpath,'(A,A,I7.7,A)') TRIM(outpath),"Rank_",rank_mpi,"/"
      WRITE(project_name,'(A,A,I7.7)') TRIM(project_name),"_",rank_mpi
     
      !------------------------------------------------------------------------------
      ! Prepare output directory via calling the c function.
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
            CALL print_err_stop(std_out,'Could not create the output directory »'//TRIM(outpath)//'«.', 1)
         END IF
      ELSE 
         WRITE(un_mon, FMT_MSG) "Reusing the output directory "//TRIM(outpath)
      END IF

     Call link_start(link_name,.True.,.True.)

     !** Worker Loop **********************************************************
     Do

        Call mpi_recv(Active, 1_mik, mpi_integer4, 0_mik, rank_mpi, &
                      MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Active didn't succseed", INT(ierr, KIND=ik))

        If (Active == -1) Exit

        CALL mpi_recv(nn, 1_mik, mpi_integer8, 0_mik, rank_mpi, &
                      MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Domain didn't succeed", INT(ierr, KIND=ik))

        Domain = Domains(nn)
        
        if (out_amount /= "PRODUCTION") then
           !** >> Recieve Job_Dir << ******************************************
           CALL mpi_recv(job_dir, 4_mik*int(mcl,mik), mpi_character, 0_mik, &
                rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
           CALL print_err_stop(std_out, "MPI_RECV on Domain path didn't succeed", INT(ierr, KIND=ik))

        Else
           job_dir = outpath

        End if
        
        if (job_dir(len(job_dir):len(job_dir)) /= "/") then
           job_dir = trim(job_dir)//"/"
        End if

        ! DEBUG INFORMATION
        If (out_amount == "DEBUG") THEN
           Write(un_lf, fmt_dbg_sep)
           Write(un_lf, fmt_MSG_xAI0)"Root pointer before exec_single_domain on proc ",rank_mpi
           Call log_tree(root, un_lf, .True.)
           Write(un_lf, fmt_dbg_sep)
        END If
        ! DEBUG INFORMATION

        !======================================================================
        Call exec_single_domain(root, nn, Domain, job_dir, Active, fh_mpi, &
             worker_rank_mpi, worker_size_mpi, worker_comm)
        !======================================================================

        ! DEBUG INFORMATION
        If (out_amount == "DEBUG") THEN
           Write(un_lf, fmt_dbg_sep)
           Write(un_lf, fmt_MSG_xAI0)"Root pointer after exec_single_domain on proc ",rank_mpi
           Call log_tree(root,un_lf,.True.)
           Write(un_lf, fmt_dbg_sep)
        END If
        ! DEBUG INFORMATION

        !======================================================================
        !== Organize Results
        !======================================================================

        !** Look for the Domain branch ****************************************
        domain_desc=''
        Write(domain_desc, FMT_TXT_AxI0)'Domain ',Domain
        
!!$        Call search_branch(trim(domain_desc), root, db, success)

!!$        !** If we want to keep all results ************************************
!!$        If (out_amount == "DEBUG") then
!!$
!!$           Write(un_lf,fmt_dbg_sep)
!!$
!!$           !** Write input files **********************************************
!!$           
!!$           !** Search input files branch ************************
!!$           Call search_branch("Input files", db, meshb, success)
!!$
!!$           If (success) then
!!$
!!$              Do ii = 1, meshb%no_branches
!!$
!!$                 tmp_un = give_new_unit()
!!$                 
!!$                 !write(tmp_fn,fmt_filename) trim(job_dir)//trim(project_name)//'_',&
!!$                 !      nn,'_',ii,'.dat'
!!$             
!!$                 Open(unit=tmp_un, file=trim(tmp_fn), action="write", &
!!$                      status="replace")
!!$                 
!!$                 Write(tmp_un,*) meshb%branches(ii)%leaves(3)%p_char
!!$
!!$                 Close(tmp_un)
!!$
!!$              End Do
!!$              
!!$           End If
!!$
!!$           Write(un_lf,fmt_dbg_sep)           
!!$
!!$        End If
!!$
!!$        !** Dump the Local domain decomposition ****************************
!!$        mesh_desc=''
!!$        Write(mesh_desc,'(A,I0)')'Local domain Decomposition of domain no ',Domain
!!$           
!!$        removed_data = 0
!!$        Call delete_branch_from_branch(trim(mesh_desc), db, removed_data)
!!$
!!$        !** Dump the mesh and results branch *******************************
!!$        mesh_desc=''
!!$        Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_',Domain
!!$        
!!$        Call search_branch(trim(mesh_desc), root, meshb, success)
!!$           
!!$        If (success) then
!!$
!!$           !** Move No of nodes, elements and constrained dofs *************
!!$           Call add_leaf_to_branch(db,3)
!!$           Call Assign_leaves(db%leaves, meshb%leaves)
!!$           
!!$           removed_data = 0
!!$           Call delete_branch_from_branch(trim(mesh_desc), db, removed_data)
!!$           
!!$        End If
!!$           
!!$        !** Dump the input files branch ************************************
!!$        Call search_branch("Input files", db, meshb, success)
!!$           
!!$        If (success) then
!!$
!!$           removed_data = 0
!!$           Call delete_branch_from_branch("Input files", db, removed_data)
!!$           
!!$        End If

        ! DEBUG INFORMATION
!!$        If (out_amount == "DEBUG") THEN
!!$           
!!$           Write(un_mon,'(2(A,I10))')"MPI rank : ",rank_mpi, &
!!$                " ; Domain number : ",Domain
!!$           
!!$           Call get_data_size(root,data_size)
!!$           
!!$           Write(un_mon,'(A,3(I10,A))')"Data size root : ", &
!!$                data_size(1), " Elements ; ",&
!!$                data_size(2) / 1024_ik," kB ; ",&
!!$                data_size(2), " Byte"
!!$           
!!$           Call get_data_size(db,data_size)
!!$           
!!$           Write(un_mon,'(A,3(I10,A))')"Data size db   : ", &
!!$                data_size(1), " Elements ; ",&
!!$                data_size(2) / 1024_ik," kB ; ",&
!!$                data_size(2), " Byte"
!!$                   
!!$           Write(un_lf,fmt_dbg_sep)
!!$           Write(un_lf,fmt_MSG_AI0)"Root pointer after compress on proc",rank_mpi
!!$           Call log_tree(root,un_lf,.True.)
!!$           Write(un_lf,fmt_dbg_sep)
!!$           
!!$        END If
        ! DEBUG INFORMATION
           
        Call MPI_ISEND(Active, 1_mik, MPI_INTEGER4, 0_mik, rank_mpi, MPI_COMM_WORLD, REQUEST, IERR)
        CALL print_err_stop(std_out, "MPI_ISEND on Active didn't succeed", INT(ierr, KIND=ik))

        Call MPI_WAIT(REQUEST, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_WAIT on request for ISEND Active didn't succeed", INT(ierr, KIND=ik))

     End Do

     CALL PetscFinalize(petsc_ierr)
     
  End If

1000 Continue
  !============================================================================
  
  ! DEBUG INFORMATION
  If (out_amount == "DEBUG") THEN
     Write(un_lf, fmt_dbg_sep)
     Write(un_lf, fmt_MSG_xAI0)"Final Root pointer proc",rank_mpi
     Call log_tree(root,un_lf,.True.)
     Write(un_lf, fmt_dbg_sep)
  END If
  ! DEBUG INFORMATION
  
  Call link_end(link_name,.True.)

1001 Continue

IF(rank_mpi == 0) THEN
   CALL meta_signing(binary)
   CALL meta_close(size_mpi)

   ! CALL meta_stop_ascii(fh=fhl  , suf=log_suf)
   CALL meta_stop_ascii(fh=fhmon, suf=mon_suf)

   IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (rank_mpi == 0)

Call MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, KIND=ik))

End Program main_struct_process
