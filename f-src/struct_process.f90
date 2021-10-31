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
!>  on : 29.10.2021
!------------------------------------------------------------------------------
Module sp_aux_routines

  Use ISO_C_BINDING
  Use Operating_System
  USE global_std
  use puredat_com
  use chain_routines
  use strings
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

    LOGICAL, PARAMETER                :: DEBUG=.TRUE.
    
    Integer(kind=mpi_ik), intent(out) :: Active

    Type(tBranch), Intent(inOut)      :: root
    Integer(kind=ik), Intent(in)      :: nn, lin_nn
    Character(LEN=*), Intent(in)      :: job_dir
    Integer(kind=mpi_ik), Intent(In), Dimension(no_streams) :: fh_mpi
    Integer(kind=mpi_ik), Intent(In)  :: rank_mpi, size_mpi, comm_mpi

    !----------------------------------------------------------------
    Integer(kind=mpi_ik)              :: ierr
    Integer(kind=mpi_ik), Dimension(MPI_STATUS_SIZE)   :: status_mpi
    Type(tBranch), pointer            :: bb, db, pb, mb, params,resb

    Character(Len=mcl)                :: desc, mesh_desc, filename
    Character(Len=mcl)                :: elt_micro
    
    Character, Dimension(4*mcl)       :: c_char_array

    Integer  (kind=c_int )            :: stat_c_int
    Integer                           :: stat
    Integer                           :: umon

    Integer, Dimension(8)             :: realt
    Character(len=mcl)                :: env_var

    Character(Len=mcl)                :: out_file

    Integer(kind=ik)                  :: m_size

    Character(len=2)                  :: ii_char
    Character(len=9)                  :: nn_char

    logical                           :: success

    Character(len=mcl)                :: timer_name, domain_desc, part_desc

    Integer(kind=mpi_ik)              :: petsc_ierr
    Type(tMat)                        :: AA, AA_org
    Type(tVec)                        :: XX
    Type(tVec), Dimension(24)         :: FF
    TYPE(tPETScViewer)                :: V
    Type(tKSP)                        :: KSP
    Integer(Kind=ik)                  :: Istart,Iend, parts, IVstart, IVend
    Integer(Kind=ik), Dimension(:), Allocatable :: nodes_in_mesh, elems_in_mesh

    Integer(kind=pd_ik), Dimension(:), Allocatable :: serial_pb
    Integer(kind=pd_ik)                            :: serial_pb_size

    Integer(Kind=pd_ik)  , Dimension(no_streams)  :: no_data
    Integer(kind=ik)                              :: nn_elems, ii, jj, kk, id
    Integer(kind=ik), Dimension(:), Allocatable   :: gnid_cref,loc_nn_cref
    Integer(kind=ik), Dimension(:,:), Allocatable :: vtk_topo, res_sizes
    Real   (kind=rk), Dimension(:), Pointer       :: displ, force
    Real   (kind=rk), Dimension(:), Allocatable   :: glob_displ, glob_force
    Real   (kind=rk), Dimension(:), Allocatable   :: zeros_R8
    Integer(kind=ik)                              :: cr, idum, tmp_i8

    Integer(kind=mpi_ik)                           :: rank_mpi_tmp, size_mpi_tmp
    Character, Dimension(:), Allocatable           :: char_arr
    
    !--------------------------------------------------------------------------

    Integer(Kind=ik), Dimension(60)    :: idxm_20, idxn_20
    Real(kind=rk),    Dimension(60,60) :: K_loc_20

    Integer(Kind=ik), Dimension(24)    :: idxm_08, idxn_08
    Real(kind=rk),    Dimension(24,24) :: K_loc_08
           
    !--------------------------------------------------------------------------

    !** Init Activity status ****
    Active = 0_mpi_ik

    write(nn_char,'(I0)')nn

    !** Get basic infos ------------------------------------------
    Call Search_branch("Input parameters", root, params, success, DEBUG)
    call pd_get(params,"No of mesh parts per subdomain",parts)
    
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

             Call execute_command_line("mkdir -p "//trim(job_dir), CMDSTAT=stat)

             If ( stat /= 0 ) Then
                Write(un_mon,*)"Could not execute syscall"
                Write(un_mon,*)"mkpir -p "//trim(job_dir)
                Write(un_mon,*)"Program halted"
                Stop
             End If

             Call Stat_Dir(c_char_array, stat_c_int)

             If ( stat_c_int /= 0 ) Then
                Write(un_mon,*)"Could not create directory"
                Write(un_mon,*)trim(job_dir)
                Write(un_mon,*)"Program halted"
                Stop
             End If

          Else

             !***********************************************************************
             !usta = give_new_unit()
             !.... Check sta file

          End If

          Write(un_lf,fmt_sep)

       End If

       !**************************************************************************
       !** Write log and monitor file
       !**************************************************************************
       Write(un_lf,fmt_msg_AI0) "Domain No. : ",nn
       Write(un_lf,fmt_msg_A  ) "Job_dir    : "//Trim(job_dir)

       Call date_and_Time(values=realt)
       Write(un_lf,"('MM ',A,I0,2('.',I0),' - ',I0,2(':',I0),',',I0)") &
            "Start time : ",realt(3),realt(2),realt(1),realt(5:8)

       Call get_environment_Variable("HOSTNAME", env_var)
       Write(un_lf,fmt_MSG_A) "Host       : "//Trim(env_var)

       umon = give_new_unit()

       Open(unit=umon, file=Trim(job_dir)//Trim(out%bsnm)//".mon",action="write", &
            status="replace")

       Write(umon,fmt_eq_sep)
       Write(umon,fmt_msg_AI0) "Domain No. : ",nn

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

       Call date_and_Time(values=realt)
       Write(un_lf,"('MM ',A,I0,2('.',I0),' - ',I0,2(':',I0),',',I0)") &
            "End time   : ",realt(3),realt(2),realt(1),realt(5:8)    

       close(umon)

       !** Look for the Domain branch ****************************************
       domain_desc=''
       Write(domain_desc,'(A,I0)')'Domain ',nn
       
       Call search_branch(trim(domain_desc), root, db, success, DEBUG)

       !** Get the no of nodes per part **************************************
       mesh_desc = ''
       Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(out%bsnm)//'_',nn
       
       Call search_branch(trim(mesh_desc), db, mb, success, DEBUG)
       call pd_get(mb, 'No of nodes in mesh',  nodes_in_mesh)

       !** Set the global matrix size ****************************************
       m_size = nodes_in_mesh(1) * 3

       Do ii = 1, parts-1

          ! Look for the Part branch *************************
          domain_desc=''
          Write(part_desc,'(A,I0)')'Part_',ii
          
          !** serialize branch with mesh part *********************************
          Call search_branch(trim(part_desc), db, pb, success)
          If (.NOT. success) Then
             Write(*,*)"Something bad and unexpected happend in exec_single_domain !!!"
             Write(*,*)"Looking for branch of part ",ii," returned ",success
             Write(*,*)"MPI proc ",rank_mpi," halted !!!"
             stop
          End If
          
          Call serialize_branch(pb,serial_pb,serial_pb_size,.TRUE.)

          Call mpi_send(serial_pb_size, 1_mpi_ik, MPI_INTEGER8, Int(ii,mpi_ik), Int(ii,mpi_ik), &
               COMM_MPI, ierr)
          Call mpi_send(serial_pb, INT(serial_pb_size,mpi_ik), MPI_INTEGER8, &
               Int(ii,mpi_ik), Int(ii,mpi_ik), COMM_MPI, ierr)

          Deallocate(serial_pb)
          
       End Do

       part_desc=''
       Write(part_desc,'(A,I0)')'Part_',parts
          
       Call search_branch(trim(part_desc), db, pb, success, DEBUG)

       !** Broadcast matrix size. TODO could also be included into part branches.
       Call mpi_bcast(m_size, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, COMM_MPI, ierr)
       
    !****************************************************************************
    !** Ranks > 0 -- Workers ****************************************************
    !****************************************************************************
    Else

       Call mpi_recv(serial_pb_size, 1_mpi_ik, mpi_integer8, 0_mpi_ik, &
            rank_mpi, COMM_MPI, status_mpi, ierr)

       if (allocated(serial_pb)) deallocate(serial_pb)
       
       Allocate(serial_pb(serial_pb_size))

       Call mpi_recv(serial_pb, INT(serial_pb_size,mpi_ik), mpi_integer8, &
            0_mpi_ik, rank_mpi, COMM_MPI, status_mpi, ierr)

       !** Deserialize part branch ****************************
       Call Start_Timer("Deserialize part branch branch")

       Allocate(pb)
       
       Call deserialize_branch(pb, serial_pb, .TRUE.)

       Call End_Timer("Deserialize part branch branch")

       Call mpi_bcast(m_size, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, COMM_MPI, ierr)
              
    End If

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       Write(un_lf,fmt_dbg_sep)
       Write(un_lf,'(A)')"part branch right after deserialization"
       Call log_tree(pb,un_lf,.TRUE.)
       Write(un_lf,fmt_dbg_sep)
       flush(un_lf)
    END If
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
    
    !**************************************************************************
    !** Setup the linear System with a constant system matrix A. Once that  ***
    !** is done setup the multiple right hand sides and solve the linear    ***
    !** system multiple times. Save the solutions to calculate effective    ***
    !** stiffness matirces.                                                 ***
    !**************************************************************************

    !** Create Stiffness matrix **************************************
    call MatCreate(COMM_MPI, AA    , petsc_ierr)
    call MatCreate(COMM_MPI, AA_org, petsc_ierr)
  
    call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)
    call MatSetFromOptions(AA,petsc_ierr)
    call MatSetUp(AA,petsc_ierr)

    call MatSetSizes(AA_org,PETSC_DECIDE,PETSC_DECIDE,m_size,m_size,petsc_ierr)
    call MatSetFromOptions(AA_org,petsc_ierr)
    call MatSetUp(AA_org,petsc_ierr)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       call MatGetOwnershipRange(AA, Istart , Iend,  petsc_ierr)
       Write(un_lf,'(A,I4,A,A6,2(A,I9))')&
            "MPI rank : ",rank_mpi,"| matrix ownership for ","A    ", ":",Istart," -- " , Iend
       call MatGetOwnershipRange(AA_org, Istart , Iend,  petsc_ierr)
              Write(un_lf,'(A,I4,A,A6,2(A,I9))')&
            "MPI rank : ",rank_mpi,"| matrix ownership for ","A_org", ":",Istart," -- " , Iend
    End If
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       

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

    if (elt_micro == "HEX08") then

       K_loc_08 = Hexe08()

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
    
    Call MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    Call MatAssemblyBegin(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    ! Computations can be done while messages are in transition
    Call MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY ,petsc_ierr)
    Call MatAssemblyEnd(AA_org, MAT_FINAL_ASSEMBLY ,petsc_ierr)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       Call PetscViewerCreate(COMM_MPI, V, petsc_ierr)
       Call PetscViewerASCIIOpen(COMM_MPI,"AA.output.1",V, petsc_ierr);
       Call PetscViewerSetFormat(V, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
       Call MatView(AA, V, petsc_ierr)
       Call PetscViewerDestroy(V, petsc_ierr)
    End If
    
    !***************************************************************************
    !** At this point the system matrix is assembled. To make it ready to be ***
    !** used, the rows and columns of the dofs with prescribed displacements ***
    !** have to be eliminated. To do that with MatZeroRowsColumns we need    ***
    !** right hand side vectors and a solution vector.                       ***
    !***************************************************************************

    Do ii = 1, 24
    
       !** Create load vectors ***************************************
       call VecCreate(COMM_MPI, FF(ii), petsc_ierr)
       call VecSetSizes(FF(ii), PETSC_DECIDE, m_size, petsc_ierr)
       call VecSetFromOptions(FF(ii), petsc_ierr)
       Call VecSet(FF(ii), 0._rk,petsc_ierr)
       
       !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       If (out_amount == "DEBUG") THEN 
          call VecGetOwnershipRange(FF(ii), IVstart, IVend, petsc_ierr)
          Write(un_lf,'(A,I4,A,A6,2(A,I9))')&
               "MPI rank : ",rank_mpi,"| vector ownership for ","F", ":",IVstart," -- " , IVend
       End If
       !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
       
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

    !** Set prescribed displacements of LC1 to solution vector *******
    Call VecSetValues(XX, bb%leaves(2)%dat_no, &
         gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)

    Call VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition
    Call VecAssemblyEnd(XX, petsc_ierr)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       call VecGetOwnershipRange(XX, IVstart, IVend, petsc_ierr)
       Write(un_lf,'(A,I4,A,A6,2(A,I9))')&
            "MPI rank : ",rank_mpi,"| vector ownership for ","X", ":",IVstart," -- " , IVend
       Call PetscViewerCreate(COMM_MPI, V, petsc_ierr)
       Call PetscViewerASCIIOpen(COMM_MPI,"FX.output.1",V, petsc_ierr);
       Call PetscViewerSetFormat(V, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
       Call VecView(XX, V, petsc_ierr)
       Call PetscViewerDestroy(V, petsc_ierr)
    End If
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !***************************************************************************
    !** At this point the right hand side vectors, filled with zeros and     ***
    !** the solution vector filled with the dirichlet boundary values of     ***
    !** load case 1 (LC1) are ready to be used.                              ***
    !** Now the righ hand side vectors have to be modified with the          ***
    !** prescribed displacements. This is currently done by multiplying A by ***
    !** X filled with the negative displacement values.                      ***
    !** (See above VecSetValues(XX, ... , -bb%leaves(2)%p_real8, ... )       ***
    !***************************************************************************
    
    !** Compute dirichlet boundary corrections of first right  *******
    !** hand side vector.                                      *******
    Call MatMult(AA,XX,FF(1), petsc_ierr);

    !** Set zero values for dofs with prescribed displacements *******
    Call VecSetValues(FF(1), bb%leaves(2)%dat_no, &
         gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
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

       !** Set prescribed displacements of LCii to solution vector *****
       Call VecSetValues(XX, bb%leaves(2)%dat_no, &
            gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
    
       Call VecAssemblyBegin(XX, petsc_ierr)
       ! Computations can be done while messages are in transition
       Call VecAssemblyEnd(XX, petsc_ierr)

       !** Compute dirichlet boundary corrections of ii th right *******
       !** hand side vector.                                     *******
       Call MatMult(AA,XX,FF(ii), petsc_ierr);

       !** Set zero values for dofs with prescribed displacements ******
       Call VecSetValues(FF(ii), bb%leaves(2)%dat_no, &
            gnid_cref, zeros_R8, INSERT_VALUES, petsc_ierr)
       Call VecAssemblyBegin(FF(ii), petsc_ierr)
       
    End Do

    !** Get Bounds branch of LC 24 **********************************
    !** bb%leaves(1)%p_int8  : Boundary displacement node global ids
    !** bb%leaves(2)%p_real8 : Boundary displacement values
    write(desc,'(A,I0)')"Boundaries_"//trim(nn_char)//"_",24
    Call search_branch(trim(desc), pb, bb, success, DEBUG)
    
    !** Set prescribed displacements of LC24 to solution vector *****
    Call VecSetValues(XX, bb%leaves(2)%dat_no, &
         gnid_cref, -bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)

    Call VecAssemblyBegin(XX, petsc_ierr)
    ! Computations can be done while messages are in transition

    !** Finalize the open assembleys *************
    Do ii = 2, 23
       Call VecAssemblyEnd(FF(ii), petsc_ierr)
    End Do

    Call VecAssemblyEnd(XX, petsc_ierr)

    !** Since we are filling XX with the prescribed displacements  ***
    !** times -1. , we have to rescale XX before using it in       ***
    !** MatZeroRowsColumns.                                        ***
    Call VecScale(XX, -1._rk, petsc_ierr)
    
    !** Apply Dirichlet Boundaries to A and the 24th right hand    ***
    !** side vector                                                ***
    call MatZeroRowsColumns(AA, bb%leaves(2)%dat_no, gnid_cref, &
         0.0_8, XX, FF(24), petsc_ierr)
    
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       Call PetscViewerCreate(COMM_MPI, V, petsc_ierr)
       Call PetscViewerASCIIOpen(COMM_MPI,"FF.output.1",V, petsc_ierr);
       Call PetscViewerSetFormat(V, PETSC_VIEWER_ASCII_DENSE, petsc_ierr)
       Call VecView(FF( 1), V, petsc_ierr)
       Call PetscViewerDestroy(V, petsc_ierr)
    End If
   
    !***************************************************************************
    !** At this point the right hand side vectors are modified with the       **
    !** prescribed displacements and the rows and columns of A representing   **
    !** dofs with prescribed displacements are filled with zeros.             **
    !** Now a linear solver context can be set up and the linear systems with **
    !** constant operator A and variable right hand side can be solved        **
    !***************************************************************************
    
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
    
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN 
       
       if ( rank_mpi == 0 ) then
          filename=''
          write(filename,'(A,I0,A)')trim(job_dir)//trim(out%bsnm)//"_",nn,"_usg.vtk"
          Call write_data_head(filename, m_size/3)
       End if

    End If
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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
       Write(part_desc,'(A,I0)')'Part_',parts
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
       
       !** Extend XX with Boundary displacements *********************
       Call VecSetValues(XX, bb%leaves(2)%dat_no, &
            gnid_cref, bb%leaves(2)%p_real8, INSERT_VALUES, petsc_ierr)
       
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
          Call mpi_send([IVstart,IVend], 2_mpi_ik, &
               MPI_Integer8,  0_mpi_ik, rank_mpi, COMM_MPI, ierr)
          
          Call mpi_send(displ, Int(IVend-IVstart,mpi_ik), &
               MPI_Real8,  0_mpi_ik, Int(rank_mpi+size_mpi,mpi_ik), &
               COMM_MPI, ierr)

          Call mpi_send(force, Int(IVend-IVstart,mpi_ik), &
               MPI_Real8,  0_mpi_ik, Int(rank_mpi+2*size_mpi,mpi_ik), &
               COMM_MPI, ierr)
       Else
          
          ! Copy rank 0 local result *********************************
          glob_displ(IVstart:IVend-1) = displ
          glob_force(IVstart:IVend-1) = force
          
          ! Recv bounds of local results *****************************
          Do ii = 1, size_mpi-1
             Call mpi_recv(res_sizes(:,ii), 2_mpi_ik, &
                  MPI_Integer8,  Int(ii, mpi_ik), Int(ii, mpi_ik), &
                  COMM_MPI, status_mpi, ierr)
          End Do
          
          ! Recv local parts of global result ************************
          Do ii = 1, size_mpi-1
             Call mpi_recv(glob_displ(res_sizes(1,ii):res_sizes(2,ii)-1), &
                  Int(res_sizes(2,ii)-res_sizes(1,ii),mpi_ik), &
                  MPI_Integer8,  Int(ii, mpi_ik), Int(ii+size_mpi, mpi_ik), &
                  COMM_MPI, status_mpi, ierr)

             Call mpi_recv(glob_force(res_sizes(1,ii):res_sizes(2,ii)-1), &
                  Int(res_sizes(2,ii)-res_sizes(1,ii),mpi_ik), &
                  MPI_Integer8,  Int(ii, mpi_ik), Int(ii+2*size_mpi, mpi_ik) , &
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

          !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
          !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 
       End If
    End Do
   
    !*****************************************************************
    !** All 24 linear system solutions are produced. Effective    ****
    !** stiffnesses can be calculated                             ****
    !*****************************************************************

    if ( rank_mpi == 0 ) then

       Deallocate(glob_displ, res_sizes, glob_force)
       
       Call start_timer(trim(timer_name), .FALSE.)
       call calc_effective_material_parameters(root, lin_nn, nn, job_dir, fh_mpi)
       Call end_timer(trim(timer_name))
       
    End if
    
!!!!!!!!!!!!!<< Development <<<<<<<<<<<<<<<<<<<<
    
    call MatDestroy(AA    , petsc_ierr)
    Call VecDestroy(XX    , petsc_ierr)
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

!==============================================================================
!> Struct Process main programm
!--
!------------------------------------------------------------------------------
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on : 29.10.2021
!>
!> \todo At >> Recieve Job_Dir << calculate the Job_Dir from the Domain number
!>       and the Domain Decomposition parameters.
!> \todo During geometry generation do not allocate the Topology field for the
!>       max possible number of elements. Instead check with count for actual
!>       number of elements
!------------------------------------------------------------------------------
Program main_struct_process

  USE global_std
  USE global_pd
  USE puredat 
  USE auxiliaries
  USE meta
  USE chain_routines
  USE Operating_System
  USE MPI
  USE ISO_C_BINDING
  USE decomp 
  USE sp_aux_routines
  USE PETSC
  USE petsc_opt
  
  Implicit None

  !-- MPI Variables -----------------------------------------------------------
  Integer(kind=mpi_ik)                               :: ierr
  Integer(kind=mpi_ik)                               :: rank_mpi, size_mpi
  Integer(kind=mpi_ik)                               :: worker_rank_mpi, worker_size_mpi
  
  Integer(kind=mpi_ik)                               :: Active

  Integer(kind=mpi_ik)                               :: request

  Integer(kind=mpi_ik)                               :: finished

  Integer(kind=mpi_ik), Dimension(no_streams)        :: fh_mpi

  Integer(kind=mpi_ik), Dimension(MPI_STATUS_SIZE)   :: status_mpi
  Integer(kind=mpi_ik), Allocatable, Dimension(:,:)  :: statuses_mpi

  Integer(kind=mpi_ik), Allocatable, Dimension(:)    :: Activity
  Integer(kind=mpi_ik), Allocatable, Dimension(:)    :: req_list
  Integer(kind=mpi_ik), Allocatable, Dimension(:)    :: members

  Integer(kind=mpi_ik)                               :: worker_comm

  !----------------------------------------------------------------------------
  Character(Len=mcl)                                :: link_name = 'struct process'
  Character, Dimension(4*mcl)                       :: c_char_array
  Integer  (kind=c_int )                            :: stat_c_int
  Integer                                           :: stat
  Type(tBranch)                                     :: root
  Type(tBranch)                                     :: phi_tree
  Type(tBranch), pointer                            :: ddc
  Type(tBranch), Pointer                            :: params, epp, meshb, db, res
  
  Integer(kind=ik)    , Allocatable, Dimension(:)   :: nn_D
 
  Character(len=mcl)                                :: fmps_epp
  Character(LEN=4*mcl)                              :: job_dir, tmp_fn
  CHARACTER                                         :: restart='N'
  Character(LEN=4*mcl), Dimension(:), Allocatable   :: domain_path

  ! Meta file variable
  CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE     :: m_rry      
  CHARACTER(LEN=mcl)                                :: infile='', int_id, file_hrdcd
  CHARACTER(LEN=mcl)                                :: cmd_arg='notempty', cmd_arg_history=''

 
  Character(LEN=mcl)             :: muCT_pd_path
  Character(LEN=mcl)             :: muCT_pd_name
  Character(Len=mcl)             :: mesh_desc, domain_desc

  Real(kind=rk)   , Dimension(3) :: pdsize

  Integer(kind=ik), Dimension(3) :: xa_d
  Integer(kind=ik), Dimension(3) :: xe_d
  
  Integer(kind=ik)               :: nn, ii, jj, kk, dc
  Integer(kind=ik)               :: Domain_number, path_count
  Integer(kind=ik)               :: alloc_stat, iun, aun, tmp_un
  Integer(kind=ik), Dimension(2) :: data_size

  Integer(kind=ik), Allocatable, Dimension(:)   :: Domains, Domain_stats, act_domains
  Integer(kind=ik)                              :: Domain
  Integer(kind=ik)                              :: llimit, parts, elo_macro
  Integer(kind=ik)                              :: no_solver, pscratch
  Character(LEN=8)                              :: elt_micro
  Character(Len=8)                              :: output
  Real(kind=rk)                                 :: strain, e_modul, nu
  
  Integer(kind=pd_ik), Dimension(:), Allocatable :: serial_root, epp_data
  Integer(kind=pd_ik)                            :: serial_root_size
  Integer(kind=pd_ik)                            :: npos,mxpt,mpnu
  Integer(kind=pd_ik)                            :: ipos, i, ip
  Integer(kind=mpi_ik)                           :: petsc_ierr

  Logical                                        :: success, fexist
  Integer(kind=pd_ik), Dimension(no_streams)     :: removed_data, dsize
  Character, Dimension(:), Allocatable           :: char_arr

  !type(Mat)                         :: A
  Integer(Kind=ik)                  :: A
  Integer(Kind=ik)                  :: Istart,Iend,ione
  
  !----------------------------------------------------------------------------
 
  Call mpi_init(ierr)
  CALL handle_err(std_out, "MPI_INIT didn't succeed", INT(ierr, KIND=ik))

  Call MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)
  CALL handle_err(std_out, "MPI_COMM_RANK couldn't be retrieved", INT(ierr, KIND=ik))
 
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
  CALL handle_err(std_out, "MPI_COMM_SIZE couldn't be retrieved", INT(ierr, KIND=ik))
 
  If (size_mpi < 2) CALL handle_err(std_out, "We need at least 2 MPI processes to execute this program.", 1)
  
  !------------------------------------------------------------------------------
  ! Rank 0 -- Init (Master) Process and broadcast init parameters 
  !------------------------------------------------------------------------------
  If (rank_mpi==0) Then

      Call Start_Timer("Init Process")
    
      !------------------------------------------------------------------------------
      ! Parse the command arguments
      !------------------------------------------------------------------------------

      IF (command_argument_count() == 0) CALL usage(1)

      DO ii=0, 10 ! Read up to 10 command arguments.

         CALL GET_COMMAND_ARGUMENT(ii, cmd_arg)

         cmd_arg_history = TRIM(cmd_arg_history)//' '//TRIM(cmd_arg)

         IF (cmd_arg(1:1) .EQ. '-') THEN
            DO jj=2, LEN_TRIM(cmd_arg)
                  SELECT CASE( cmd_arg(jj:jj) )
                     CASE('y')
                        restart = 'Y'
                     ! It is highly recommended not to use this option on the cluster  :-)
                     CASE('h')       ! CASE ('-h', '--help')
                        CALL usage(1)
                  END SELECT
            END DO
         END IF

         ! Simple check whether the cmd arg can be a meta file since no other flag is longer than 3 characters.
         ! May crash when the flags are extended improperly :-)
         ! 
         ! the check against the meta suffix happens at the CALL meta_append routine.
         IF (LEN_TRIM(cmd_arg) .GT. 5_ik) infile=cmd_arg

         IF(cmd_arg == '') EXIT           
      END DO

      IF(TRIM(infile) == '') CALL handle_err(std_out, 'No input file given via command argument.', 1)

      in%full = TRIM(infile)

      !------------------------------------------------------------------------------
      ! Check and open the input file; Modify the Meta-Filename / Basename
      !------------------------------------------------------------------------------
      CALL meta_append(restart, m_rry)
      
      !------------------------------------------------------------------------------
      ! Structure output directories and namings
      !------------------------------------------------------------------------------
      out%path = in%path(1:LEN_TRIM(in%path)-1)//"_reg/"
      out%bsnm = TRIM(in%bsnm)

      !------------------------------------------------------------------------------
      ! Spawn a log file and a results file
      !------------------------------------------------------------------------------
      CALL meta_add_ascii(fh=fhl , suf=log_suf, st='start', restart=restart)
      CALL meta_add_ascii(fh=fhmo, suf=mon_suf, st='start', restart=restart)
      ! CALL meta_add_ascii(fh=fhr, suf=res_suf, st='start', restart=restart)

      ! CALL meta_io (fhmo, 'BASE_PATH'   , '(-)'  , m_rry, chars = out%path   , wl=.TRUE.)
      ! CALL meta_io (fhmo, 'PROJECT_NAME', '(-)'  , m_rry, chars = out%bsnm   , wl=.TRUE.)



      ! outpath      = TRIM(in%path)//"_reg/"
      ! project_name = TRIM(in%bsnm)

      Allocate(params)
      Call raise_tree("Input parameters",params)

      !------------------------------------------------------------------------------
      ! Read input parameters
      !------------------------------------------------------------------------------
      ! CALL meta_io (fhmo, 'MCT_PD_PRO_PATH'  , '(-)'  , m_rry,    chars = muCT_pd_path, wl=.TRUE.)
      ! CALL meta_io (fhmo, 'MCT_PD_PRO_NAME'  , '(-)'  , m_rry,    chars = muCT_pd_name, wl=.TRUE.)
      CALL meta_io (fhmo, 'MICRO_ELMNT_TYPE' , '(-)'  , m_rry,    chars = elt_micro   , wl=.TRUE.)
      CALL meta_io (fhmo, 'DBG_LVL'          , '(-)'  , m_rry,    chars = out_amount  , wl=.TRUE.)
      CALL meta_io (fhmo, 'OUT_FMT'          , '(-)'  , m_rry,    chars = output      , wl=.TRUE.)
      CALL meta_io (fhmo, 'SIZE_DOMAIN'      , '(mm)' , m_rry, real_1D3 = pdsize      , wl=.TRUE.)
      CALL meta_io (fhmo, 'LO_BNDS_DMN_RANGE', '(-)'  , m_rry,  int_1D3 = xa_d        , wl=.TRUE.)
      CALL meta_io (fhmo, 'UP_BNDS_DMN_RANGE', '(-)'  , m_rry,  int_1D3 = xe_d        , wl=.TRUE.)
      CALL meta_io (fhmo, 'BINARIZE_LO'      , '(-)'  , m_rry,  int_0D  = llimit      , wl=.TRUE.)
      CALL meta_io (fhmo, 'MESH_PER_SUB_DMN' , '(-)'  , m_rry,  int_0D  = parts       , wl=.TRUE.)
      CALL meta_io (fhmo, 'RVE_STRAIN'       , '(mm)' , m_rry, real_0D  = strain      , wl=.TRUE.)
      CALL meta_io (fhmo, 'YOUNG_MODULUS'    , '(MPa)', m_rry, real_0D  = e_modul     , wl=.TRUE.)
      CALL meta_io (fhmo, 'POISSON_RATIO'    , '(-)'  , m_rry, real_0D  = nu          , wl=.TRUE.)
      CALL meta_io (fhmo, 'MACRO_ELMNT_ORDER', '(-)'  , m_rry,  int_0D  = elo_macro   , wl=.TRUE.)

      ! Error handling
      IF ( (xa_d(1) > xe_d(1)) .OR. (xa_d(1) > xe_d(1)) .or. (xa_d(1) > xe_d(1)) ) THEN
         CALL handle_err(std_out, 'Input parameter error: Start value of domain range larger than end value.', 1)
      END IF

      IF (MOD(size_mpi-1, parts) .NE. 0) CALL handle_err(std_out, 'Please provide more domains than processors.', 1)

CALL skip(std_out)
WRITE(*,*) "muCT_pd_path: ", TRIM(in%path)
WRITE(*,*) "muCT_pd_name: ", TRIM(in%bsnm)
CALL skip(std_out)

      CALL add_leaf_to_branch(params, "muCT puredat pro_path"                , mcl            , str_to_char(in%p_n_bsnm))
      CALL add_leaf_to_branch(params, "muCT puredat pro_name"                , mcl            , str_to_char(in%bsnm))
      CALL add_leaf_to_branch(params, 'Element type  on micro scale'         , len(elt_micro) , str_to_char(elt_micro))     
      CALL add_leaf_to_branch(params, 'Output amount'                        , len(out_amount), str_to_char(out_amount))
      CALL add_leaf_to_branch(params, 'Restart'                              , len(restart)   , str_to_char(restart))
      CALL add_leaf_to_branch(params, 'Output Format'                        , len(output)    , str_to_char(output))
      CALL add_leaf_to_branch(params, "Physical domain size"                 , 3_ik           , pdsize)
      CALL add_leaf_to_branch(params, "Lower bounds of selected domain range", 3_ik           , xa_d)
      CALL add_leaf_to_branch(params, "Upper bounds of selected domain range", 3_ik           , xe_d)     
      CALL add_leaf_to_branch(params, 'Lower limit of iso value'             , 1_ik           , [llimit])     
      CALL add_leaf_to_branch(params, 'No of mesh parts per subdomain'       , 1              , [parts])
      CALL add_leaf_to_branch(params, 'Average strain on RVE'                , 1              , [strain])   
      CALL add_leaf_to_branch(params, "Young_s modulus"                      , 1              , [e_modul])
      CALL add_leaf_to_branch(params, "Poisson_s ratio"                      , 1              , [nu])
      CALL add_leaf_to_branch(params, 'Element order on macro scale'         , 1              , [elo_macro])
     

     
     !** Prepare output directory ****************************************
     c_char_array(1:len(Trim(out%path)//Char(0))) = str_to_char(Trim(out%path)//Char(0))

     Call Stat_Dir(c_char_array, stat_c_int)

     If ( stat_c_int /= 0 ) Then

        Call execute_command_line("mkdir -p "//trim(out%path),CMDSTAT=stat)

        If ( stat /= 0 ) Then
           Write(un_mon,*)"Could not execute syscall"
           Write(un_mon,*)"mkpir -p "//trim(out%path)
           Write(un_mon,*)"Program halted"
           success = .FALSE.
        End If

        Call Stat_Dir(c_char_array, stat_c_int)

        If ( stat_c_int /= 0 ) Then
           Write(un_mon,*)"Could not create directory"
           Write(un_mon,*)trim(out%path)
           Write(un_mon,*)"Program halted"
           success = .FALSE.
        End If

        success = .TRUE.
        
     Else If ( (stat_c_int == 0) .AND. (Restart == "Y") ) Then

        Write(un_mon,FMT_MSG_A)"Reusing the output directory"
        Write(un_mon,FMT_MSG_A)trim(out%path)
        success = .TRUE.
        
     Else
        
        Write(std_out, FMT_ERR_A)"The output directory"
        Write(std_out, FMT_ERR_A)trim(out%path)
        Write(std_out, FMT_ERR_A)"apparently exists already with restart not equal to Y !!!"
        Write(std_out, FMT_ERR_A)"Please check your struct-process-parameters.sh file   !!!"
        write(std_out, FMT_STOP)
        success = .FALSE.
        
     End If

     If (.NOT. success) CALL handle_err(std_out, "Something went wrong during init of the output dir.", 1)

     Call link_start(link_name,.TRUE.,.FALSE., success)
     
     If (.NOT. success) CALL handle_err(std_out, "Something went wrong during link_start", 1)

     !** Set project name and path of puredat root ***
     pro_path = in%path
     pro_name = in%bsnm
        
     If ( restart == "N" ) then 

        !** Raise root branch *************************************************
        Call raise_tree(Trim(out%bsnm),root)

        Call include_branch_into_branch(s_b=params, t_b=root, blind=.TRUE.)
     
        !** Load puredat tree of micro-CT data and calculate the global
        !** parameters of the domain decomposition
        pro_path = in%path
        pro_name = in%bsnm

        phi_tree = read_tree()

        !** Set project name and path of global domain decomposition     
        pro_path = in%path
        pro_name = in%bsnm

        allocate(ddc)
        ddc = calc_general_ddc_params(pdsize, phi_tree)
        
        call include_branch_into_branch(s_b=ddc, t_b=root, blind=.TRUE.)
     
     Else if ( Restart == "Y" ) Then

        !** Check whether there is already a project header *******************

      !   inquire(file=Trim(pro_path)//Trim(pro_name)//'.head',exist=fexist)
        
        inquire(file=Trim(in%p_n_bsnm)//'.head',exist=fexist)
        If (.not.fexist) then
         Call End_Timer("Init Process")
         CALL handle_err(std_out, "Restart on job without PureDat root header ! This case is not supported.", 1)         
        End If
        
        !** Read root branch ******************************
        root = read_tree()

        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        If (out_amount == "DEBUG") THEN 
           Write(un_lf,fmt_dbg_sep)
           Write(un_lf,'(A)')"root right after restart read"
           Call log_tree(root,un_lf,.FALSE.)
           Write(un_lf,fmt_dbg_sep)
           flush(un_lf)
        END If
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !** !!! This calling sequence is only valid since "Averaged Material 
        !** !!! Properties" only contains r8 data added at the end of the 
        !** !!! r8-stream. More correct would be a routine that ensures data
        !** !!! integrity and compresses potentially missing stream data in
        !** !!! an efficient way.
        call delete_branch_from_branch("Averaged Material Properties", &
                                       root, dsize)
        Call get_stream_size(root, dsize)
        root%streams%dim_st = dsize
        root%streams%ii_st = dsize+1
        
        Call read_streams(root)
     
        Call connect_pointers(root%streams, root)

        Call search_branch("Global domain decomposition", root, ddc, success)
        
        If (.not. success) then
           Write(un_mon,FMT_ERR_A)"Found no branch named 'Global domain decomposition' !"
           Write(un_mon,FMT_ERR_A)"with restart option set to Y                        !"
           Write(un_mon,FMT_STOP)
           STOP
        End If

        Call search_branch("Input parameters", root, params, success)
        
        If (.not. success) then
           Write(un_mon,FMT_ERR_A)"Found no branch named 'Input parameters' !"
           Write(un_mon,FMT_ERR_A)"with restart option set to Y             !"
           Write(un_mon,FMT_STOP)
           STOP
        End If

        !** Reset Output amount and Restart in loaded param branch ************
        params%leaves(17)%p_char = str_to_char(out_amount)
        params%leaves(18)%p_char = "Y"

        !**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !** TODO Check that existing parameters match with given ones !!!!!
        !**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        If (out_amount == "DEBUG") THEN 
           Write(un_lf,fmt_dbg_sep)
           Write(un_lf,'(A)')"root right after restart and deletion of avg mat props branch"
           Call log_tree(root,un_lf,.FALSE.)
           Write(un_lf,fmt_dbg_sep)
           flush(un_lf)
        END If
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     
        
     Else
        
        Write(un_mon,FMT_ERR_A)"Only Y or N are support as input for the restart option !"
        Write(un_mon,FMT_STOP)
        STOP
        
     End If
        
     Call pd_get(ddc,"nn_D",nn_D)

     !** Allocate and init field for selected domain numbers ******************
     Domain_number = (xe_d(1)-xa_d(1)+1) * &
                     (xe_d(2)-xa_d(2)+1) * &
                     (xe_d(3)-xa_d(3)+1)
     
     Allocate(Domains(Domain_number),stat=alloc_stat)
     Call alloc_err("Domains",alloc_stat)
     
     Allocate(Domain_stats(Domain_number),stat=alloc_stat)
     Call alloc_err("Domain_stats",alloc_stat)
     
     Allocate(domain_path(0:Domain_number))
     domain_path = ''

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
          Domain_number * 24*24, Domain_number * 24*24, Domain_number        , &
          Domain_number *  6*24, Domain_number *  6*24, Domain_number *  6* 6, &
          Domain_number        , Domain_number *  6* 6, Domain_number        , &
          Domain_number        , Domain_number *     3, Domain_number *     9, &
          Domain_number *  6* 6, Domain_number        , Domain_number *     3, &
          Domain_number *     9, Domain_number *  6* 6, Domain_number           ], &
          branch = res)
     
     res%leaves(:)%pstat = -1
     
     Call set_bounds_in_branch(root, root%streams)
     
     !End If
     
     !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     If (out_amount == "DEBUG") THEN 
        Write(un_lf,fmt_dbg_sep)
        Write(un_lf,'(A)')"root right before serialisation"
        Call log_tree(root,un_lf,.True.)
        Write(un_lf,fmt_dbg_sep)
        flush(un_lf)
     END If
     !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     !** serialize root branch ************************************************
     Call serialize_branch(root,serial_root,serial_root_size,.TRUE.)
     
     !** Init restart and result files ****************************************
     aun = give_new_unit()
     
     If (restart == "Y") then

        Open(aun, file=trim(out%path)//"/"//trim(out%bsnm)//"_Activity.raw", &
             action="read", status="old", access="stream")

        read(aun)Domain_stats

        close(aun)
        
        Open(aun, file=trim(out%path)//"/"//trim(out%bsnm)//"_Activity.raw", &
             action="write", status="old", access="stream")
        
     Else

        !** Init State tracker file *******************************************
        Open(aun, file=trim(out%path)//"/"//trim(out%bsnm)//"_Activity.raw", &
             action="write", status="replace", access="stream")
        
        Domain_stats = 1
        write(aun)Domain_stats
        flush(aun)

     End If

     !** Init Domain Cross Reference and domain paths *************************
     dc = 0
     domain_path(0) = in%p_n_bsnm//"_results"

     nn = 1
     Do kk = xa_d(3), xe_d(3)
        Do jj = xa_d(2), xe_d(2)
           Do ii = xa_d(1), xe_d(1)
              
              dc         = dc + 1_mpi_ik
              path_count = dc / 2_mpi_ik
              Write(domain_path(dc),'(A,"/",I0)')Trim(domain_path(path_count)),dc

              Domains(nn) = ii + jj * nn_D(1)         + &
                                 kk * nn_D(1)*nn_D(2)
              nn = nn + 1_mpi_ik

           End Do
        End Do
     End Do

     !** Generate Activity_List ***********************************************
     Allocate(Activity(size_mpi-1),stat=alloc_stat)
     Call alloc_err("Activity_List",alloc_stat)

     Activity=1

     Allocate(act_domains(size_mpi-1),stat=alloc_stat)
     Call alloc_err("act_domains",alloc_stat)

     act_domains = 0
     
     If ( (Domain_Number*parts < size_mpi-1) .OR. &
          (count( Domain_stats < 10 ) == 0) .OR. &
          ((count( Domain_stats < 10 )*parts) < size_mpi-1) ) Then

        If (Domain_Number*parts < size_mpi-1) then
           Write(un_mon,FMT_ERR_A)"Domain_Number < size_mpi-1"
        Else If ((count( Domain_stats < 10 )*parts) < size_mpi-1) then
           Write(un_mon,FMT_ERR_A)"Remaining Domain_Number < Number of Solution Master"
           Write(un_mon,FMT_ERR_AI0)"Remaining Domain_Number    :", count( Domain_stats < 10 )
           Write(un_mon,FMT_ERR_AI0)"Number of solution masters :", (size_mpi-1)/parts
        Else
           Write(un_mon,FMT_ERR_A)"Restart on fully finished job"
        End If
        Write(un_mon,FMT_ERR_A)"This case is not supported"
        
        Call mpi_bcast(pro_path, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik,&
                       MPI_COMM_WORLD, ierr)

        Call mpi_bcast(pro_name, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik,&
             MPI_COMM_WORLD, ierr)
        
        !** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
        Call mpi_bcast(-1_ik, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik,&
                       MPI_COMM_WORLD, ierr)

        Call End_Timer("Init Process")
 
        Goto 1000

     End If

     Call End_Timer("Init Process")

     !*************************************************************************
     !** Start Workers ********************************************************

     Call Start_Timer("Broadcast Init Params")
     
     Call mpi_bcast(pro_path, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

     Call mpi_bcast(pro_name, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

     write(un_lf,FMT_MSG_AI0)"Broadcasting serialized root of size [Byte] ", serial_root_size*8
     
     Call mpi_bcast(serial_root_size, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)
     
     Call mpi_bcast(serial_root, INT(serial_root_size,mpi_ik), MPI_INTEGER8, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)

     Call End_Timer("Broadcast Init Params")

     !** Execute collective mpi_comm_split. Since mpi_comm_world rank 0 is
     !** the head master worker_comm is not needed and it should not be in
     !** any worker group and communicator. With MPI_UNDEFINED passed as
     !** color worker_comm gets the value MPI_COMM_NULL
     Call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, &
          rank_mpi, worker_comm, ierr)
     CALL handle_err(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", INT(ierr, KIND=ik))
    
  !****************************************************************************
  !** Ranks > 0 -- Worker slaves **********************************************
  !****************************************************************************
  Else
     
     !** Broadcast recieve init parameters ************************************
     Call Start_Timer("Broadcast Init Params")
     
     Call mpi_bcast(out%path         , INT(mcl,mpi_ik), MPI_CHAR    , 0_mpi_ik, MPI_COMM_WORLD, ierr)
     Call mpi_bcast(out%bsnm    , INT(mcl,mpi_ik), MPI_CHAR    , 0_mpi_ik, MPI_COMM_WORLD, ierr)
     Call mpi_bcast(serial_root_size, 1_mpi_ik       , MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)

     !** Serial_root_size == -1 ==> Signal that Domain_Number < size_mpi-1 ****
     If ( serial_root_size == -1 ) then

        Call End_Timer("Broadcast Init Params")
        
        Goto 1001

     End If
     
     Allocate(serial_root(serial_root_size))
     
     Call mpi_bcast(serial_root, INT(serial_root_size,mpi_ik), MPI_INTEGER8, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)

     Call End_Timer("Broadcast Init Params")

     !** Deserialize root branch **********************************************
     Call Start_Timer("Deserialize root branch")

     Call deserialize_branch(root, serial_root, .TRUE.)
     Call assign_pd_root (root)
     Call set_bounds_in_branch(root,root%streams)

     Call pd_get(root%branches(1),"Output amount", char_arr)
     out_amount = char_to_str(char_arr)
     deallocate(char_arr)

     Call pd_get(root%branches(1),"Restart", char_arr)
     restart = char_to_str(char_arr)
     deallocate(char_arr)

     Call End_Timer("Deserialize root branch")

     !** Init Domain Cross Reference ******************************************
     Call pd_get(root%branches(1),"Lower bounds of selected domain range",xa_d,3)
     Call pd_get(root%branches(1),"Upper bounds of selected domain range",xe_d,3)
     
     Domain_number = (xe_d(1)-xa_d(1)+1) * &
                     (xe_d(2)-xa_d(2)+1) * &
                     (xe_d(3)-xa_d(3)+1)

     
     Allocate(Domains(Domain_number),stat=alloc_stat)
     Call alloc_err("Domains",alloc_stat)

     Call pd_get(root%branches(2),"nn_D",nn_D)
     
     nn = 1
     Do kk = xa_d(3), xe_d(3)
        Do jj = xa_d(2), xe_d(2)
           Do ii = xa_d(1), xe_d(1)
              
              Domains(nn) = ii + jj * nn_D(1)         + &
                                 kk * nn_D(1)*nn_D(2)
              nn = nn + 1_mpi_ik

           End Do
        End Do
     End Do

     Call pd_get(root%branches(1),"No of mesh parts per subdomain",parts)

     !*************************************************************************
     !** All Worker Ranks -- Init worker Communicators ************************
     !*************************************************************************
     Call MPI_Comm_split(MPI_COMM_WORLD, Int((rank_mpi-1)/parts,mpi_ik), &
                         rank_mpi, worker_comm, ierr)
     CALL handle_err(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", INT(ierr, KIND=ik))
     
     Call MPI_COMM_RANK(WORKER_COMM, worker_rank_mpi, ierr)
     CALL handle_err(std_out, "MPI_COMM_RANK couldn't retrieve worker_rank_mpi", INT(ierr, KIND=ik))

 
     Call MPI_COMM_SIZE(WORKER_COMM, worker_size_mpi, ierr)
     CALL handle_err(std_out, "MPI_COMM_SIZE couldn't retrieve worker_size_mpi", INT(ierr, KIND=ik))

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

     !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     If (out_amount == "DEBUG") THEN
        write(un_mon,*)"On rank zero, stream sizes: ",dsize
     End If
     !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
     nn = 1
     ii = 1

     !** Supply all worker masters  with their first work package *************
     !** ii is incremented by ii = ii + parts                     *************
     Do While (ii <= (size_mpi-1_mpi_ik))

        if (nn > Domain_number) exit
        
        If ( Domain_stats(nn) /= 1 ) then
           nn = nn + 1_mpi_ik
           cycle
        End If

        act_domains(ii) = nn

        Do jj = ii, ii + parts-1
           
           !** Activity = 1 (Set above during init) ***
           Call mpi_send(Activity(jj), 1_mpi_ik, mpi_integer, Int(jj,mpi_ik), Int(jj,mpi_ik), &
             MPI_COMM_WORLD,ierr)
           CALL handle_err(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

           Call mpi_send(nn, 1_mpi_ik, mpi_integer8, Int(jj,mpi_ik), Int(jj,mpi_ik), &
                MPI_COMM_WORLD,ierr)
           CALL handle_err(std_out, "MPI_SEND of Domain number didn't succeed", INT(ierr, KIND=ik))
           
           if (out_amount /= "PRODUCTION") then
              Call mpi_send(domain_path(nn), Int(4_mpi_ik*mcl,mpi_ik), &
                   MPI_CHARACTER, Int(jj,mpi_ik), Int(jj,mpi_ik), MPI_COMM_WORLD,ierr)
              CALL handle_err(std_out, "MPI_SEND of Domain path didn't succeed", INT(ierr, KIND=ik))
           End if
        End Do
        
        !** Log to global stdout **********************************************
        Write(un_mon,'(2(A,I10))')"MPI rank : ",ii, &
                                  " ; Domain number : ",Domains(nn)
        flush(un_mon)
        
        nn = nn + 1_mpi_ik
        
        Call MPI_IRECV(Activity(ii), 1_mpi_ik, MPI_INTEGER, Int(ii,mpi_ik), Int(ii,mpi_ik), &
             MPI_COMM_WORLD, REQ_LIST(ii), IERR)
        CALL handle_err(std_out, "MPI_IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

        ii = ii + Int(parts,mpi_ik)

     End Do

     If ( restart == "N" ) then
        !** Add leaf with analyzed cube numbers to root ***********************
        Call add_leaf_to_branch(root, "Domain Numbers", Domain_number, Domains)
        Call set_bounds_in_branch(root, root%streams)
     End If
     
     !** Write Root header and input parameters *******************************
     Call Start_Timer("Write Root Branch")
     call store_parallel_branch(root, FH_MPI)
     Call Write_Tree(root)
     Call End_Timer("Write Root Branch")
     
     Call MPI_WAITANY(size_mpi-1_mpi_ik,req_list,finished,status_mpi,ierr)
     CALL handle_err(std_out, "MPI_WAITANY on req_list for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

     ii = finished
     Domain_stats(act_domains(ii)) = Activity(ii)
     write(aun,pos=(act_domains(ii)-1)*8+1) Activity(ii)
     flush(aun)
     act_domains(ii) = nn
     
     Do While (nn <= Domain_number)

        If ( Domain_stats(nn) /= 1 ) then
           if (nn > Domain_number) exit
           nn = nn + 1_mpi_ik
           cycle
        End If

        Do jj = ii, ii + parts-1

           Activity(jj) = 1_mpi_ik
           
           Call mpi_send(Activity(jj), 1_mpi_ik, mpi_integer, Int(jj,mpi_ik), Int(jj,mpi_ik), &
             MPI_COMM_WORLD,ierr)
           CALL handle_err(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))

           Call mpi_send(nn, 1_mpi_ik, mpi_integer8, Int(jj,mpi_ik), Int(jj,mpi_ik), &
                MPI_COMM_WORLD,ierr)
           CALL handle_err(std_out, "MPI_SEND of Domain number didn't succeed", INT(ierr, KIND=ik))

           if (out_amount /= "PRODUCTION") then
              Call mpi_send(domain_path(nn), Int(4_mpi_ik*mcl,mpi_ik), &
                   MPI_CHARACTER, Int(jj,mpi_ik), Int(jj,mpi_ik), MPI_COMM_WORLD,ierr)
              CALL handle_err(std_out, "MPI_SEND of Domain path didn't succeed", INT(ierr, KIND=ik))
           End if
        End Do
        
        !** Log to global stdout **********************************************
        Write(un_mon,'(2(A,I10))')"MPI rank : ",ii, " ; Domain number : ",Domains(nn)
        flush(un_mon)
        
        nn = nn + 1_mpi_ik
        
        Call MPI_IRECV(Activity(ii), 1_mpi_ik, MPI_INTEGER, Int(ii,mpi_ik), &
                       Int(ii,mpi_ik), MPI_COMM_WORLD, REQ_LIST(ii), IERR)
        CALL handle_err(std_out, "MPI_IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))


        Call MPI_WAITANY(size_mpi-1_mpi_ik,req_list, finished, status_mpi, ierr)
        CALL handle_err(std_out, "MPI_WAITANY on req_list for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

        ii = finished

        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        If (out_amount == "DEBUG") THEN
           Write(un_mon,*)"Domain ",ii, Domains(act_domains(ii))," finished"
        End If
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        Domain_stats(act_domains(ii)) = Activity(ii)
        write(aun,pos=(act_domains(ii)-1)*8+1) Activity(ii)
        flush(aun)
        act_domains(ii) = nn

     End Do

     !** Write last data element to ensure correct file size ******************
     Call MPI_FILE_WRITE_AT(FH_MPI(5), &
         Int(res%leaves(18)%lbound-1+Domain_number, MPI_OFFSET_KIND), &
         0_pd_rk, Int(1,pd_mpi_ik), MPI_Real8, status_mpi, ierr)
     
     Call MPI_WAITALL(size_mpi-1_mpi_ik, req_list, statuses_mpi, ierr)
     CALL handle_err(std_out, "MPI_WAITANY on req_list for IRECV of Activity(ii) didn't succeed", INT(ierr, KIND=ik))

     !** TODO refactor domain_cross reference from size size_mpi-1 to *********
     !** (size_mpi-1)/parts ***************************************************
     Do ii = 1, size_mpi-1, parts
        write(aun,pos=(act_domains(ii)-1)*8+1) Activity(ii)
     End Do
     
     flush(aun)

     Activity = -1
     
     Do ii = 1_mpi_ik, size_mpi-1_mpi_ik
        Call mpi_send(Activity(ii), 1_mpi_ik, mpi_integer, Int(ii,mpi_ik), &
                      Int(ii,mpi_ik), MPI_COMM_WORLD,ierr)
        CALL handle_err(std_out, "MPI_SEND of activity didn't succeed", INT(ierr, KIND=ik))
     End Do
   
  !****************************************************************************
  !** Ranks > 0 -- Workers ****************************************************
  !****************************************************************************
  Else 

     !** Extend out%bsnm and out%path with rank number ******************
     Write(out%path,'(A,A,I7.7,A)')Trim(out%path),"Rank_",rank_mpi,"/"
     Write(out%bsnm,'(A,A,I7.7)')Trim(out%bsnm),"_",rank_mpi
     
     !** Prepare Rank output directory ****************************************
     c_char_array(1:len(Trim(out%path)//Char(0))) = str_to_char(Trim(out%path)//Char(0))

     Call Stat_Dir(c_char_array, stat_c_int)

     If ( stat_c_int /= 0 ) Then

        Call execute_command_line("mkdir -p "//trim(out%path),CMDSTAT=stat)

        If ( stat /= 0 ) Then
           Write(un_mon,*)"Could not execute syscall"
           Write(un_mon,*)"mkpir -p "//trim(out%path)
           Write(un_mon,*)"Program halted"
           Stop
        End If

        Call Stat_Dir(c_char_array, stat_c_int)

        If ( stat_c_int /= 0 ) Then
           Write(un_mon,*)"Could not create directory"
           Write(un_mon,*)trim(out%path)
           Write(un_mon,*)"Program halted"
           Stop
        End If
        
     Else If ( (stat_c_int == 0) .AND. (Restart == "Y") ) Then

        Write(un_mon,FMT_MSG_A)"Reusing the output directory"
        Write(un_mon,FMT_MSG_A)trim(out%path)
        
     Else
        
        Write(std_out,FMT_ERR_A)"The output directory"
        Write(std_out,FMT_ERR_A)trim(out%path)
        Write(std_out,FMT_ERR_A)"apparently exists already with restart not equal to Y !!!"
        Write(std_out,FMT_ERR_A)"Please check your struct-process-parameters.sh file   !!!"
        write(std_out,FMT_STOP)
        Goto 1001
        
     End If

     Call link_start(link_name,.True.,.True.)

     !** Worker Loop **********************************************************
     Do

        Call mpi_recv(Active, 1_mpi_ik, mpi_integer, 0_mpi_ik, rank_mpi, &
                      MPI_COMM_WORLD, status_mpi, ierr)
        CALL handle_err(std_out, "MPI_RECV on Active didn't succseed", INT(ierr, KIND=ik))

        If (Active == -1) Exit

        CALL mpi_recv(nn, 1_mpi_ik, mpi_integer8, 0_mpi_ik, rank_mpi, &
                      MPI_COMM_WORLD, status_mpi, ierr)
        CALL handle_err(std_out, "MPI_RECV on Domain didn't succeed", INT(ierr, KIND=ik))

        Domain = Domains(nn)
        
        if (out_amount /= "PRODUCTION") then
           !** >> Recieve Job_Dir << ******************************************
           CALL mpi_recv(job_dir, 4_mpi_ik*int(mcl,mpi_ik), mpi_character, 0_mpi_ik, &
                rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
           CALL handle_err(std_out, "MPI_RECV on Domain path didn't succeed", INT(ierr, KIND=ik))

        Else
           job_dir = out%path

        End if
        
        if (job_dir(len(job_dir):len(job_dir)) /= "/") then
           job_dir = trim(job_dir)//"/"
        End if

        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        If (out_amount == "DEBUG") THEN
           Write(un_lf, fmt_dbg_sep)
           Write(un_lf, fmt_MSG_AI0)"Root pointer before exec_single_domain on proc",rank_mpi
           Call log_tree(root,un_lf,.True.)
           Write(un_lf, fmt_dbg_sep)
        END If
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !======================================================================
        Call exec_single_domain(root, nn, Domain, job_dir, Active, fh_mpi, &
             worker_rank_mpi, worker_size_mpi, worker_comm)
        !======================================================================
        
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        If (out_amount == "DEBUG") THEN
           Write(un_lf, fmt_dbg_sep)
           Write(un_lf, fmt_MSG_AI0)"Root pointer after exec_single_domain on proc",rank_mpi
           Call log_tree(root,un_lf,.True.)
           Write(un_lf, fmt_dbg_sep)
        END If
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !======================================================================
        !== Organize Results
        !======================================================================

        !** Look for the Domain branch ****************************************
        domain_desc=''
        Write(domain_desc, '(A,I0)')'Domain ',Domain
        
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
!!$                 !write(tmp_fn,fmt_filename) trim(job_dir)//trim(out%bsnm)//'_',&
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
!!$        Write(mesh_desc,'(A,I0)')'Mesh info of '//trim(out%bsnm)//'_',Domain
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

        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           
        Call MPI_ISEND(Active, 1_mpi_ik, MPI_INTEGER, 0_mpi_ik, rank_mpi, &
             MPI_COMM_WORLD, REQUEST, IERR)
        CALL handle_err(std_out, "MPI_ISEND on Active didn't succeed", INT(ierr, KIND=ik))

        Call MPI_WAIT(REQUEST, status_mpi, ierr)
        CALL handle_err(std_out, "MPI_WAIT on request for ISEND Active didn't succeed", INT(ierr, KIND=ik))

     End Do

     call PetscFinalize(petsc_ierr)
     
  End If

1000 Continue
  !============================================================================
  
  !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  If (out_amount == "DEBUG") THEN
     Write(un_lf,fmt_dbg_sep)
     Write(un_lf,fmt_MSG_AI0)"Final Root pointer proc",rank_mpi
     Call log_tree(root,un_lf,.True.)
     Write(un_lf,fmt_dbg_sep)
  END If
  !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  Call link_end(link_name,.True.)

1001 Continue



IF(rank_mpi == 0) THEN
   !------------------------------------------------------------------------------
   ! Assign out=in and define the app-name
   ! 1. Call Close the meta file and decide whether to alter the basename.
   ! 2. Call Close the log  file.
   !------------------------------------------------------------------------------
   ! Hardcode a program specific altered meta app name
   ! This segment must be considered best practice, as it does not 
   ! hide the principle functionality of the setup.
   !------------------------------------------------------------------------------
   out = in
   out%app = 'ddtc' 
   CALL meta_close    (m_rry)
   CALL meta_add_ascii(fh=un_lf, suf=log_suf, st='stop')
   ! CALL meta_add_ascii(fh=fhr, suf=res_suf, st='stop')
END IF ! (rank_mpi == 0)

Call MPI_FINALIZE(ierr)
CALL handle_err(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, KIND=ik))

End Program main_struct_process
