!==============================================================================
!> \file mod_puredat.f90
!> Modules of the data handling puredat library.
!>
!> The file holds modules which are necessary to use the puredat data format
!>
!> \author Ralf Schneider
!> \date 22.01.2010
!>
!> \todo Character stream copy in and copy out von strings realisieren
!> \todo Warning in log_tree wenn pointer of leaf nicht korrekt connected

!==============================================================================
!> Global integer and real kinds for the puredat data handling library
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat_precision

  Implicit none
  
  Integer, Parameter :: pd_ik = 8         !** Puredat Integer kind parameter
  Integer, Parameter :: pd_rk = 8         !** Puredat Real    kind parameter
  Integer, Parameter :: pd_mpi_ik = 4     !** Puredat Integer MPI kind parameter
  
End Module puredat_precision

!==============================================================================
!> Global constants and parameters for the puredat data handling library
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat_constants

  Implicit none
  
  !> Number of currently used stream variables in puredat_streams
  !>
  !> The total number of currently used stream variables which is the  
  !> number of arrays defined in puredat_streams independently from 
  !> their data type
  Integer, Parameter :: no_streams = 7

  !> Maximum character length used in puredat library
  Integer, Parameter :: pd_mcl = 512
  !> Maximum Character Length in pd_ik elements
  Integer, Parameter :: pd_ce  = 512/8

End Module puredat_constants

!==============================================================================
!> Global variables for the puredat data handling library
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat_globals
  
  Use puredat_constants

  Implicit None

  !============================================================================
  !== files and paths
  !> puredat project path
  !>
  !> Path to puredat project files which means header-, stream- and 
  !> log-files
  Character(len=pd_mcl) :: pro_path
  !> puredat project name
  !>
  !> Base name of the puredat project files which are ...
  Character(len=pd_mcl) :: pro_name
  
  !> puredat monitor unit
  Integer               :: pd_umon  != OUTPUT_UNIT

End Module puredat_globals

!==============================================================================
!> Derived datatypes for puredat data handling
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat_types

  use puredat_precision
  use puredat_constants

  implicit none

  !============================================================================
  !> Type: Stream allocatable stream arrays
  !>
  !> Type for the stream arrays within a tBranch structure. <BR>
  !> The type holds also the global informations about the streams and the 
  !>stream files
  Type tStreams

     !> Stream dimensions
     !>
     !> Array of dimension puredat_constants::no_streams which holds the 
     !> dimensions of the puredat streams. If a stream is not in use the
     !> corresponding element of dim_st = 0
     Integer(Kind=pd_ik)  , Dimension(no_streams) :: dim_st = 0

     !> Stream position
     !>
     !> Array of dimension puredat_constants::no_streams which holds the 
     !> current position in global stream coordinates of the streams.
     !> If a stream is not in use the corresponding element of ii_st = 0
     Integer(Kind=pd_ik)  , Dimension(no_streams) :: ii_st

     Integer(Kind=pd_ik)                          :: no_leaves    = 0
     Integer(Kind=pd_ik)                          :: no_branches  = 0

     Character(len=pd_mcl), Dimension(no_streams) :: stream_files = ""
     Logical              , Dimension(no_streams) :: ifopen       = .FALSE.
     Integer              , Dimension(no_streams) :: units        = -1

     !> Global stream variables
     Integer(kind=1)      , Dimension(:), pointer :: int1_st  => null()
     Integer(kind=2)      , Dimension(:), pointer :: int2_st  => null()
     Integer(kind=4)      , Dimension(:), pointer :: int4_st  => null()
     Integer(kind=8)      , Dimension(:), pointer :: int8_st  => null()
     !
     Real(kind=8)         , Dimension(:), pointer :: real8_st => null()
     !
     Character            , Dimension(:), pointer :: char_st  => null()
     !> Global stream variable for logical data
     Logical(Kind=1)      , Dimension(:), pointer :: log_st   => null()

  End type tStreams

  !============================================================================
  !> Type: Stream chunk pointer
  !>
  !> Type for the leafs in the puredat data description tree. <BR>
  !> The type holds the informations of the actual data chunk position within 
  !> the stream arrays.
  Type tLeaf

     !> ASCII description of what the data are
     Character(Len=pd_mcl)             :: desc
     !> Number of data
     Integer(Kind=pd_ik)               :: dat_no
     !> Data Type
     !>
     !> Currently the reference is <br>
     !> 1: 1 Byte Integer data <BR>
     !> 2: 2 Byte Integer data <BR>
     !> 3: 4 Byte Integer data <BR>
     !> 4: 8 Byte Integer data <BR>
     !> 5: 8 Byte Floating point data <BR>
     !> 6: 1 Byte Character data <BR>
     !> 7: Logical data <BR>
     Integer(Kind=1)                   :: dat_ty

     !> Lower bound index in stream array
     Integer(Kind=pd_ik)               :: lbound

     !> Upper bound index in stream array
     Integer(Kind=pd_ik)               :: ubound

     Integer(kind=pd_ik)               :: pstat = 0

     !> Data chunk pointer to x Byte y data
     Integer(Kind=1), Dimension(:),Pointer   :: p_int1  => null()
     Integer(Kind=2), Dimension(:),Pointer   :: p_int2  => null()
     Integer(Kind=4), Dimension(:),Pointer   :: p_int4  => null()
     Integer(Kind=8), Dimension(:),Pointer   :: p_int8  => null()

     Real   (Kind=8), Dimension(:),Pointer   :: p_real8 => null()
     
     Character      , Dimension(:),Pointer   :: p_char  => null()
     
     Logical(Kind=1), Dimension(:),Pointer   :: p_log   => null()

  End Type tLeaf

  !============================================================================
  !> Type: Branches in the puredat data structure tree
  !>
  !> Recursive data type for branches of the puredat data structure tree.
  !> It can contain leaves as well as other branches
  Type tBranch

     !> ASCII description of what the data are
     Character(Len=pd_mcl)                   :: desc 
     !> Number of childdren of type tBranch
     Integer(Kind=pd_ik)                     :: no_branches = 0
     !> Number of childdren of type tLeaf
     Integer(Kind=pd_ik)                     :: no_leaves = 0

     !> Pointer to branch children
     Type(tBranch), Dimension(:), Pointer    :: branches => null()
     !> Pointer to leaf children
     Type(tLeaf)  , Dimension(:), Pointer    :: leaves   => null()

     !> Streams
     Type(tStreams), Allocatable             :: streams

  End Type tBranch

End Module puredat_types

!==============================================================================
!> Variables, Functions and subroutines for global puredat data communication
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat_com

  use puredat_globals
  use puredat_types  

  implicit none

  Type(tBranch), Pointer, Private, Save :: pd_root => null()
  Logical               , Private, Save :: root_assigned = .FALSE.

contains

  !============================================================================
  !> Subroutine which assigns the global puredat root pointer
  Function pd_root_assigned() Result(assigned)

    Logical :: assigned

    assigned = root_assigned
    
  End Function pd_root_assigned

  
  !============================================================================
  !> Subroutine which assigns the global puredat root pointer
  Subroutine assign_pd_root(rb)

    Type(tBranch), target, intent(in) :: rb

    if (.not. root_assigned) then
    
       pd_root => rb
       root_assigned = .TRUE.
       !Write(pd_umon,PDF_M_A)"Assigned  pd_root to "//trim(pd_root%desc)
       
    else

       Write(pd_umon,*)"!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!-!-!"
       Write(pd_umon,*)"The root pointer is allready assigned to a branch with description !!"
       Write(pd_umon,*)trim(pd_root%desc)
       Write(pd_umon,*)"No reassignment is done                                            !!"
       Write(pd_umon,*)"!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!-!-!"
       
       !Write(pd_umon,*)"This routine can only be used once !"
       !Write(pd_umon,*)"Please check your implementation !!!"
       !Write(pd_umon,*)"Program halted !!!!!!!!!!!!!!!!!!!!!"

    End if

  End Subroutine assign_pd_root
  
  !============================================================================
  !> Subroutine which retrieves the global puredat root pointer
  Subroutine get_pd_root(rb)

    Type(tBranch), pointer, intent(out) :: rb

    if (root_assigned) then
    
       rb => pd_root
       !Write(pd_umon,PDF_M_A)"Retrieved pd_root as "//trim(pd_root%desc)

    else

       Write(pd_umon,*)"!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!"
       Write(pd_umon,*)"The root pointer is not assigned !!!!"
       Write(pd_umon,*)"Please check your implementation !!!!"
       Write(pd_umon,*)"Program halted                   !!!!"
       Write(pd_umon,*)"!-!-!-!-!-!-!-!-!_!-!-!-!-!-!-!-!-!-!"
       Stop
       
    End if

  End Subroutine get_pd_root

  !============================================================================
  !> Subroutine which nullifies the puredat root pointer
  Subroutine nullify_pd_root()

    pd_root => null()
    root_assigned = .FALSE.
    
  End Subroutine nullify_pd_root

End Module puredat_com

!==============================================================================
!> Function and subroutines for puredat data handling
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module puredat

USE global_std
USE puredat_types
USE puredat_com
USE mpi

Implicit None

  !============================================================================
  !== Private routines
  Private alloc_error
  Private file_err
  Private char_to_str
  Private str_to_char
  Private copy_leaves_to_streams
  
  !============================================================================
  !== Interfaces
  !> Interface: pd_store
  !> \author Ralf Schneider
  !> \date 22.01.2010
  Interface pd_store

     Module Procedure pd_store_1
     Module Procedure pd_store_2
     Module Procedure pd_store_3
     Module Procedure pd_store_4
     Module Procedure pd_store_4_2D
     Module Procedure pd_store_5
     Module Procedure pd_store_5_2D
     Module Procedure pd_store_6
     Module Procedure pd_store_6_str

  End Interface pd_store

  !> Interface: add_leaf_to_branch
  !> \author Ralf Schneider
  !> \date 22.01.2010
  Interface add_leaf_to_branch

     Module Procedure add_empty_leaf_to_branch
     Module Procedure add_filled_leaf_to_branch
     Module Procedure add_filled_leaf_to_branch_4
     Module Procedure add_filled_leaf_to_branch_5
     Module Procedure add_filled_leaf_to_branch_6

  End Interface add_leaf_to_branch

  !> Interface: pd_load_leaf
  !> \author Ralf Schneider
  !> \date 08.06.2012
  Interface pd_load_leaf

     Module Procedure pd_load_leaf_1
     Module Procedure pd_load_leaf_2
     Module Procedure pd_load_leaf_3
     Module Procedure pd_load_leaf_4
     Module Procedure pd_load_leaf_5
     Module Procedure pd_load_leaf_6

     Module Procedure pd_load_leaf_4_2D
     Module Procedure pd_load_leaf_5_2D

  End Interface pd_load_leaf

  !> Interface: pd_get
  !> \author Ralf Schneider
  !> \date 08.06.2012
  Interface pd_get

     Module Procedure pd_get_4_scalar
     Module Procedure pd_get_5_scalar

     Module Procedure pd_get_4
     Module Procedure pd_get_4_vector
     Module Procedure pd_get_5
     Module Procedure pd_get_6

     Module Procedure pd_get_leaf_4_const
     Module Procedure pd_get_leaf_4_2D_const
     Module Procedure pd_get_leaf_5_const

     Module Procedure pd_get_leaf_pointer
     
  End Interface pd_get

  !> Interface: pd_read
  !> \author Ralf Schneider
  !> \date 08.06.2012
  Interface pd_read_leaf

     Module Procedure pd_read_leaf_1
     Module Procedure pd_read_leaf_2
     Module Procedure pd_read_leaf_3
     Module Procedure pd_read_leaf_4
     Module Procedure pd_read_leaf_4_2D
     Module Procedure pd_read_leaf_4_scal
     Module Procedure pd_read_leaf_5
     Module Procedure pd_read_leaf_5_2D
     Module Procedure pd_read_leaf_6

  End Interface pd_read_leaf

  !> Interface: open_stream_files
  !> \author Ralf Schneider
  !> \date 26.07.2012
  Interface open_stream_files

     Module Procedure open_stream_files_from_tree
     Module Procedure open_stream_files_from_streams
     Module Procedure open_stream_files_from_streams_mpi

  End Interface open_stream_files

  !> Interface: close_stream_files
  !> \author Ralf Schneider
  !> \date 26.07.2012
  Interface close_stream_files

     Module Procedure close_stream_files_from_tree
     Module Procedure close_stream_files_from_streams

  End Interface close_stream_files

  !> Interface: read_stream_files
  !> \author Ralf Schneider
  !> \date 26.07.2012
  Interface read_streams

     Module Procedure read_streams_from_branch
     Module Procedure read_streams_from_streams

  End Interface read_streams
  
  !> Interface: search_branch
  !> \author Ralf Schneider
  !> \date 05.10.2018
  Interface search_branch

     Module Procedure Search_branch_rec
     Module Procedure Search_branch_wrn

  End Interface search_branch
  
  !============================================================================
  !== Character constants
  !> Character constant for branch separation in the puredat ascii header file
  Character(Len=*), Parameter :: fmt_bsep="('<==branch==>')"
  !> Character constant for leaf separation in the puredat ascii header file
  Character(Len=*), Parameter :: fmt_lsep="('<--leaf-->')"

Contains

  !============================================================================
  !> Subroutine that initializes the output to std-out according to the
  !> environment variable PRD_STDOUT
  Subroutine pd_init_std_out()

    Character(len=pd_mcl) :: env_var
    Logical            :: opened, exist
    
    Call get_environment_Variable("PRD_STDOUT", env_var)

    !** If PRD_STDOUT is set to a value ***************************************
    If (Len_Trim(env_var) > 0) then

       Inquire(file=trim(env_var), opened = opened)
       Inquire(file=trim(env_var), exist  = exist)

       IF (exist .AND. opened) Then
          
          Inquire(file=trim(env_var), number=pd_umon)

       Else if (exist .AND. (.NOT.opened)) then

          pd_umon = pd_give_new_unit()
          Open(unit=pd_umon,  file=trim(env_var), Action='Write', &
               status='old', position='Append' )

       Else if (.NOT. exist) then

          pd_umon = pd_give_new_unit()
          Open(unit=pd_umon,  file=trim(env_var), Action='Write', &
               status='NEW' )
          
       End IF
       
    End If
   
  End Subroutine pd_init_std_out
  
  !============================================================================
  !> \name Puredat operators
  !> @{
  !> Module procedures for operator definition for the puredat types
  
  !============================================================================
  !> Subroutine which assigns two tLeaf structures by full copy
  Subroutine Assign_leaves(leavesl, leavesr)

    Type(tLeaf), Intent(out), Dimension(:) :: leavesl
    Type(tLeaf), Intent(in) , Dimension(:) :: leavesr

    Integer(kind=pd_ik)                      :: ii
    
    If (size(leavesl) /= size(leavesr)) then

       Write(pd_umon, PDF_E_AI0)"In 'Assign_leaves' size(leavesl) = ",size(leavesl), &
                                " /= size(leavesr) = ",size(leavesr)
       Write(pd_umon, PDF_E_A  )"No assignment is done !!"
       Goto 1000
       
    End If

    Do ii = 1, size(leavesl)

       leavesl(ii)%desc   = leavesr(ii)%desc
       leavesl(ii)%dat_no = leavesr(ii)%dat_no
       leavesl(ii)%dat_ty = leavesr(ii)%dat_ty
       leavesl(ii)%lbound = leavesr(ii)%lbound
       leavesl(ii)%ubound = leavesr(ii)%ubound

       If (leavesl(ii)%dat_no>0) then

          leavesl(ii)%pstat = 1

          Select Case ( leavesl(ii)%dat_ty )
          Case(1)
             Allocate(leavesl(ii)%p_int1(leavesl(ii)%dat_no))
             leavesl(ii)%p_int1 = leavesr(ii)%p_int1
          Case(2)
             Allocate(leavesl(ii)%p_int2(leavesl(ii)%dat_no))
             leavesl(ii)%p_int2 = leavesr(ii)%p_int2
          Case(3)
             Allocate(leavesl(ii)%p_int4(leavesl(ii)%dat_no))
             leavesl(ii)%p_int4 = leavesr(ii)%p_int4
          Case(4)
             Allocate(leavesl(ii)%p_int8(leavesl(ii)%dat_no))
             leavesl(ii)%p_int8 = leavesr(ii)%p_int8
          Case(5)
             Allocate(leavesl(ii)%p_real8(leavesl(ii)%dat_no))
             leavesl(ii)%p_real8 = leavesr(ii)%p_real8
          Case(6)
             Allocate(leavesl(ii)%p_char(leavesl(ii)%dat_no))
             leavesl(ii)%p_char = leavesr(ii)%p_char
          Case(7)
             Allocate(leavesl(ii)%p_log(leavesl(ii)%dat_no))
             leavesl(ii)%p_log = leavesr(ii)%p_log

          End Select
          
       End If
       
    End Do
    
1000   continue

  End Subroutine Assign_leaves
  
  !============================================================================
  !> Subroutine which assigns two tBranch structures
  !>
  !> The routine DOES NOT perform a full copy of a recursive tBranch structure
  !> Only the atributes of the top level are copied from the right hand side to
  !> the left hand side the branches attribute is targeted by a pointer
  !> assignment.
  Subroutine Assign_branches(branchl, branchr)

    Type(tBranch), Intent(out)          :: branchl
    Type(tBranch), Intent(in)           :: branchr

    branchl%desc        =  branchr%desc
    branchl%no_branches =  branchr%no_branches
    branchl%no_leaves   =  branchr%no_leaves

    If (allocated(branchr%streams) .AND. allocated(branchl%streams)) then
       Write(pd_umon,*)"In Subroutine Assign_branches"
       Write(pd_umon,*)"Component streams is allocated for both left and"
       Write(pd_umon,*)"right hand side branches"
       write(pd_umon,*)"This usecase is not supported at the moment"
       Write(pd_umon,*)"PROGRAM STOPPED !!"
       stop
    End If

    If (allocated(branchr%streams)) then
       allocate(branchl%streams)
       branchl%streams = branchr%streams  
    End If

    If (branchr%no_branches > 0) then
       If (associated(branchl%branches)) Deallocate(branchl%branches)
       !allocate(branchl%branches(branchr%no_branches))
       branchl%branches    => branchr%branches
    Else
       nullify(branchl%branches)
    End If

    If (branchr%no_leaves > 0) then
       If (associated(branchl%leaves)) Deallocate(branchl%leaves)
       !Allocate(branchl%leaves(branchr%no_leaves))
       branchl%leaves      => branchr%leaves
    Else
       nullify(branchl%leaves)
    End If

  End Subroutine Assign_branches
  !> @} 
  != End of memeber group "Puredat operators"

  !############################################################################
  !> \name Puredat structure routines
  !> @{
  !> Routines for initalising, altering and cleaning up puredat data-location
  !> trees defined by tBranch and tLeaf

  !============================================================================
  !> Subroutine which returns the data size currently hold by all leafs in a
  !> tBranch structure
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_data_size_rec
  Subroutine get_data_size(br,dsize)

    Type(tBranch), Intent(in)                     :: br
    Integer(kind=pd_ik), intent(Out),Dimension(2) :: dsize

    dsize = 0

    Call get_data_size_rec(br, dsize)
    
    
  End Subroutine get_data_size
  
  !============================================================================
  !> Subroutine which returns the data size currently hold by all leafs in a
  !> tBranch structure
  recursive Subroutine get_data_size_rec(br,dsize)

    Type(tBranch), Intent(in)        :: br
    Integer(kind=pd_ik),Dimension(2) :: dsize

    Integer(kind=pd_ik)              :: ii
    
    Do ii = 1, br%no_leaves

       dsize(1) = dsize(1) + br%leaves(ii)%dat_no
       
       Select Case (br%leaves(ii)%dat_ty)

       Case (1)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no
       Case (2)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no*2_pd_ik
       Case (3)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no*4_pd_ik
       Case (4)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no*8_pd_ik
       Case (5)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no*8_pd_ik
       Case (6)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no
       Case (7)
          dsize(2) = dsize(2) + br%leaves(ii)%dat_no
       End Select

    End Do

    Do ii = 1, br%no_branches

       call get_data_size_rec(br%branches(ii),dsize)

    End Do
    
  End Subroutine get_data_size_rec

  !============================================================================
  !> Subroutine which returns the data size currently hold by all leafs in a
  !> tBranch structure
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_stream_size_rec
  Subroutine get_stream_size(br,dsize)

    Type(tBranch), Intent(in)                              :: br
    Integer(kind=pd_ik), intent(Out),Dimension(no_streams) :: dsize

    dsize = 0

    Call get_stream_size_rec(br, dsize)
    
    
  End Subroutine get_stream_size
  
  !============================================================================
  !> Subroutine which returns the data size currently hold by all leafs in a
  !> tBranch structure
  recursive Subroutine get_stream_size_rec(br,dsize)

    Type(tBranch), Intent(in)                 :: br
    Integer(kind=pd_ik),Dimension(no_streams) :: dsize

    Integer(kind=pd_ik)                       :: ii
    
    Do ii = 1, br%no_leaves
       dsize(br%leaves(ii)%dat_ty) = dsize(br%leaves(ii)%dat_ty) + br%leaves(ii)%dat_no     
    End Do

    Do ii = 1, br%no_branches

       call get_stream_size_rec(br%branches(ii),dsize)

    End Do
    
  End Subroutine get_stream_size_rec
  
  !============================================================================
  !> Subroutine which adds a leaf of type dat_ty to a tBranch structure 
  !>
  !> The subroutine adds a leaf to a tBranch structure at the end of its leaves
  !> array and fills the attributes of the new leaf desc, dat_ty and dat_no
  !> as given by the input parameters.<br>
  !> !! There is no data asignment done.<br>
  !> !! There is no assignment of the lbound and ubound attributes done in the
  !> new leaf.<br>
  !> !! The pstat attribute is set to 0.
  Subroutine add_filled_leaf_to_branch(t_b, desc, dat_ty, dat_no )

    Type(tBranch), Intent(inout)                  :: t_b

    Character(Len=*)   , Intent(in)               :: desc
    Integer(Kind=1)    , Intent(in)               :: dat_ty 
    Integer(Kind=pd_ik), Intent(in)               :: dat_no

    !> Pointer to leaf children
    Type(tLeaf), Dimension(:), Pointer   :: leaves

    Integer(kind=pd_ik)                  :: ii

    leaves     => t_b%leaves
    t_b%leaves => Null()

    allocate(t_b%leaves(t_b%no_leaves+1))

    Do ii = 1, t_b%no_leaves
       t_b%leaves(ii)=leaves(ii)
    End Do

    if (t_b%no_leaves > 0) deallocate(leaves)

    t_b%no_leaves = t_b%no_leaves + 1

    t_b%leaves(t_b%no_leaves)%desc   = trim(desc)
    t_b%leaves(t_b%no_leaves)%dat_ty = dat_ty
    t_b%leaves(t_b%no_leaves)%dat_no = dat_no
    t_b%leaves(t_b%no_leaves)%lbound = 0_pd_ik
    t_b%leaves(t_b%no_leaves)%ubound = 0_pd_ik

    t_b%leaves(t_b%no_leaves)%p_int1 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int2 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int4 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_real8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_char => NULL()
    t_b%leaves(t_b%no_leaves)%p_log  => NULL()

    t_b%leaves(t_b%no_leaves)%pstat = 0

  End Subroutine add_filled_leaf_to_branch

  !============================================================================
  !> Subroutine which adds a leaf to a tBranch structure
  !>
  !> The subroutine adds a leaf of type 4 (Int_8) to a tBranch structure at 
  !> the end of its leaves array and fills the attributes of the new leaf 
  !> desc and dat_no as given by the input parameters.<br>
  !> The p_int8 data pointer is allocated with size dat_no and the values array
  !> passed in as the last parameter is copied to the allocated memory. <br>
  !> !! There is no assignment of the lbound and ubound attributes done in the
  !> new leaf as the data values reside directly in the tLeaf structure and not
  !> in a tStreams structure belonging to one of the Leaf's parent tBranch
  !> structures. <br>
  !> !! The fourth element of the pstat attribute arrayis set to 1 to indicate,
  !> that the data resides directly in the leaf.
  Subroutine add_filled_leaf_to_branch_4(t_b, desc, dat_no, values)

    Type(tBranch), Intent(inout)                  :: t_b

    Character(Len=*)   , Intent(in)               :: desc
    Integer(Kind=1)    , Parameter                :: dat_ty = 4
    Integer(Kind=pd_ik), Intent(in)               :: dat_no

    Integer(Kind=8)    , Intent(in), Dimension(:) :: values

    !> Pointer to leaf children
    Type(tLeaf), Dimension(:), Pointer   :: leaves

    Integer(kind=pd_ik)                  :: ii
    Integer                              :: alloc_stat

    leaves     => t_b%leaves
    t_b%leaves => Null()

    allocate(t_b%leaves(t_b%no_leaves+1))

    Do ii = 1, t_b%no_leaves
       t_b%leaves(ii)=leaves(ii)
    End Do

    if (t_b%no_leaves > 0) deallocate(leaves)

    t_b%no_leaves = t_b%no_leaves + 1

    t_b%leaves(t_b%no_leaves)%desc    = trim(desc)
    t_b%leaves(t_b%no_leaves)%dat_ty  = dat_ty
    t_b%leaves(t_b%no_leaves)%dat_no  = dat_no
    t_b%leaves(t_b%no_leaves)%lbound  = 0_pd_ik
    t_b%leaves(t_b%no_leaves)%ubound  = 0_pd_ik

    t_b%leaves(t_b%no_leaves)%p_int1  => NULL()
    t_b%leaves(t_b%no_leaves)%p_int2  => NULL()
    t_b%leaves(t_b%no_leaves)%p_int4  => NULL()
    t_b%leaves(t_b%no_leaves)%p_int8  => NULL()
    t_b%leaves(t_b%no_leaves)%p_real8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_char  => NULL()
    t_b%leaves(t_b%no_leaves)%p_log   => NULL()

    t_b%leaves(t_b%no_leaves)%pstat = 0

    allocate(t_b%leaves(t_b%no_leaves)%p_int8(dat_no), stat=alloc_stat)
    call alloc_error(alloc_stat, "t_b%leaves(t_b%no_leaves)%p_int8(dat_no)", &
         "add_filled_leaf_to_branch_4") 
    t_b%leaves(t_b%no_leaves)%p_int8 = values
    
    t_b%leaves(t_b%no_leaves)%pstat = 1

  End Subroutine add_filled_leaf_to_branch_4

  !============================================================================
  !> Subroutine which adds a leaf to a tBranch structure
  !>
  !> The subroutine adds a leaf of type 5 (Real8) to a tBranch structure at 
  !> the end of its leaves array and fills the attributes of the new leaf 
  !> desc and dat_no as given by the input parameters.<br>
  !> The p_int8 data pointer is allocated with size dat_no and the values array
  !> passed in as the last parameter is copied to the allocated memory. <br>
  !> !! There is no assignment of the lbound and ubound attributes done in the
  !> new leaf as the data values reside directly in the tLeaf structure and not
  !> in a tStreams structure belonging to one of the Leaf's parent tBranch
  !> structures. <br>
  !> !! The fifth element of the pstat attribute array is set to 1 to indicate,
  !> that the data resides directly in the leaf.
  Subroutine add_filled_leaf_to_branch_5(t_b, desc, dat_no, values)

    Type(tBranch), Intent(inout)                  :: t_b

    Character(Len=*)   , Intent(in)               :: desc
    Integer(Kind=1)    , Parameter                :: dat_ty=5
    Integer(Kind=pd_ik), Intent(in)               :: dat_no

    Real(Kind=8)       , Intent(in), Dimension(:) :: values

    !> Pointer to leaf children
    Type(tLeaf), Dimension(:), Pointer   :: leaves

    Integer(kind=pd_ik)                  :: ii
    Integer                              :: alloc_stat

    leaves     => t_b%leaves
    t_b%leaves => Null()

    allocate(t_b%leaves(t_b%no_leaves+1))

    Do ii = 1, t_b%no_leaves
       t_b%leaves(ii)=leaves(ii)
    End Do

    if (t_b%no_leaves > 0) deallocate(leaves)

    t_b%no_leaves = t_b%no_leaves + 1

    t_b%leaves(t_b%no_leaves)%desc   = trim(desc)
    t_b%leaves(t_b%no_leaves)%dat_ty = dat_ty
    t_b%leaves(t_b%no_leaves)%dat_no = dat_no
    t_b%leaves(t_b%no_leaves)%lbound = 0_pd_ik
    t_b%leaves(t_b%no_leaves)%ubound = 0_pd_ik

    t_b%leaves(t_b%no_leaves)%p_int1 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int2 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int4 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_real8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_char => NULL()
    t_b%leaves(t_b%no_leaves)%p_log  => NULL()

    t_b%leaves(t_b%no_leaves)%pstat = 0

    allocate(t_b%leaves(t_b%no_leaves)%p_real8(dat_no), stat=alloc_stat)
    call alloc_error(alloc_stat, "t_b%leaves(t_b%no_leaves)%p_int8(dat_no)", &
         "add_filled_leaf_to_branch_5") 
    t_b%leaves(t_b%no_leaves)%p_real8 = values
    
    t_b%leaves(t_b%no_leaves)%pstat = 1


  End Subroutine add_filled_leaf_to_branch_5

  !============================================================================
  !> Subroutine which adds a leaf to a tBranch structure
  !>
  !> The subroutine adds a leaf of type 6 (char) to a tBranch structure at 
  !> the end of its leaves array and fills the attributes of the new leaf 
  !> desc and dat_no as given by the input parameters.<br>
  !> The p_int8 data pointer is allocated with size dat_no and the values array
  !> passed in as the last parameter is copied to the allocated memory. <br>
  !> !! There is no assignment of the lbound and ubound attributes done in the
  !> new leaf as the data values reside directly in the tLeaf structure and not
  !> in a tStreams structure belonging to one of the Leaf's parent tBranch
  !> structures. <br>
  !> !! The sixth element of the pstat attribute arrayis set to 1 to indicate,
  !> that the data resides directly in the leaf.
  Subroutine add_filled_leaf_to_branch_6(t_b, desc, dat_no, values)

    Type(tBranch), Intent(inout)                  :: t_b

    Character(Len=*)   , Intent(in)               :: desc
    Integer(Kind=1)    , PARAMETER                :: dat_ty = 6
    Integer(Kind=pd_ik), Intent(in)               :: dat_no

    character          , Intent(in), Dimension(:) :: values

    !> Pointer to leaf children
    Type(tLeaf), Dimension(:), Pointer   :: leaves

    Integer(kind=pd_ik)                  :: ii
    Integer                              :: alloc_stat

    leaves     => t_b%leaves
    t_b%leaves => Null()

    allocate(t_b%leaves(t_b%no_leaves+1))

    Do ii = 1, t_b%no_leaves
       t_b%leaves(ii)=leaves(ii)
    End Do

    if (t_b%no_leaves > 0) deallocate(leaves)

    t_b%no_leaves = t_b%no_leaves + 1

    t_b%leaves(t_b%no_leaves)%desc   = trim(desc)
    t_b%leaves(t_b%no_leaves)%dat_ty = dat_ty
    t_b%leaves(t_b%no_leaves)%dat_no = dat_no
    t_b%leaves(t_b%no_leaves)%lbound = 0_pd_ik
    t_b%leaves(t_b%no_leaves)%ubound = 0_pd_ik

    t_b%leaves(t_b%no_leaves)%p_int1 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int2 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int4 => NULL()
    t_b%leaves(t_b%no_leaves)%p_int8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_real8 => NULL()
    t_b%leaves(t_b%no_leaves)%p_char => NULL()
    t_b%leaves(t_b%no_leaves)%p_log  => NULL()

    t_b%leaves(t_b%no_leaves)%pstat = 0

    allocate(t_b%leaves(t_b%no_leaves)%p_char(dat_no), stat=alloc_stat)
    call alloc_error(alloc_stat, "t_b%leaves(t_b%no_leaves)%p_char(dat_no)", &
         "add_filled_leaf_to_branch_6") 
    t_b%leaves(t_b%no_leaves)%p_char = values
    
    t_b%leaves(t_b%no_leaves)%pstat = 1

  End Subroutine add_filled_leaf_to_branch_6

  !============================================================================
  !> Subroutine which adds an empty leaf to a tBranch structure
  !>
  !> The subroutine adds an empty  leaf to a tBranch structure at 
  !> the end of its leaves array. Empty means the attributes of the new leaf 
  !> are initialized to the empty string in the case of desc and to 0 in the
  !> case of dat_ty, dat_no, lbound and ubound. <br>
  !> !! All data pointers are initialized to NULL(). <br>
  !> !! The pstat attribute is set to 0.
  Subroutine add_empty_leaf_to_branch(t_b,n)

    Type(tBranch), Intent(inout)             :: t_b
    Integer(kind=pd_ik), intent(in),optional :: n
    
    !> Pointer to leaf children
    Type(tLeaf), Dimension(:), Pointer   :: leaves

    Integer(kind=pd_ik)                  :: ii, loc_n

    If (present(n)) then
       loc_n = n
    Else
       loc_n = 1
    End If

    leaves     => t_b%leaves
    t_b%leaves => Null()

    allocate(t_b%leaves(t_b%no_leaves+loc_n))

    Do ii = 1, t_b%no_leaves
       t_b%leaves(ii)=leaves(ii)
    End Do

    if (associated(leaves)) deallocate(leaves)

    Do ii = t_b%no_leaves+1, t_b%no_leaves+loc_n
       t_b%leaves(ii)%desc   = ''
       t_b%leaves(ii)%dat_ty = 0_1
       t_b%leaves(ii)%dat_no = 0_pd_ik
       t_b%leaves(ii)%lbound = 0_pd_ik
       t_b%leaves(ii)%ubound = 0_pd_ik
       
       t_b%leaves(ii)%p_int1 => NULL()
       t_b%leaves(ii)%p_int2 => NULL()
       t_b%leaves(ii)%p_int4 => NULL()
       t_b%leaves(ii)%p_int8 => NULL()
       t_b%leaves(ii)%p_real8 => NULL()
       t_b%leaves(ii)%p_char => NULL()
       t_b%leaves(ii)%p_log  => NULL()

       t_b%leaves(ii)%pstat = 0
    End Do

    t_b%no_leaves = t_b%no_leaves + loc_n
    
  End Subroutine add_empty_leaf_to_branch

  !============================================================================
  !> Subroutine which adds a branch to a tBranch structure
  !>
  !> The subroutines adds a new element at the end of t_b's branches array.
  !> If n_b is present, it points to the newly created last element of t_b's
  !> branches array.<br>
  !> !! Only the pointer assignment between n_b and the new element is done.
  !> No initialisation is performed. Neither on the new element nor on n_b.
  Subroutine add_branch_to_branch(t_b,n_b)

    Type(tBranch), Intent(inout)                       :: t_b
    Type(tBranch), Intent(  out), Pointer, optional    :: n_b

    !> Pointer to branch children
    Type(tBranch), Dimension(:), Pointer :: branches

    Integer(kind=pd_ik)                  :: ii

    branches     => t_b%branches
    t_b%branches => Null()

    allocate(t_b%branches(t_b%no_branches+1))

    Do ii = 1, t_b%no_branches
       call Assign_branches(t_b%branches(ii),branches(ii))
    End Do

    If (associated(branches)) deallocate(branches)

    t_b%no_branches = t_b%no_branches + 1
    if (present(n_b)) n_b => t_b%branches(t_b%no_branches)

  End Subroutine add_branch_to_branch

  !============================================================================
  !> Subroutine which removes a branch from a tBranch structure
  !>
  !> The Subroutine removes all elements with the given description from a 
  !> branch's branches list. !! No recursive search is performed for desc !!
  !> The branches are recursively deallocted by destroy_tree. No stream size
  !> integrity checks are performed since the parent branch can not be
  !> determined. If correct stream sizes of the parent are needed execution of
  !> set_bounds_in_branch is necessary outside this routine.
  Subroutine delete_branch_from_branch(desc, br, rem_data)

    Character(len=*), Intent(in)               :: desc
    Type(tBranch), intent(InOut)               :: br

    Logical                                    :: success, skip
    Type(tBranch), Dimension(:), Pointer       :: branches
    Integer(kind=pd_ik), Dimension(no_streams) :: rem_data

    Integer(kind=pd_ik)                        :: ii, jj, kk
    
    branches => Null()

    Do ii = 1, br%no_branches

       if (trim(br%branches(ii)%desc) == trim(desc)) then

          skip = .FALSE.
          branches    => br%branches
          br%branches => Null()

          If (br%no_branches > 1) then

             br%no_branches = br%no_branches - 1
             
             allocate(br%branches(br%no_branches))

             kk = 1
             Do jj = 1, br%no_branches
                If ( (trim(branches(kk)%desc) == trim(desc)) .AND. &
                     .not. skip ) then

                   Call destroy_tree(branches(kk), rem_data)

                   kk = kk + 1
                   skip = .TRUE.

                End If
                
                call Assign_branches(br%branches(jj),branches(kk))
                kk = kk + 1
                
             End Do
             
          Else
             
             br%no_branches = br%no_branches - 1
             
          End If
          
          success = .TRUE.
          Exit
          
       End if

    End Do
    
    If (associated(branches)) deallocate(branches)
        
    If (.not. success) then
       Write(pd_umon,PDF_M_A)"Branch with desc        ",trim(desc)
       Write(pd_umon,PDF_M_A)"was not found in branch ",trim(br%desc)
       Write(pd_umon,PDF_M_A)"Nothing deleted"
    End If

  End Subroutine delete_branch_from_branch

  !============================================================================
  !> Subroutine which initializes a new Root Branch
  !> 
  !> \param[in] desc short description of the tree content
  !>
  !> \param[out] The initialized tree
  !>
  !> <b>tBranch is initalised with:</b> \n\n
  !> tBranch::desc   = desc             \n
  !> tBranch::streams::ii_st  = 1_pd_ik          \n
  !> tBranch::streams::dim_st = 0_pd_ik          \n
  !> tBranch::streams::stream_files(1) = Trim(pro_path)//Trim(pro_name)//'.int1.st'  \n
  !> tBranch::streams::stream_files(2) = Trim(pro_path)//Trim(pro_name)//'.int2.st'  \n
  !> tBranch::streams::stream_files(3) = Trim(pro_path)//Trim(pro_name)//'.int4.st'  \n
  !> tBranch::streams::stream_files(4) = Trim(pro_path)//Trim(pro_name)//'.int8.st'  \n
  !> tBranch::streams::stream_files(5) = Trim(pro_path)//Trim(pro_name)//'.real8.st'  \n
  !> tBranch::streams::stream_files(6) = Trim(pro_path)//Trim(pro_name)//'.char.st'  \n
  !> tBranch::streams::stream_files(7) = Trim(pro_path)//Trim(pro_name)//'.log.st'   \n
  !> tBranch::streams::ifopen          = .FALSE.                                     \n
  !> tBranch::streams::units           = -1
  Subroutine raise_tree(desc,tree)

    Character(len=*), Intent(in)  :: desc
    Type(tBranch)   , intent(out) :: tree

    Integer                       :: alloc_stat

    !* Set description ********************************************************
    tree%desc = desc

    !* Init tBranch components ************************************************
    tree%no_branches = 0
    tree%no_leaves   = 0
    tree%branches    => NULL()
    tree%leaves      => NULL()

    !* Streams component ******************************************************
    IF (.not. allocated(tree%streams)) then
       allocate(tree%streams,stat=alloc_stat)
       Call alloc_error(alloc_stat,'tree%streams' , 'raise_tree')
    Else
       Call cons_error("raise_tree", "Streams component is already allocated")
    End IF

    !* stream counters ********************************************************
    tree%streams%ii_st = 1

    !* Init stream dimensions *************************************************
    tree%streams%dim_st = 0_pd_ik

    !* Init stream files logic ************************************************
    tree%streams%stream_files(1) = Trim(pro_path)//Trim(pro_name)//'.int1.st'
    tree%streams%stream_files(2) = Trim(pro_path)//Trim(pro_name)//'.int2.st'
    tree%streams%stream_files(3) = Trim(pro_path)//Trim(pro_name)//'.int4.st'
    tree%streams%stream_files(4) = Trim(pro_path)//Trim(pro_name)//'.int8.st'
    tree%streams%stream_files(5) = Trim(pro_path)//Trim(pro_name)//'.real8.st'
    tree%streams%stream_files(6) = Trim(pro_path)//Trim(pro_name)//'.char.st'
    tree%streams%stream_files(7) = Trim(pro_path)//Trim(pro_name)//'.log.st'

    tree%streams%ifopen = .FALSE.
    tree%streams%units  = -1

    tree%streams%int1_st  => null()
    tree%streams%int2_st  => null()
    tree%streams%int4_st  => null()
    tree%streams%int8_st  => null()
    tree%streams%real8_st => null()
    tree%streams%char_st  => null()
    tree%streams%log_st   => null()

  End Subroutine raise_tree

  !============================================================================
  !> Subroutine for tree deallocation
  !>
  !> The subroutine recursively deallocates a complete puredat data-location
  !> tree set up by tBranch and tLeaf elements by recursively cycling 
  !> the complete tree structure.<br>
  !> !! This routine does not take care about the parent branch and stream
  !> structures !!
  Recursive Subroutine destroy_tree(branch, no_data)

    Type(tBranch), Intent(InOut)                              :: branch
    Integer(kind=pd_ik), Dimension(no_streams), Intent(InOut) :: no_data

    Integer :: ii, alloc_stat

    If (Allocated(Branch%streams)) then
       IF (associated(Branch%streams%int1_st)) then
          deallocate(Branch%streams%int1_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%int1_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%int2_st)) then
          deallocate(Branch%streams%int2_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%int2_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%int4_st)) then
          deallocate(Branch%streams%int4_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%int4_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%int8_st)) then
          deallocate(Branch%streams%int8_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%int8_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%real8_st)) then
          deallocate(Branch%streams%real8_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%real8_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%char_st)) then
          deallocate(Branch%streams%char_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%char_st', 'destroy_tree')
       End IF
       IF (associated(Branch%streams%log_st)) then
          deallocate(Branch%streams%log_st,stat=alloc_stat)
          Call dealloc_error(alloc_stat, &
               'Branch%streams%log_st', 'destroy_tree')
       End IF

       deallocate(Branch%streams,stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'Branch%streams', 'destroy_tree')

    End If

    If (Associated(Branch%leaves)) Then

       Do ii = 1, Branch%no_leaves

          If (Associated(Branch%leaves(ii)%p_int1) .AND. &
               (Branch%leaves(ii)%pstat == 1) ) then
             no_data(1) = no_data(1) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_int1,stat=alloc_stat)
          End If
          If (Associated(Branch%leaves(ii)%p_int2) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(2) = no_data(2) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_int2,stat=alloc_stat)
          ENd If
          If (Associated(Branch%leaves(ii)%p_int4) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(3) = no_data(3) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_int4,stat=alloc_stat)
          End If
          If (Associated(Branch%leaves(ii)%p_int8) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(4) = no_data(4) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_int8,stat=alloc_stat)
          End If
          If (Associated(Branch%leaves(ii)%p_real8) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(5) = no_data(5) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_real8,stat=alloc_stat)
          End If
          If (Associated(Branch%leaves(ii)%p_char) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(6) = no_data(6) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_char,stat=alloc_stat)
          End If
          If (Associated(Branch%leaves(ii)%p_log) .AND. &
               (Branch%leaves(ii)%pstat == 1) )  then
             no_data(7) = no_data(7) + Branch%leaves(ii)%dat_no
             deallocate(Branch%leaves(ii)%p_log,stat=alloc_stat)
          End If
       End Do

       Deallocate(branch%leaves,stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'Branch%leaves', 'destroy_tree')

    End If

    Do ii = 1, Branch%no_branches
       Call destroy_tree(Branch%branches(ii), no_data)
    End Do

    If (Associated(Branch%branches)) Then
       Deallocate(Branch%branches,stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'Branch%branches', 'destroy_tree')
    End If

    branch%desc        = ""
    branch%no_leaves   = 0
    branch%no_branches = 0

  End Subroutine destroy_tree

  !============================================================================
  !> Subroutine which destroies a streams component of a tBranch structure
  !>
  !> Subroutine for deallocating a streams component of a tBranch structure 
  !>
  !> \param[in] streams The streams component to be deallocated
  Subroutine destroy_streams(streams)

    Type(tStreams), Allocatable :: streams
    Integer                     :: alloc_stat

    If((streams%dim_st(1) > 0) .AND. associated(streams%int1_st)) then
       deallocate(streams%int1_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%int1_st', 'destroy_streams')
       streams%int1_st=>null()
    End If
    If((streams%dim_st(2) > 0) .AND. associated(streams%int2_st)) then
       deallocate(streams%int2_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%int2_st', 'destroy_streams')
       streams%int2_st=>null()
    End If
    If((streams%dim_st(3) > 0) .AND. associated(streams%int4_st)) then
       deallocate(streams%int4_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%int4_st', 'destroy_streams')
       streams%int4_st=>null()
    End If
    If((streams%dim_st(4) > 0) .AND. associated(streams%int8_st)) then 
       deallocate(streams%int8_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%int8_st', 'destroy_streams')
       streams%int8_st=>null()
    End If
    If((streams%dim_st(5) > 0) .AND. associated(streams%real8_st)) then
       deallocate(streams%real8_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%real8_st', 'destroy_streams')
       streams%real8_st=>null()
    End If
    If((streams%dim_st(6) > 0) .AND. associated(streams%char_st)) then
       deallocate(streams%char_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%char_st', 'destroy_streams')
       streams%char_st=>null()
    End If
    If((streams%dim_st(7) > 0) .AND. associated(streams%log_st)) then
       deallocate(streams%log_st, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams%log_st', 'destroy_streams')
       streams%log_st=>null()
    End If

    if(allocated(streams)) then
       deallocate(streams, stat=alloc_stat)
       Call dealloc_error(alloc_stat, 'streams', 'destroy_streams')
    End if

  End Subroutine destroy_streams

  !============================================================================
  !> Subroutine which raises a new branch form scratch
  !>
  !> Subroutine for allocating a new branch with a specified number of leaves 
  !> and child branches from which all content is initialised to zero
  !>
  !> \param[in] desc short description of the content
  !> \param[in] no_branches Number of child branches to be allocated
  !> \param[in] no_leaves Number of child leaves to be allocated
  !>
  !> \param[inout] The raised branch
  Subroutine raise_branch(desc, no_branches, no_leaves, branch)

    Character(Len=*)   , Intent(In)     :: desc
    Integer(Kind=pd_ik), Intent(In)     :: no_branches
    Integer(Kind=pd_ik), Intent(In)     :: no_leaves
    Type(tBranch)      , Intent(InOut)  :: branch

    Integer :: alloc_stat, ii

    branch%desc = desc

    branch%no_branches = no_branches
    branch%no_leaves   = no_leaves

    Allocate(branch%branches(no_branches), Stat=alloc_stat)
    Call alloc_error(alloc_stat, 'branch%branches' ,&
         'raise_branches_from_scratch', no_branches)

    Do ii = 1, no_branches 
       branch%branches(ii)%desc = ''
       branch%branches(ii)%no_branches = 0 
       branch%branches(ii)%no_leaves   = 0 
       Nullify(branch%branches(ii)%branches)
       Nullify(branch%branches(ii)%leaves)
    End Do

    Allocate(branch%leaves(no_leaves), Stat=alloc_stat)
    Call alloc_error(alloc_stat,'branch%leaves' ,'raise_leaves', no_leaves)

    Do ii = 1, no_leaves
       branch%leaves(ii)%desc = ''
       branch%leaves(ii)%dat_no = 0
       branch%leaves(ii)%dat_ty = 0
       branch%leaves(ii)%lbound = 0
       branch%leaves(ii)%ubound = 0
       Nullify(branch%leaves(ii)%p_int1)
       Nullify(branch%leaves(ii)%p_int2)
       Nullify(branch%leaves(ii)%p_int4)
       Nullify(branch%leaves(ii)%p_int8)
       Nullify(branch%leaves(ii)%p_real8)
       Nullify(branch%leaves(ii)%p_char)
       Nullify(branch%leaves(ii)%p_log)
       branch%leaves(ii)%pstat = 0
    End Do

  End Subroutine raise_branch

  !============================================================================
  !> Subroutine for raising leaves to a tree
  !>
  !> This subroutine initalises the leaves of a branch
  !>
  !> \param[in] no_leaves Number of child leaves to be initialised
  !> \param[in] desc   Array of size (no_leaves) with leaf descriptions
  !> \param[in] dat_ty Array of size (no_leaves) with leaf data type numbers
  !> \param[in] dat_no Array of size (no_leaves) with te number of leaf data 
  !>                   elements
  !>
  !> \param[in] lb ! OPTIONAL ! Array of size (no_leaves) with lower stream 
  !>                            bounds
  !> \param[in] ub ! OPTIONAL ! Array of size (no_leaves) with upper stream 
  !>                            bounds
  !>
  !> \param branch Branch of Type(tBranch) whose leaves should be initialised
  !>
  !> If lb and ub are not given, they are initalised to zero. The stream bounds 
  !> of the data chunks can be set afterwards for a complete Type(tBranch) by 
  !> calling "set_bounds_in_branch"
  Subroutine raise_leaves(no_leaves, desc, dat_ty, dat_no, lb, ub, branch)

    Integer         , Intent(in)                       :: no_leaves

    Character(Len=*), Intent(in), Dimension(no_leaves) :: desc

    Integer(Kind=1) , Intent(in), Dimension(no_leaves) :: dat_ty 

    Integer(Kind=pd_ik), Intent(in), Dimension(no_leaves) :: dat_no

    Integer(Kind=pd_ik), Intent(in), Dimension(no_leaves), optional :: lb
    Integer(Kind=pd_ik), Intent(in), Dimension(no_leaves), optional :: ub

    Type(tBranch), Intent(InOut) :: branch

    Integer :: ii

    Do ii = 1, no_leaves

       branch%leaves(ii)%desc = desc(ii)

       branch%leaves(ii)%dat_ty = dat_ty(ii)
       branch%leaves(ii)%dat_no = dat_no(ii)

    End Do

    if (present(lb) .and. present(ub)) then
       Do ii = 1, no_leaves    
          branch%leaves(ii)%lbound = lb(ii)
          branch%leaves(ii)%ubound = ub(ii)
       End Do
    Else

       Do ii = 1, no_leaves    
          branch%leaves(ii)%lbound = 0
          branch%leaves(ii)%ubound = 0
       End Do

    End if

  End Subroutine raise_leaves

  !============================================================================
  !> Subroutine: Recursive subroutine to connect stream chunk pointers
  !>
  !> This recursive subroutine connects all stream chunk pointers specified
  !> by the leaves of type tLeaf in a tBranch structure to their stream targets
  Recursive Subroutine connect_pointers(streams, branch)

    Type(tStreams)   , Intent(inout)         :: streams
    Type(tBranch)    , Intent(inout)         :: branch

    Integer :: ii

    Do ii =1, branch%no_branches
       Call connect_pointers(streams, branch%branches(ii))
    End Do

    Do ii = 1, branch%no_leaves

       If (branch%leaves(ii)%pstat == 0) then
       
          Select Case(branch%leaves(ii)%dat_ty)
             
          Case(1)
             branch%leaves(ii)%p_int1 => &
                  streams%int1_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(2)
             branch%leaves(ii)%p_int2 => &
                  streams%int2_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(3)
             branch%leaves(ii)%p_int4 => &
                  streams%int4_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(4)
             branch%leaves(ii)%p_int8 => &
                  streams%int8_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(5)
             branch%leaves(ii)%p_real8 => &
                  streams%real8_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(6)
             branch%leaves(ii)%p_char => &
                  streams%char_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
          Case(7)
             branch%leaves(ii)%p_log => &
                  streams%log_st(branch%leaves(ii)%lbound:branch%leaves(ii)%ubound)
             
          Case default
             Write(pd_umon,*)'!!!!! Bad data type in leaf descriptor !!!!!!'
             Write(pd_umon,*)'!! See last leaf descriptor in header file !!'
             Write(pd_umon,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             Write(pd_umon,*)'!! ========= PROGRAM TERMINATED ========== !!'
             Stop
          End Select

       Else if (branch%leaves(ii)%pstat == 1) then

          Write(pd_umon,PDF_W_A)'In connect_pointers, leaf '
          Write(pd_umon,PDF_W_A)trim(branch%leaves(ii)%desc)
          Write(pd_umon,PDF_W_A)'has pstat = 1. No pointer reasignment is done !'
          
       End If
          
    End Do

  End Subroutine connect_pointers

  !============================================================================
  !> Subroutine which includes a branch into another branch
  !>
  !> The include of the data is done based on files. The include of the branch
  !> structure is done based on pointers. 
  !> !! If the source branch is destroyed afterwards the included branch is  !!
  !> !! also lost. No copy of the childs is performed.                       !!
  !>
  !> s_b       : Has to be a root branch which holds all data of its childs
  !>
  !> t_b       : Has to be a root branch without t_streams present. With 
  !>             t_streams present it has to be a child of the root branch 
  !>             holding t_streams
  !> t_streams : If present t_steams has to be the stream component of t_b's
  !>             parent. Otherwise the storage place gets seperated from its
  !>             tree. Integrity checks of structure and storage place are 
  !>             not performed !!!
  !>
  !> \todo Add exception if stream files are not opened
  Subroutine include_branch_into_branch(s_b, t_b, s_streams, t_streams, &
                                        clean_target_files, blind)

    Type(tBranch) , Intent(inout)                        :: s_b
    Type(tBranch) , Intent(inout)                        :: t_b
    Type(tStreams), Intent(inout), Allocatable, Optional :: s_streams 
    Type(tStreams), Intent(inout), Allocatable, Optional :: t_streams 
    Type(tStreams)               , Allocatable           :: src_streams
    Type(tStreams)                                       :: trg_streams

    Logical, Intent(in), optional                        :: clean_target_files
    Logical, intent(in), Optional                        :: blind

    Logical                                 :: s_root = .True.
    Logical                                 :: t_root = .True.
    Logical                                 :: ctf    = .False.

    Logical                                 :: loc_blind

    !> Pointer to branch children
    Type(tBranch), Dimension(:), Pointer :: branches

    Integer(kind=pd_ik)                  :: ii

    !** Handle optional parameters ********************************************

    if (present(t_streams)) then

       if (.not. allocated(t_streams)) then
          Write(pd_umon,*)"In Subroutine include_branch_into_branch"
          Write(pd_umon,*)"An inconsistency was detected "
          Write(pd_umon,*)"t_streams is present without being allocated"
          stop
       End If

       trg_streams = t_streams
       t_root   = .False.
    ELse

       if (.not. allocated(t_b%streams)) then
          Write(pd_umon,*)"In Subroutine include_branch_into_branch"
          Write(pd_umon,*)"An inconsistency was detected "
          Write(pd_umon,*)"t_b has no allocated stream component with t_streams"
          write(pd_umon,*)"not being present"
          stop
       End If

       trg_streams = t_b%streams
    End if

    allocate(src_streams)
    if (present(s_streams)) then

       if (.not. allocated(s_streams)) then
          Write(pd_umon,*)"In Subroutine include_branch_into_branch"
          Write(pd_umon,*)"An inconsistency was detected "
          Write(pd_umon,*)"s_streams is present without being allocated"
          stop
       End If

       src_streams = s_streams
       s_root   = .False.
    ELse

       if (.not. allocated(s_b%streams)) then
          Write(pd_umon,*)"In Subroutine include_branch_into_branch"
          Write(pd_umon,*)"An inconsistency was detected "
          Write(pd_umon,*)"s_b has no allocated stream component with s_streams"
          write(pd_umon,*)"not being present"
          stop
       End If

       src_streams = s_b%streams
    End if

    if (present(clean_target_files)) ctf = clean_target_files

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    !** Include branch structure **********************************************

    branches     => t_b%branches
    t_b%branches => Null()

    allocate(t_b%branches(t_b%no_branches+1))

    Do ii = 1, t_b%no_branches
       call Assign_branches(t_b%branches(ii),branches(ii))
    End Do

    if (t_b%no_branches>0) deallocate(branches)

    t_b%no_branches = t_b%no_branches + 1

    t_b%branches(t_b%no_branches)%desc        =  s_b%desc
    t_b%branches(t_b%no_branches)%no_branches =  s_b%no_branches
    t_b%branches(t_b%no_branches)%no_leaves   =  s_b%no_leaves

    If (s_b%no_branches > 0) then
       t_b%branches(t_b%no_branches)%branches    => s_b%branches
    Else
       nullify(t_b%branches(t_b%no_branches)%branches)
    End If

    If (s_b%no_leaves > 0) then
       t_b%branches(t_b%no_branches)%leaves      => s_b%leaves
    Else
       nullify(t_b%branches(t_b%no_branches)%leaves)
    End If

    !** Include data **********************************************************
    If (.not.loc_blind) then
       call read_streams(src_streams)
       Call open_stream_files(trg_streams,"write","old","append")
    End If

    !** Usecase 1 : ***********************************************************
    !** s_b is a root branch. This means we can copy its streams to trg_streams
    if (s_root) then

       If (.not.loc_blind) then
          If (src_streams%dim_st(1) > 0) write(trg_streams%units(1))src_streams%int1_st
          If (src_streams%dim_st(2) > 0) write(trg_streams%units(2))src_streams%int2_st
          If (src_streams%dim_st(3) > 0) write(trg_streams%units(3))src_streams%int4_st
          If (src_streams%dim_st(4) > 0) write(trg_streams%units(4))src_streams%int8_st
          If (src_streams%dim_st(5) > 0) write(trg_streams%units(5))src_streams%real8_st
          If (src_streams%dim_st(6) > 0) write(trg_streams%units(6))src_streams%char_st
       End If

       call shift_bounds_in_branch(t_b%branches(t_b%no_branches), trg_streams, src_streams)

       !** Usecase 2 : ***********************************************************
       !** s_b is not a root branch. This means we have to cycle it recursively 
       !** and copy its data leaf by leaf
    Else

       call connect_pointers(src_streams,s_b)
       call store_branch(s_b,trg_streams,loc_blind)

    End if

    if (ctf .AND. (.NOT.loc_blind)) call close_stream_files(trg_streams,ctf)
    call destroy_streams(src_streams)

    if (present(t_streams)) then      
       t_streams = trg_streams
    ELse
       t_b%streams = trg_streams
    End if

  End Subroutine include_branch_into_branch

  !****************************************************************************
  !> Subroutine for storing all the data of a branch in memory to stream files
  Recursive Subroutine store_branch(branch,streams,blind)

    Type(tBranch) , Intent(InOut) :: branch
    Type(tStreams), Intent(InOut) :: streams

    Integer                        :: ii
    Logical, intent(in), Optional  :: blind
    Logical                        :: loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    Do ii = 1, branch%no_branches
       Call store_branch(branch%branches(ii),streams,loc_blind)
    End Do

    Do ii = 1, branch%no_leaves

       If (branch%leaves(ii)%dat_ty == 1) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_int1,loc_blind)
       End If
       If (branch%leaves(ii)%dat_ty == 2) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_int2,loc_blind)
       End If
       If (branch%leaves(ii)%dat_ty == 3) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_int4,loc_blind)
       End If
       If (branch%leaves(ii)%dat_ty == 4) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_int8,loc_blind)
       End If
       If (branch%leaves(ii)%dat_ty == 5) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_real8,loc_blind)
       End If
       If (branch%leaves(ii)%dat_ty == 6) then
          call pd_store(streams,branch,trim(branch%leaves(ii)%desc),&
               branch%leaves(ii)%p_char,loc_blind)
       End If

    End Do

  End Subroutine store_branch

  !****************************************************************************
  !> Subroutine for setting data chunk bounds
  !>
  !> This subroutine sets the lbound and ubound members of all tLeaf members 
  !> in a tBranch structure.
  recursive subroutine reset_bounds_in_branch(branch, t_streams)

    Type(tBranch) , Intent(InOut)  :: branch
    Type(tStreams), Intent(InOut)  :: t_streams
   
    Integer(kind=pd_ik)            :: ii
    
    Do ii = 1, branch%no_branches 

       Call reset_bounds_in_branch(branch%branches(ii), t_streams)

    End Do

    Do ii = 1, branch%no_leaves 

       branch%leaves(ii)%lbound = t_streams%ii_st(branch%leaves(ii)%dat_ty)
       branch%leaves(ii)%ubound = branch%leaves(ii)%lbound + branch%leaves(ii)%dat_no - 1
       t_streams%ii_st(branch%leaves(ii)%dat_ty) = branch%leaves(ii)%ubound + 1

    End Do

    t_streams%dim_st  = t_streams%ii_st - 1
    
  End subroutine reset_bounds_in_branch
  
  !****************************************************************************
  !> Subroutine for setting data chunk bounds
  !>
  !> This subroutine sets the lbound and ubound members of tLeaf members in a
  !> tBranch structure. <br>
  !> If leafs are encounterd with lbound and ubound not equal zero the bounds
  !> are not reset since it is assumed that these leaves are already valid
  !> parts of the target streams t_streams.
  recursive subroutine set_bounds_in_branch(branch, t_streams)

    Type(tBranch) , Intent(InOut)  :: branch
    Type(tStreams), Intent(InOut)  :: t_streams
   
    Integer(kind=pd_ik)            :: ii
    
    Do ii = 1, branch%no_branches 

       Call set_bounds_in_branch(branch%branches(ii), t_streams)

    End Do

    Do ii = 1, branch%no_leaves 
       
       If (  (branch%leaves(ii)%lbound == 0) .AND. &
             (branch%leaves(ii)%ubound == 0)        ) then

          branch%leaves(ii)%lbound = t_streams%ii_st(branch%leaves(ii)%dat_ty)
          branch%leaves(ii)%ubound = branch%leaves(ii)%lbound + branch%leaves(ii)%dat_no - 1
          t_streams%ii_st(branch%leaves(ii)%dat_ty) = branch%leaves(ii)%ubound + 1

       End If
                         
    End Do

    t_streams%dim_st  = t_streams%ii_st - 1
    
  End subroutine set_bounds_in_branch

  !****************************************************************************
  !> Subroutine for shifting data chunk bounds
  !>
  !> This subroutine shifts the lbound and ubound members of all tLeaf
  !> members in a tBranch structure
  subroutine shift_bounds_in_branch(branch, t_streams, s_streams)

    Type(tBranch) , Intent(InOut)  :: branch
    Type(tStreams), Intent(InOut)  :: t_streams
    Type(tStreams), Intent(In   )  :: s_streams

    Integer(kind=pd_ik), dimension(no_streams) :: offset

    offset = t_streams%dim_st

    Call shift_bounds_rec(branch, offset)
    t_streams%dim_st = t_streams%dim_st + s_streams%dim_st
    t_streams%ii_st  = t_streams%dim_st + 1

  End subroutine shift_bounds_in_branch

  !****************************************************************************
  !> Subroutine for setting data chunk bounds
  !>
  !> This subroutine sets the lbound and ubound members of all tLeaf
  !> members in a tBranch structure
  recursive subroutine shift_bounds_rec(branch, offset)

    Type(tBranch)                             , Intent(InOut) :: branch
    Integer(kind=pd_ik), dimension(no_streams), intent(in)    :: offset

    Integer                                                   :: ii

    Do ii = 1, branch%no_branches 

       Call shift_bounds_rec(branch%branches(ii), offset)

    End Do

    Do ii = 1, branch%no_leaves 

       if ( (branch%leaves(ii)%dat_ty > 0) .AND. &
            (branch%leaves(ii)%dat_ty <= no_streams) ) then
          branch%leaves(ii)%lbound = &
               offset(branch%leaves(ii)%dat_ty)+branch%leaves(ii)%lbound
          branch%leaves(ii)%ubound = &
               offset(branch%leaves(ii)%dat_ty)+branch%leaves(ii)%ubound
       end if
       
    End Do

  End subroutine shift_bounds_rec

  !****************************************************************************
  !> Subroutine for counting the branches and leaves in a tBranch structure
  Recursive Subroutine count_elems(branch, no_branches, no_leaves)
    
    Type(tBranch), Intent(In)          :: branch
    Integer(Kind=pd_ik), Intent(InOut) :: no_branches, no_leaves 

    Integer(Kind=pd_ik)                :: ii
    
    no_branches = no_branches + branch%no_branches
    no_leaves   = no_leaves   + branch%no_leaves  

    Do ii = 1, branch%no_branches
       Call count_elems(branch%branches(ii), no_branches, no_leaves)
    End Do

  End Subroutine count_elems

  !****************************************************************************
  !> Subroutine for homogenization of a tBranch structure to streams in memory
  !>
  !> This subroutine copys the leaf local data of a tBranch structure to
  !> streams in memory (e.g. to prepare I/O).
  !> This is done for all leafs with pstat = 1 and pstat = 0.
  !> if pstat = 0, which means that the data are not actually part of the leaf
  !> but reside somewhere else in a tStream component, a warning
  !> is issued that the pointer to the original data may be lost since the data
  !> are copied and the pointer is relocated to the location of the copied data
  !> in the streams containing the homogenized data of the branch.<br>
  !> Remark: The streams of the streams input parameter can not already be 
  !> allocated on input !
  Subroutine homogenize_branch(branch, streams)

    Type(tBranch) , Intent(InOut) :: branch
    Type(tStreams), Intent(InOut) :: streams

    Integer(Kind=pd_ik), Dimension(no_streams) :: data_size

    Integer(Kind=pd_ik)           :: ii
    
    If ( associated(streams%int1_st ) .OR. &
         associated(streams%int2_st ) .OR. &
         associated(streams%int4_st ) .OR. &
         associated(streams%int8_st ) .OR. &
         associated(streams%real8_st) .OR. &
         associated(streams%char_st ) .OR. &
         associated(streams%log_st  )       ) then

       Write(pd_umon,PDF_E_A)"In homogenize_branch : "
       Write(pd_umon,PDF_E_A)"One or more of the stream pointers are allready allocated"
       Write(pd_umon,PDF_E_STOP)
       STOP

    End If

    Call get_stream_size(branch,data_size)

    Allocate(streams%int1_st (data_size(1)))
    Allocate(streams%int2_st (data_size(2)))
    Allocate(streams%int4_st (data_size(3)))
    Allocate(streams%int8_st (data_size(4)))
    Allocate(streams%real8_st(data_size(5)))
    Allocate(streams%char_st (data_size(6)))
    Allocate(streams%log_st  (data_size(7)))

    Do ii = 1, branch%no_branches

       Call copy_leaves_to_streams(branch%branches(ii),streams)

    End Do

    Do ii = 1, branch%no_leaves 

       If ( (branch%leaves(ii)%dat_no >  0) .And. &
            (branch%leaves(ii)%pstat  >= 0) ) then

          Select Case (branch%leaves(ii)%dat_ty)

          Case(1)

             streams%int1_st(&
                  streams%ii_st(1) : streams%ii_st(1) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int1

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int1)
             Else if (associated(branch%leaves(ii)%p_int1) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int1 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int1 => &
                  streams%int1_st(streams%ii_st(1) : streams%ii_st(1) + branch%leaves(ii)%dat_no - 1)

          Case(2)

             streams%int2_st(&
                  streams%ii_st(2) : streams%ii_st(2) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int2

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int2)
             Else if (associated(branch%leaves(ii)%p_int2) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int2 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int2 => &
                  streams%int2_st(streams%ii_st(2) : streams%ii_st(2) + branch%leaves(ii)%dat_no - 1)

          Case(3)

             streams%int4_st(&
                  streams%ii_st(3) : streams%ii_st(3) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int4

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int4)
             Else if (associated(branch%leaves(ii)%p_int4) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int4 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int4 => &
                  streams%int4_st(streams%ii_st(3) : streams%ii_st(3) + branch%leaves(ii)%dat_no - 1)

          Case(4)

             streams%int8_st(&
                  streams%ii_st(4) : streams%ii_st(4) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int8

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int8)
             Else if (associated(branch%leaves(ii)%p_int8) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int8 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int8 => &
                  streams%int8_st(streams%ii_st(4) : streams%ii_st(4) + branch%leaves(ii)%dat_no - 1)

          Case(5)

             streams%real8_st(&
                  streams%ii_st(5) : streams%ii_st(5) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_real8

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_real8)
             Else if (associated(branch%leaves(ii)%p_real8) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_real8 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_real8 => &
                  streams%real8_st(streams%ii_st(5) : streams%ii_st(5) + branch%leaves(ii)%dat_no - 1)

          Case(6)

             streams%char_st(&
                  streams%ii_st(6) : streams%ii_st(6) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_char

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_char)
             Else if (associated(branch%leaves(ii)%p_char) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_char is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_char => &
                  streams%char_st(streams%ii_st(6) : streams%ii_st(6) + branch%leaves(ii)%dat_no - 1)

          Case(7)

             streams%log_st(&
                  streams%ii_st(7) : streams%ii_st(7) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_log

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_log)
             Else if (associated(branch%leaves(ii)%p_log) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_log is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_log => &
                  streams%log_st(streams%ii_st(7) : streams%ii_st(7) + branch%leaves(ii)%dat_no - 1)

          End Select

          branch%leaves(ii)%pstat = 0

          streams%ii_st(branch%leaves(ii)%dat_ty) = &
               streams%ii_st(branch%leaves(ii)%dat_ty) + branch%leaves(ii)%dat_no
          streams%dim_st = streams%ii_st - 1
          
       End If
       
    End Do
    
  End Subroutine homogenize_branch
  
  !****************************************************************************
  !> Subroutine for homogenization of a tBranch structure to streams in memory
  !>
  !> This subroutine copys the leaf local data of a tBranch structure to
  !> streams in memory (e.g. to prepare I/O). This routine is meant to be
  !> called by homogenize_branch only. 
  Recursive Subroutine copy_leaves_to_streams(branch, streams)

    Type(tBranch) , Intent(InOut) :: branch
    Type(tStreams), Intent(InOut) :: streams

    Integer(Kind=pd_ik)           :: ii

    Do ii = 1, branch%no_branches

       Call copy_leaves_to_streams(branch%branches(ii),streams)

    End Do

    Do ii = 1, branch%no_leaves 

       If ( (branch%leaves(ii)%dat_no > 0 ) .AND. &
            (branch%leaves(ii)%pstat  >= 0) ) then

          Select Case (branch%leaves(ii)%dat_ty)

          Case(1)

             streams%int1_st(&
                  streams%ii_st(1) : streams%ii_st(1) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int1

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int1)
             Else if (associated(branch%leaves(ii)%p_int1) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int1 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int1 => &
                  streams%int1_st(streams%ii_st(1) : streams%ii_st(1) + branch%leaves(ii)%dat_no - 1)

          Case(2)

             streams%int2_st(&
                  streams%ii_st(2) : streams%ii_st(2) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int2

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int2)
             Else if (associated(branch%leaves(ii)%p_int2) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int2 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int2 => &
                  streams%int2_st(streams%ii_st(2) : streams%ii_st(2) + branch%leaves(ii)%dat_no - 1)

          Case(3)

             streams%int4_st(&
                  streams%ii_st(3) : streams%ii_st(3) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int4

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int4)
             Else if (associated(branch%leaves(ii)%p_int4) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int4 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int4 => &
                  streams%int4_st(streams%ii_st(3) : streams%ii_st(3) + branch%leaves(ii)%dat_no - 1)

          Case(4)

             streams%int8_st(&
                  streams%ii_st(4) : streams%ii_st(4) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_int8

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_int8)
             Else if (associated(branch%leaves(ii)%p_int8) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_int8 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_int8 => &
                  streams%int8_st(streams%ii_st(4) : streams%ii_st(4) + branch%leaves(ii)%dat_no - 1)

          Case(5)

             streams%real8_st(&
                  streams%ii_st(5) : streams%ii_st(5) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_real8

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_real8)
             Else if (associated(branch%leaves(ii)%p_real8) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_real8 is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_real8 => &
                  streams%real8_st(streams%ii_st(5) : streams%ii_st(5) + branch%leaves(ii)%dat_no - 1)

          Case(6)

             streams%char_st(&
                  streams%ii_st(6) : streams%ii_st(6) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_char

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_char)
             Else if (associated(branch%leaves(ii)%p_char) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_char is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_char => &
                  streams%char_st(streams%ii_st(6) : streams%ii_st(6) + branch%leaves(ii)%dat_no - 1)

          Case(7)

             streams%log_st(&
                  streams%ii_st(7) : streams%ii_st(7) + branch%leaves(ii)%dat_no - 1) &
                  = branch%leaves(ii)%p_log

             If ( branch%leaves(ii)%pstat == 1 ) then
                deallocate(branch%leaves(ii)%p_log)
             Else if (associated(branch%leaves(ii)%p_log) .AND.  &
                  (branch%leaves(ii)%pstat == 0) ) then
                Write(pd_umon,PDF_W_A)"In 'copy_leaves_to_streams' possible memory leak detected !!!"
                Write(pd_umon,PDF_W_A)"p_log is associated with pstat == 0 !!!!!!!!!!!!!!!!!!!!!!!!"
             End If

             branch%leaves(ii)%p_log => &
                  streams%log_st(streams%ii_st(7) : streams%ii_st(7) + branch%leaves(ii)%dat_no - 1)

          End Select

          branch%leaves(ii)%pstat = 0

          streams%ii_st(branch%leaves(ii)%dat_ty) = &
               streams%ii_st(branch%leaves(ii)%dat_ty) + branch%leaves(ii)%dat_no
          streams%dim_st = streams%ii_st - 1
          
       End If
       
    End Do
    
  End Subroutine copy_leaves_to_streams
  !> @} 
  !# End of memeber group "Puredat structure routines" ########################

  !############################################################################
  !> \name Puredat stream file routines
  !> @{
  !> Module procedures which operate on the stream files of a tree

  !============================================================================
  !> Subroutine to reset stream file names
  Subroutine set_stream_fienames(streams)

    Type(tstreams) :: streams
    
    !* Init stream files logic ************************************************
    streams%stream_files(1) = Trim(pro_path)//Trim(pro_name)//'.int1.st'
    streams%stream_files(2) = Trim(pro_path)//Trim(pro_name)//'.int2.st'
    streams%stream_files(3) = Trim(pro_path)//Trim(pro_name)//'.int4.st'
    streams%stream_files(4) = Trim(pro_path)//Trim(pro_name)//'.int8.st'
    streams%stream_files(5) = Trim(pro_path)//Trim(pro_name)//'.real8.st'
    streams%stream_files(6) = Trim(pro_path)//Trim(pro_name)//'.char.st'
    streams%stream_files(7) = Trim(pro_path)//Trim(pro_name)//'.log.st'

  End Subroutine set_stream_fienames
  
  !============================================================================
  !> Subroutine which connects stream files
  !> 
  !> The subroutine opens all stream files that belong to  puredat tree 
  !> structure it takes the followig parametres IN SMALL LETERS for action and 
  !> status :
  !> ACTION : "read", "write"
  !> STATUS : "old", "new", "replace" 
  !> They are passed tro the corresponding prameters of the Fortran open 
  !> statement for
  Subroutine open_stream_files_from_tree(tree, action, status, position)

    Type(tBranch)   , Intent(InOut) :: tree
    Character(len=*), intent(In)    :: action, status

    Character(len=*), intent(In),optional :: position
    Character(len=pd_mcl)                 :: lpos

    Character(len=10)               :: faction

    logical :: fexist
    Integer :: ii, funit

    If (present(position)) then
       lpos = position
    Else
       lpos = "REWIND"
    End If

    If (.NOT. allocated(tree%streams)) then
       call alloc_error(-1, "tree%streams", "open_stream_files_from_tree")
    End If

    Do ii = 1, no_streams

       Inquire(file=trim(tree%streams%stream_files(ii)), exist=fexist, &
            number=funit, action=faction) 

       !** File exists ********************************************************
       If (fexist) Then

          !** File is connected to a unit *************************************
          If (funit > -1) then

             !** Check tree%units integrity ***********************************
             if (tree%streams%units(ii) /= funit) then

                Write(pd_umon,*)"In Subroutine open_stream_files_from_tree'"
                Write(pd_umon,*)"An inconsistency was detected in tree%units"
                Write(pd_umon,*)"Expected unit no. is     :",tree%streams%units(ii)
                write(pd_umon,*)"Inqure returned unit no. :",funit
                Write(pd_umon,*)"PROGRAM STOPPED !!"
                stop

             End if

             !** Check whether the file is opened for the desired action ******
             if (action /= trim(faction)) then

                close(tree%streams%units(ii))
                Open(unit=tree%streams%units(ii), &
                     file=tree%streams%stream_files(ii), status=status, &
                     action=action, access='stream', form='unformatted',&
                     position=trim(lpos))

             End if

             !** File is not connected to a unit *********************************
          Else

             tree%streams%units(ii) = pd_give_new_unit()
             Open(unit=tree%streams%units(ii), &
                  file=tree%streams%stream_files(ii), status=status, &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             tree%streams%ifopen(ii) = .TRUE.

          End If

          !** File does not exist ************************************************
       Else 

          If ((status == 'new') .or. (status == 'replace')) then

             tree%streams%units(ii) = pd_give_new_unit()
             Open(unit=tree%streams%units(ii), &
                  file=tree%streams%stream_files(ii), status=status, &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             tree%streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine open_stream_files_from_tree' file-open action ", action
                Write(pd_umon,*)"is used on the empty file :"
                write(pd_umon,*)trim(tree%streams%stream_files(ii))
             End If

          Else If ( (status =='old') .AND. ( trim(lpos) == 'append' ) ) then

             tree%streams%units(ii) = pd_give_new_unit()
             Open(unit=tree%streams%units(ii), &
                  file=tree%streams%stream_files(ii) , status='new', &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             tree%streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine 'open_stream_files_from_tree'"
                Write(pd_umon,*)"file-open status  ", status," ws used on a    "
                Write(pd_umon,*)"non existent file so a new one was created !! "
             End If

          Else If ( (status =='old') .AND. (tree%streams%dim_st(ii) > 0) ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_tree' file-open status ", status
             Write(pd_umon,*)"is inconsistent with a non existent file and stream dimension > 0 !!"
             Write(pd_umon,*)"Please check whether there's file Number ",ii," missing "
             write(pd_umon,*)trim(tree%streams%stream_files(ii))
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          Else if ( (status /= 'new') .and. (status /= 'old') .and. (status /= 'replace') ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_tree' file-open status ", status
             Write(pd_umon,*)"is not a valid parameter !!"
             Write(pd_umon,*)"Only 'new', 'old' and 'replace' are supported"
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          End If

       End If

    End Do

  End Subroutine open_stream_files_from_tree

  !============================================================================
  !> Subroutine which connects stream files
  !> 
  !> The subroutine opens all stream files that belong to a puredat tree 
  !> structure it takes the followig parametres IN SMALL LETERS for action and 
  !> status :
  !> ACTION : "read", "write"
  !> STATUS : "old", "new", "replace" 
  !> They are passed to the corresponding prameters of the Fortran open 
  !> statement
  Subroutine open_stream_files_from_streams(streams, action, status, position)

    Type(tStreams)  , Intent(InOut)              :: streams
    Character(len=*), intent(In)                 :: action, status

    Character(len=*), intent(In),optional        :: position
    Character(len=pd_mcl)                        :: lpos

    Character(len=10)                            :: faction

    logical :: fexist
    Integer :: ii, funit

    If (present(position)) then
       lpos = position
    Else
       lpos = "REWIND"
    End If

    Do ii = 1, no_streams

       Inquire(file=trim(streams%stream_files(ii)), exist=fexist, &
            number=funit, action=faction) 

       !** File exists ********************************************************
       If (fexist) Then

          !** File is connected to a unit *************************************
          If (funit > -1) then

             !** Check units integrity ***********************************
             if (streams%units(ii) /= funit) then

                Write(pd_umon,*)"In Subroutine open_stream_files_from_streams'"
                Write(pd_umon,*)"An inconsistency was detected in units"
                Write(pd_umon,*)"Expected unit no. is     :",streams%units(ii)
                write(pd_umon,*)"Inqure returned unit no. :",funit
                Write(pd_umon,*)"PROGRAM STOPPED !!"
                stop

             End if

             !** Check whether the file is opened for the desired action ******
             if (action /= trim(faction)) then

                close(streams%units(ii))
                Open(unit=streams%units(ii), &
                     file=streams%stream_files(ii), status=status, &
                     action=action, access='stream', form='unformatted',&
                     position=trim(lpos))

             End if

             !** File is not connected to a unit *********************************
          Else

             streams%units(ii) = pd_give_new_unit()
             Open(unit=streams%units(ii), &
                  file=streams%stream_files(ii), status=status, &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             streams%ifopen(ii) = .TRUE.

          End If

          !** File does not exist ************************************************
       Else 

          If ((status == 'new') .or. (status == 'replace')) then

             streams%units(ii) = pd_give_new_unit()
             Open(unit=streams%units(ii), &
                  file=streams%stream_files(ii), status=status, &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open action ", action
                Write(pd_umon,*)"is used on an empty file !!"
             End If


          Else If ( (status =='old') .AND. ( trim(lpos) == 'append' ) ) then

             streams%units(ii) = pd_give_new_unit()
             Open(unit=streams%units(ii), &
                  file=streams%stream_files(ii) , status='new', &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine 'open_stream_files_from_streams'"
                Write(pd_umon,*)"file-open status  ", status," ws used on a    "
                Write(pd_umon,*)"non existent file so a new one was created !! "
             End If

          Else If ( (status =='old') .AND. (streams%dim_st(ii) > 0) ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open status ", status
             Write(pd_umon,*)"is inconsistent with a non existent file and stream dimension > 0 !!"
             Write(pd_umon,*)"Please check whether there's file Number ",ii," missing "
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          Else if ( (status /= 'new') .and. (status /= 'old') .and. (status /= 'replace') ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open status ", status
             Write(pd_umon,*)"is not a valid parameter !!"
             Write(pd_umon,*)"Only 'new', 'old' and 'replace' are supported"
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          End If

       End If

    End Do

  End Subroutine open_stream_files_from_streams

  !============================================================================
  !> Subroutine which connects stream files
  !> 
  !> The subroutine opens all stream files that belong to a puredat tree 
  !> structure with MPI_FILE_OPEN.
  !> It takes the followig parametres IN SMALL LETERS for action and 
  !> status :
  !> ACTION : "read", "write"
  !> STATUS : "old", "new", "replace" 
  !> They are transformed to the corresponding prameters of the MPI open 
  !> statement
  !> Remark: Subroutine is in pre alpha state !!!
  Subroutine open_stream_files_from_streams_mpi(streams, action, status, fh_mpi, position)

    Type(tStreams)  , Intent(InOut)              :: streams
    Character(len=*), intent(In)                 :: action, status

    Integer(kind=pd_mpi_ik), Intent(InOut), Dimension(no_streams) :: fh_mpi
    Integer(kind=pd_mpi_ik)                      :: ierr
    
    Character(len=*), intent(In),optional        :: position
    Character(len=pd_mcl)                        :: lpos

    Character(len=10)                            :: faction

    logical :: fexist
    Integer :: ii, funit

    If (present(position)) then
       !lpos = position
       Write(*,*)"Sorry, the position parameter is not yet implemented"
    Else
       lpos = "REWIND"
    End If

    Do ii = 1, no_streams

!!$       Inquire(file=trim(streams%stream_files(ii)), exist=fexist, &
!!$            number=funit, action=faction) 
       fexist = .FALSE.
       
       !** File exists ********************************************************
       If (fexist) Then

          !** File is connected to a unit *************************************
          If (funit > -1) then

             !** Check units integrity ***********************************
             if (streams%units(ii) /= funit) then

                Write(pd_umon,*)"In Subroutine open_stream_files_from_streams'"
                Write(pd_umon,*)"An inconsistency was detected in units"
                Write(pd_umon,*)"Expected unit no. is     :",streams%units(ii)
                write(pd_umon,*)"Inqure returned unit no. :",funit
                Write(pd_umon,*)"PROGRAM STOPPED !!"
                stop

             End if

             !** Check whether the file is opened for the desired action ******
             if (action /= trim(faction)) then

                close(streams%units(ii))
                Open(unit=streams%units(ii), &
                     file=streams%stream_files(ii), status=status, &
                     action=action, access='stream', form='unformatted',&
                     position=trim(lpos))

             End if

             !** File is not connected to a unit *********************************
          Else

             streams%units(ii) = pd_give_new_unit()
             Open(unit=streams%units(ii), &
                  file=streams%stream_files(ii), status=status, &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             streams%ifopen(ii) = .TRUE.

          End If

          !** File does not exist ************************************************
       Else 

          If ((status == 'new') .or. (status == 'replace')) then

             Call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(streams%stream_files(ii)), &
                  MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, FH_MPI(ii), ierr)

             If (ierr /= 0) call file_err(trim(streams%stream_files(ii)), Int(ierr,pd_ik), &
                  "MPI_FILE_OPEN", "open_stream_files_from_streams_mpi")
             
             Select Case (ii)
             
             Case (1)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_INTEGER1, MPI_INTEGER1, &
                     "native", MPI_INFO_NULL, IERR)
                
             Case (2)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_INTEGER2, MPI_INTEGER2, &
                     "native", MPI_INFO_NULL, IERR)

             Case (3)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_INTEGER4, MPI_INTEGER4, &
                     "native", MPI_INFO_NULL, IERR)

             Case (4)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_INTEGER8, MPI_INTEGER8, &
                     "native", MPI_INFO_NULL, IERR)

             Case (5)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_Real8   , MPI_Real8   , &
                     "native", MPI_INFO_NULL, IERR)

             Case (6)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_Character, MPI_Character, &
                     "native", MPI_INFO_NULL, IERR)

             Case (7)
                Call MPI_FILE_SET_VIEW(FH_MPI(ii), 0_MPI_OFFSET_KIND, &
                     MPI_Character, MPI_Character, &
                     "native", MPI_INFO_NULL, IERR)   

             End Select
          
             If (ierr /= 0) call file_err(trim(streams%stream_files(ii)), Int(ierr,pd_ik), &
                   "MPI_FILE_SET_VIEW", "open_stream_files_from_streams_mpi")
             
             streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open action ", action
                Write(pd_umon,*)"is used on an empty file !!"
             End If


          Else If ( (status =='old') .AND. ( trim(lpos) == 'append' ) ) then

             streams%units(ii) = pd_give_new_unit()
             Open(unit=streams%units(ii), &
                  file=streams%stream_files(ii) , status='new', &
                  action=action, access='stream', form='unformatted', &
                  position=trim(lpos))
             streams%ifopen(ii) = .TRUE.

             If (action == 'read') then
                Write(pd_umon,*)"In Subroutine 'open_stream_files_from_streams'"
                Write(pd_umon,*)"file-open status  ", status," ws used on a    "
                Write(pd_umon,*)"non existent file so a new one was created !! "
             End If

          Else If ( (status =='old') .AND. (streams%dim_st(ii) > 0) ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open status ", status
             Write(pd_umon,*)"is inconsistent with a non existent file and stream dimension > 0 !!"
             Write(pd_umon,*)"Please check whether there's file Number ",ii," missing "
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          Else if ( (status /= 'new') .and. (status /= 'old') .and. (status /= 'replace') ) then

             Write(pd_umon,*)"In Subroutine open_stream_files_from_streams' file-open status ", status
             Write(pd_umon,*)"is not a valid parameter !!"
             Write(pd_umon,*)"Only 'new', 'old' and 'replace' are supported"
             Write(pd_umon,*)"PROGRAM STOPPED !!"
             stop

          End If

       End If

    End Do

  End Subroutine open_stream_files_from_streams_mpi

  !============================================================================
  !> Subroutine which disconnects stream files
  !> 
  !> The subroutine closes all stream files that belong to a puredat tree 
  !> structure
  Subroutine close_stream_files_from_tree(tree,clean)

    Type(tBranch), Intent(InOut)           :: tree
    Logical      , intent(in)   , optional :: clean

    Logical :: loc_clean

    Integer :: ii

    If (present(clean)) then
       loc_clean = clean
    Else
       loc_clean = .FALSE.
    end If

    If (.NOT. allocated(tree%streams)) then
       call alloc_error(-1, "streams", "close_stream_files_from_tree")
    End If

    Do ii = 1, no_streams

       If (tree%streams%ifopen(ii) .AND. (.NOT.loc_clean)) then
          Close(tree%streams%units(ii))
          tree%streams%units(ii) = -1
          tree%streams%ifopen(ii) = .FALSE.
       Else if (tree%streams%ifopen(ii) .AND. (loc_clean)) then
          If (tree%streams%dim_st(ii) > 0) then
             Close(tree%streams%units(ii))
          Else
             Close(tree%streams%units(ii), status="delete")
          End If
          tree%streams%units(ii) = -1
          tree%streams%ifopen(ii) = .FALSE.
       End If

    End Do

  End Subroutine close_stream_files_from_tree

  !============================================================================
  !> Subroutine which disconnects stream files
  !> 
  !> The subroutine closes all stream files that belong to a puredat tree 
  !> structure
  Subroutine close_stream_files_from_streams(streams,clean)

    Type(tStreams), Intent(InOut)              :: streams
    Logical       , intent(in)   , optional    :: clean

    Logical :: loc_clean

    Integer :: ii

    If (present(clean)) then
       loc_clean = clean
    Else
       loc_clean = .FALSE.
    end If

    Do ii = 1, no_streams

       If (streams%ifopen(ii) .AND. (.NOT.loc_clean)) then
          Close(streams%units(ii))
          streams%units(ii) = -1
          streams%ifopen(ii) = .FALSE.
       Else if (streams%ifopen(ii) .AND. (loc_clean)) then
          If (streams%dim_st(ii) > 0) then
             Close(streams%units(ii))
          Else
             Close(streams%units(ii), status="delete")
          End If
          streams%units(ii) = -1
          streams%ifopen(ii) = .FALSE.
       End If

    End Do

  End Subroutine close_stream_files_from_streams

  !============================================================================
  !> Function which writes tStreaam comonents to disk
  Subroutine Write_Streams(streams)

    Type(tStreams), intent(In) :: streams

    If ( streams%dim_st(1) > 0 ) then
       Write(streams%units(1)) streams%int1_st
    End If

    If ( streams%dim_st(2) > 0 ) then
       Write(streams%units(2)) streams%int2_st
    End If

    If ( streams%dim_st(3) > 0 ) then
       Write(streams%units(3)) streams%int4_st
    End If

    If ( streams%dim_st(4) > 0 ) then
       Write(streams%units(4)) streams%int8_st
    End If
    
    If ( streams%dim_st(5) > 0 ) then
       Write(streams%units(5)) streams%real8_st
    End If
    
    If ( streams%dim_st(6) > 0 ) then
       Write(streams%units(6)) streams%char_st
    End If
    
    If ( streams%dim_st(7) > 0 ) then
       Write(streams%units(7)) streams%log_st
    End If
    
  End Subroutine Write_Streams
  !> @} 
  !# End of memeber group "Puredat stream file routines" ######################

  !############################################################################
  !> \name Puredat load routines
  !> @{
  !> Module procedures for accessing data from a puredat tree structure 

  !============================================================================
  !> Function which retrieves integer 8 leaf data to a scalar
  !>
  !> Function which retrieves integer 8 leaf data from the corresponding leaf
  !> pointer to a scalar of type Integer(kind=8)
  Subroutine pd_get_4_scalar(branch, desc, values)

    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=8)  , intent(out)                        :: values

    Integer              :: ii, no
    Logical              :: desc_found

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          values     = branch%leaves(ii)%p_int8(1)
          desc_found = .TRUE.
          no         = ii 

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'( A)')"Something bad and unexpected happend during retrival&
            & of leaf data pointer"
       Write(pd_umon,'( A)')trim(desc)
       Write(pd_umon,'(3A)')"In puredat function pd_get_4_scalar(",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'( A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

    If (branch%leaves(no)%dat_no > 1) then
       Write(pd_umon,'( A)')"WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       Write(pd_umon,'(3A)')"In puredat function pd_get_4_scalar(",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"leaf%dat_no > 1 for scalar retrival of leaf data"
       Write(pd_umon,'( A)')"You will only get the first data element as a scalar"
    End If

  End Subroutine pd_get_4_scalar
 
  !============================================================================
  !> Function which retrieves real 8 leaf data to a scalar
  !>
  !> Function which retrieves real 8 leaf data from the corresponding leaf
  !> pointer to a scalar of type Real(kind=8)
  Subroutine pd_get_5_scalar(branch, desc, values)

    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Real(kind=8)     , intent(out)                        :: values

    Integer              :: ii, no
    Logical              :: desc_found

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          values     = branch%leaves(ii)%p_real8(1)
          desc_found = .TRUE.
          no         = ii 

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'( A)')"Something bad and unexpected happend during retrival&
            & of leaf data pointer"
       Write(pd_umon,'( A)')trim(desc)
       Write(pd_umon,'(3A)')"In puredat function pd_get_5_scalar(",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'( A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

    If (branch%leaves(no)%dat_no > 1) then
       Write(pd_umon,'( A)')"WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       Write(pd_umon,'(5A)')"In puredat function pd_get_5_scalar(",TRIM(desc)," from ",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"leaf%dat_no > 1 for scalar retrival of leaf data"
       Write(pd_umon,'( A)')"You will only get the first data element as a scalar"
    End If

  End Subroutine pd_get_5_scalar

  !============================================================================
  !> Function which retrieves integer 8 leaf data to an allocatable array
  !>
  !> Function which retrieves integer 8 leaf data from the corresponding leaf
  !> pointer to an allocatable array of rank 1 and type Integer(kind=8). The 
  !> leaf is searched !!non recursively!! in the given branch by full 
  !> comparison of trim(branch%leaves(ii)%desc) == trim(desc).
  Subroutine pd_get_4(branch, desc, values)

    Type(tBranch)    , Intent(in)                             :: branch
    Character(Len=*) , Intent(in)                             :: desc

    Integer(kind=8) , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found

    desc_found=.FALSE.
    
    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          Call alloc_error(alloc_stat,'values', 'pd_get_4', branch%leaves(ii)%dat_no)
          values = branch%leaves(ii)%p_int8
          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during retrival&
            & of leaf data pointer"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat function pd_get_4(branch,desc)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_get_4

  !============================================================================
  !> Function which retrieves integer 8 leaf data to a constant array
  !>
  !> Function which retrieves integer 8 leaf data from the corresponding leaf
  !> pointer to a constant array of rank 1 and type Integer(kind=8). The 
  !> leaf is searched !!non recursively!! in the given branch by full 
  !> comparison of trim(branch%leaves(ii)%desc) == trim(desc).
  Subroutine pd_get_4_vector(branch, desc, values, size)

    Type(tBranch)    , Intent(in)                  :: branch
    Character(Len=*) , Intent(in)                  :: desc
    Integer(kind=8) ,  Intent(in)                  :: size
    
    Integer(kind=8) , Intent(out), Dimension(size) :: values

    Integer              :: ii
    Logical              :: desc_found

    desc_found=.FALSE.
    
    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then
          
          If (size /= branch%leaves(ii)%dat_no)  then
             Write(pd_umon,'( A)')"Something bad and unexpected happend during retrival&
                  & of leaf data pointer"
             Write(pd_umon,'( A)')trim(desc)
             Write(pd_umon,'(3A)')"In puredat function pd_get_4(",TRIM(branch%desc),")"
             Write(pd_umon,'( A)')"The specified leaf size is not equal to the size  !!!"
             Write(pd_umon,'( A)')"of the passed actual argument                     !!!"
             Write(pd_umon,'( A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             stop
          End If
          
          values = branch%leaves(ii)%p_int8
          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during retrival&
            & of leaf data pointer"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat function pd_get_4(branch,desc)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_get_4_vector
  
  !============================================================================
  !> Function which retrieves real 8 leaf data to an allocatable array
  !>
  !> Function which retrieves real 8 leaf data from the corresponding leaf
  !> pointer to an allocatable array of rank 1 and type Real(kind=8). The leaf
  !> is searched !!non recursively!! in the given branch by full comparison
  !> of trim(branch%leaves(ii)%desc) == trim(desc).
  Subroutine pd_get_5(branch, desc, values)

    Type(tBranch)    , Intent(in)                             :: branch
    Character(Len=*) , Intent(in)                             :: desc

    Real(kind=8)     , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found

    desc_found=.FALSE.
    
    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          Call alloc_error(alloc_stat,'values' , 'pd_get_5', branch%leaves(ii)%dat_no)
          values = branch%leaves(ii)%p_real8
          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during retrival&
            & of leaf data values"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat function pd_get_5(branch,desc)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_get_5

  !============================================================================
  !> Function which retrieves character leaf data to an allocatable array
  !>
  !> Function which retrieves character leaf data from the corresponding leaf
  !> pointer to an allocatable array of rank 1 and type Character. The leaf
  !> is searched !!non recursively!! in the given branch by full comparison
  !> of trim(branch%leaves(ii)%desc) == trim(desc).
  Subroutine pd_get_6(branch, desc, values)

    Type(tBranch)    , Intent(in)                             :: branch
    Character(Len=*) , Intent(in)                             :: desc

    Character        , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          Call alloc_error(alloc_stat,'values' ,&
               'pd_get_5', branch%leaves(ii)%dat_no)
          values = branch%leaves(ii)%p_char
          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'( A)')"Something bad and unexpected happend during retrival&
            & of leaf data values"
       Write(pd_umon,'( A)')trim(desc)
       Write(pd_umon,'(3A)')"In puredat function pd_get_6(",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'( A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_get_6

  !============================================================================
  !> Function which retrieves integer 8 leaf data to a constant array
  !>
  !> Function which retrieves integer 8 leaf data from the corresponding leaf
  !> pointer to a constant array of rank 1, type Integer(kind=8) and size
  !> leaf%dat_no
  Subroutine pd_get_leaf_4_const(leaf, values)

    Type(tLeaf)     , Intent(in)                          :: leaf
    Integer(kind=8) , Intent(out), Dimension(leaf%dat_no) :: values

    values = leaf%p_int8
          
  End Subroutine pd_get_leaf_4_const

  !============================================================================
  !> Function which retrieves integer 8 leaf data to a constant array
  !>
  !> Function which retrieves integer 8 leaf data from the corresponding leaf
  !> pointer to a constant array of rank 1, type Integer(kind=8) and size
  !> leaf%dat_no
  Subroutine pd_get_leaf_4_2D_const(leaf, values)

    Type(tLeaf)     , Intent(in)                   :: leaf
    Integer(kind=8) , Intent(out), Dimension(:,:)  :: values

    values = reshape(leaf%p_int8,shape(values))
          
  End Subroutine pd_get_leaf_4_2D_const
  
  !============================================================================
  !> Function which retrieves real 8 leaf data to a constant array
  !>
  !> Function which retrieves real 8 leaf data from the corresponding leaf
  !> pointer to a constant array of rank 1, type Real(kind=8) and size
  !> leaf%dat_no
  Subroutine pd_get_leaf_5_const(leaf, values)

    Type(tLeaf)     , Intent(in)                       :: leaf
    Real(kind=8) , Intent(out), Dimension(leaf%dat_no) :: values

    values = leaf%p_real8
          
  End Subroutine pd_get_leaf_5_const

  !============================================================================
  !> Function which retrieves a leaf pointer
  !>
  !> Function which retrieves a leaf pointer by non recursive comparison of
  !> desc to branch%leaves(ii)%desc.
  Subroutine pd_get_leaf_pointer(branch, desc, leaf_p)

    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Type(tLeaf)      , intent(out), Pointer               :: leaf_p

    Integer              :: ii, no
    Logical              :: desc_found

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          leaf_p     => branch%leaves(ii)
          desc_found = .TRUE.
          no         = ii 

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'( A)')"Something bad and unexpected happend during retrival&
            & of leaf data pointer"
       Write(pd_umon,'( A)')trim(desc)
       Write(pd_umon,'(3A)')"In puredat function pd_get_4_scalar(",TRIM(branch%desc),")"
       Write(pd_umon,'( A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'( A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_get_leaf_pointer

  !> @} 
  !# End of memeber group "Puredat load routines" #############################

  !############################################################################
  !> \name Puredat read routines
  !> @{
  !> Module procedures for reading a puredat tree structure and its data from
  !> disc

  !============================================================================
  !> subroutine which reads a tree from header file
  Subroutine sub_read_tree(tree,success)

    Type(tBranch)        , Intent(Out) :: tree
    Logical , Optional   , Intent(out) :: success

    Integer                            :: io_stat=0
    Integer                            :: un_head
    Logical                            :: fexist
    !---------------------------------------------------------------------

    call raise_tree('',tree)

    Inquire(file=Trim(pro_path)//Trim(pro_name)//'.head', exist=fexist)

    If (fexist) then

       un_head = pd_give_new_unit()
       Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.head', &
            status='old', action='read',iostat=io_stat)
       call file_err(Trim(pro_path)//Trim(pro_name)//'.head',io_stat)

       if (present(success)) then
          success = .TRUE.
          call read_branch_ws(un_head,tree,success)
       Else
!          call read_branch(un_head,tree)
          write(pd_umon,PDF_E_A)"call read_branch(un_head,tree) is not implemented"
          write(pd_umon,PDF_E_STOP)
          stop
       End if

       Close(un_head)

    Else

       call file_err(Trim(pro_path)//Trim(pro_name)//'.head',50000)

    End if

  End Subroutine sub_read_tree

  !============================================================================
  !> Function which reads a tree from header file
  Function read_tree(streams) Result(tree)

    Type(tStreams), Intent(inout), optional :: streams

    Type(tBranch)         :: tree
    Integer               :: io_stat=0
    Integer               :: un_head, alloc_stat
    Logical               :: fexist
    Integer(kind=pd_ik)   :: fsize, pos
    character, dimension(:), Allocatable :: head
    
    !---------------------------------------------------------------------

    call raise_tree('',tree)

    Inquire(file=Trim(pro_path)//Trim(pro_name)//'.head', &
         exist=fexist, size=fsize)

    If (fexist) then

       un_head = pd_give_new_unit()
       Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.head', &
            status='old', action='read',iostat=io_stat, access="stream")
       call file_err(Trim(pro_path)//Trim(pro_name)//'.head',io_stat)

       Allocate(head(fsize), Stat=alloc_stat)
       Call alloc_error(alloc_stat,'head', 'read_tree', fsize)

       read(un_head) head

       Close(un_head)
       pos = 1
       if (present(streams)) then
          streams%no_branches = 1
          call read_branch(head,tree,fsize,pos,streams)
       Else
          call read_branch(head,tree,fsize,pos)
       End if

    Else

       call file_err(Trim(pro_path)//Trim(pro_name)//'.head',50000)

    End if

  End Function read_tree

  !*********************************************************
  !** Subroutine for reading a branch form header file *****
  Recursive Subroutine read_branch(head,branch,size,pos,streams)

    character, dimension(:), intent(in)     :: head
    Type(tBranch), Intent(inOut)            :: branch
    Integer(kind=pd_ik), Intent(inOut)      :: pos
    Integer(kind=pd_ik), Intent(in)         :: size
    Type(tStreams), Intent(inout), Optional :: streams

    Integer               :: alloc_stat, ii
    Character(Len=pd_mcl) :: tmp_char
    Logical               :: streams_allocated

    Character(len=pd_mcl)    :: line
    
    line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
    pos=pos+len_trim(line)+1

    line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, branch%desc

    line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, branch%no_branches

    line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, branch%no_leaves

    if (present(streams)) then
       streams%no_branches = streams%no_branches + branch%no_branches
       streams%no_leaves   = streams%no_leaves   + branch%no_leaves
    ENd if

    Allocate(branch%branches(branch%no_branches), Stat=alloc_stat)
    Call alloc_error(alloc_stat,'branch%branches' ,&
         'read_branch', branch%no_branches)

    Allocate(branch%leaves(branch%no_leaves), Stat=alloc_stat)
    Call alloc_error(alloc_stat,'branch%leaves' ,&
         'read_branch', branch%no_leaves)

    line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, streams_allocated

    If (streams_allocated) then
       If (.NOT.allocated(branch%streams)) then         
          allocate(branch%streams, Stat=alloc_stat)
          Call alloc_error(alloc_stat,'branch%streams' ,&
               'read_branch', branch%no_branches)
       End If

       Do ii = 1, no_streams
          line = char_to_str(head(pos:min(pos+pd_mcl-1,size)))
          pos=pos+len_trim(line)+1
          Read(line,*)tmp_char, branch%streams%dim_st(ii)
       End Do
       branch%streams%ii_st = branch%streams%dim_st + 1
       
    End If

    Do ii = 1, branch%no_leaves
       branch%leaves(ii) = read_leaf(head,size,pos)
    End Do

    if (present(streams)) then
       Do ii = 1, branch%no_branches
          call read_branch(head,branch%branches(ii),size,pos,streams)
       End Do
    Else
       Do ii = 1, branch%no_branches
          call read_branch(head,branch%branches(ii),size,pos)
       End Do
    End if

  End Subroutine read_branch

  !*********************************************************
  !** Subroutine for reading a branch form header file *****
  Recursive Subroutine read_branch_ws(un_head,branch,success)

    Integer      , Intent(in)     :: un_head
    Type(tBranch), Intent(inOut)  :: branch
    Logical      , Intent(inOut)  :: success

    Integer               :: alloc_stat, ii, io_stat
    Character(Len=pd_mcl) :: tmp_char
    Logical               :: streams_allocated

    Read(un_head,*,iostat=io_stat)
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, branch%desc
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, branch%no_branches
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, branch%no_leaves
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if

    Allocate(branch%branches(branch%no_branches), Stat=alloc_stat)
    Call alloc_error(alloc_stat,'branch%branches' ,&
         'read_branch', branch%no_branches)

    Allocate(branch%leaves(branch%no_leaves), Stat=alloc_stat)
    Call alloc_error(alloc_stat,'branch%leaves' ,&
         'read_branch', branch%no_leaves)

    Read(un_head,*,iostat=io_stat)tmp_char, streams_allocated
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if

    If (streams_allocated) then
       If (.NOT.allocated(branch%streams)) then         
          allocate(branch%streams, Stat=alloc_stat)
          Call alloc_error(alloc_stat,'branch%streams' ,&
               'read_branch', branch%no_branches)
       End If

       Do ii = 1, no_streams
          Read(un_head,*,iostat=io_stat)tmp_char, branch%streams%dim_st(ii)
          if(io_stat /= 0) then
             success=.FALSE.
             goto 1000
          End if
       End Do
       branch%streams%ii_st = branch%streams%dim_st + 1

    End If

    Do ii = 1, branch%no_leaves
       branch%leaves(ii) = read_leaf_ws(un_head,success)
    End Do

    if(success) then
       Do ii = 1, branch%no_branches
          call read_branch_ws(un_head,branch%branches(ii),success)
          If (.not.success) exit
       End Do
    End if

1000 continue

  End Subroutine read_branch_ws

  !*********************************************************
  !** Subroutine for reading a leaf form header file *******
  Function read_leaf_ws(un_head,success) Result(leaf)

    Integer      , Intent(in)      :: un_head
    Logical      , intent(out)     :: success
    Type(tLeaf)                    :: leaf

    Character(Len=pd_mcl) :: tmp_char
    integer               :: io_stat

    Read(un_head,*,iostat=io_stat)
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, leaf%desc
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, leaf%dat_no
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, leaf%dat_ty
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, leaf%lbound
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    Read(un_head,*,iostat=io_stat)tmp_char, leaf%ubound
    if(io_stat /= 0) then
       success=.FALSE.
       goto 1000
    End if
    
    leaf%p_int1  => null()
    leaf%p_int2  => null()
    leaf%p_int4  => null()
    leaf%p_int8  => null()
    leaf%p_real8 => null()
    leaf%p_char  => null()
    leaf%p_log   => null()

    leaf%pstat = 0

1000 Continue

  End Function read_leaf_ws

  !*********************************************************
  !** Subroutine for reading a leaf form header file *******
  Function read_leaf(head,size,pos) Result(leaf)

    Character, Dimension(:), Intent(in) :: head
    Integer(kind=pd_ik), Intent(in)     :: size
    Integer(kind=pd_ik), Intent(inout)  :: pos
    Type(tLeaf)                         :: leaf
    
    Character(Len=pd_mcl) :: tmp_char, line

    Integer(kind=pd_ik)   :: end

    end = min(pos+6*pd_mcl,size)

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, leaf%desc

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, leaf%dat_no

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, leaf%dat_ty

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, leaf%lbound

    line = char_to_str(head(pos:end))
    pos=pos+len_trim(line)+1
    Read(line,*)tmp_char, leaf%ubound

    leaf%p_int1  => null()
    leaf%p_int2  => null()
    leaf%p_int4  => null()
    leaf%p_int8  => null()
    leaf%p_real8 => null()
    leaf%p_char  => null()
    leaf%p_log   => null()

    leaf%pstat = 0

  End Function read_leaf

  !============================================================================
  !> Function which retrieves all stream data from stream files
  subroutine read_streams_from_branch(branch)

    Type(tBranch)    , Intent(inout)                      :: branch
    Integer                                               :: alloc_stat

    CALL OPEN_STREAM_FILES(branch, "read", "old")

    !**************************************************************************
    If (branch%streams%dim_st(1) > 0) then

       If (associated(branch%streams%int1_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%int1_st allready associated")
       End If
       allocate(branch%streams%int1_st(branch%streams%dim_st(1)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "int1_st", "read_streams", &
                        branch%streams%dim_st(1))
       
       Read(branch%streams%units(1))branch%streams%int1_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(2) > 0) then

       If (associated(branch%streams%int2_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%int2_st allready associated")
       End If
       allocate(branch%streams%int2_st(branch%streams%dim_st(2)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "int2_st", "read_streams", &
                        branch%streams%dim_st(2))
       
       Read(branch%streams%units(2))branch%streams%int2_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(3) > 0) then

       If (associated(branch%streams%int4_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%int4_st allready associated")
       End If
       allocate(branch%streams%int4_st(branch%streams%dim_st(3)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "int4_st", "read_streams", &
                        branch%streams%dim_st(3))
       
       Read(branch%streams%units(3))branch%streams%int4_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(4) > 0) then

       If (associated(branch%streams%int8_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%int8_st allready associated")
       End If
       allocate(branch%streams%int8_st(branch%streams%dim_st(4)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "int8_st", "read_streams", &
                        branch%streams%dim_st(4))

       Read(branch%streams%units(4))branch%streams%int8_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(5) > 0) then

       If (associated(branch%streams%real8_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%real8_st allready associated")
       End If
       allocate(branch%streams%real8_st(branch%streams%dim_st(5)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "real8_st", "read_streams", &
                        branch%streams%dim_st(5))
       
       Read(branch%streams%units(5))branch%streams%real8_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(6) > 0) then

       If (associated(branch%streams%char_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%char_st allready associated")
       End If
       allocate(branch%streams%char_st(branch%streams%dim_st(6)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "char_st", "read_streams", &
                        branch%streams%dim_st(6))
       
       Read(branch%streams%units(6))branch%streams%char_st

    End If
    !**************************************************************************
    If (branch%streams%dim_st(7) > 0) then

       If (associated(branch%streams%log_st)) then
          call cons_error("read_streams_from_branch", &
               "branch%streams%log_st allready associated")
       End If
       allocate(branch%streams%log_st(branch%streams%dim_st(7)),&
                stat=alloc_stat)
       call alloc_error(alloc_stat, "log_st", "read_streams", &
                        branch%streams%dim_st(7))
       
       Read(branch%streams%units(7))branch%streams%log_st

    End If

    call close_stream_files(branch)

  End subroutine read_streams_from_branch

  !============================================================================
  !> Function which retrieves all stream data from stream files
  subroutine read_streams_from_streams(streams)

    Type(tStreams)    , Intent(inout)                      :: streams
    Integer                                               :: alloc_stat

    CALL OPEN_STREAM_FILES(streams, "read", "old")

    !**************************************************************************
    If (streams%dim_st(1) > 0) then

       If (.not.associated(streams%int1_st)) then
          allocate(streams%int1_st(streams%dim_st(1)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "int1_st", "read_streams", &
               streams%dim_st(1))
       End If
       Read(streams%units(1))streams%int1_st

    End If
    !**************************************************************************
    If (streams%dim_st(2) > 0) then

       If (.not.associated(streams%int2_st)) then
          allocate(streams%int2_st(streams%dim_st(2)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "int2_st", "read_streams", &
               streams%dim_st(2))
       End If
       Read(streams%units(2))streams%int2_st

    End If
    !**************************************************************************
    If (streams%dim_st(3) > 0) then

       If (.not.associated(streams%int4_st)) then
          allocate(streams%int4_st(streams%dim_st(3)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "int4_st", "read_streams", &
               streams%dim_st(3))
       End If
       Read(streams%units(3))streams%int4_st

    End If
    !**************************************************************************
    If (streams%dim_st(4) > 0) then

       If (.not.associated(streams%int8_st)) then
          allocate(streams%int8_st(streams%dim_st(4)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "int8_st", "read_streams", &
               streams%dim_st(4))
       End If
       Read(streams%units(4))streams%int8_st

    End If
    !**************************************************************************
    If (streams%dim_st(5) > 0) then

       If (.not.associated(streams%real8_st)) then
          allocate(streams%real8_st(streams%dim_st(5)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "real8_st", "read_streams", &
               streams%dim_st(5))
       End If
       Read(streams%units(5))streams%real8_st

    End If
    !**************************************************************************
    If (streams%dim_st(6) > 0) then

       If (.not.associated(streams%char_st)) then
          allocate(streams%char_st(streams%dim_st(6)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "char_st", "read_streams", &
               streams%dim_st(6))
       End If
       Read(streams%units(6))streams%char_st

    End If
    !**************************************************************************
    If (streams%dim_st(7) > 0) then

       If (.not.associated(streams%log_st)) then
          allocate(streams%log_st(streams%dim_st(7)),&
               stat=alloc_stat)
          call alloc_error(alloc_stat, "log_st", "read_streams", &
               streams%dim_st(7))
       End If
       Read(streams%units(7))streams%log_st

    End If

    call close_stream_files(streams)

  End subroutine read_streams_from_streams

  !============================================================================
  !> Function which retrieves integer 1 leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 1 (int1) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.

  Subroutine pd_load_leaf_1(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=1)  , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_1", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(1),pos=(branch%leaves(ii)%lbound-1)*1+1)&
               values

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_1&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_1

  !============================================================================
  !> Function which retrieves integer 2 leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 2 (int2) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.
  Subroutine pd_load_leaf_2(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=2)  , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_2", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(2),pos=(branch%leaves(ii)%lbound-1)*2+1)&
               values

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_2&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_2

  !============================================================================
  !> Function which retrieves integer 4 leaf data from stream file
  !>
  !> The Subroutine retrieves leaf data of type 3 (int4) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.
  Subroutine pd_load_leaf_3(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=4)  , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_3", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(3),pos=(branch%leaves(ii)%lbound-1)*4+1)&
               values

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_3&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_3

  !============================================================================
  !> Function which retrieves integer 8 leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 4 (int8) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.
  Subroutine pd_load_leaf_4(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=8)  , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_4", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(4),pos=(branch%leaves(ii)%lbound-1)*8+1)&
               values

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_4&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_4

  !============================================================================
  !> Function which retrieves real 8 leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 5 (real8) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.
  Subroutine pd_load_leaf_5(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Real(kind=8)     , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_5", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(5),pos=(branch%leaves(ii)%lbound-1)*8+1 )&
               values(1:branch%leaves(ii)%dat_no)

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_5&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_5

  !============================================================================
  !> Function which retrieves character leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 6 (char) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is allocated to the correct size and the read is directly performed
  !> on values.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.  
  Subroutine pd_load_leaf_6(streams,branch,desc, values)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Character        , Intent(out), Allocatable, Dimension(:) :: values

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          Allocate(values(branch%leaves(ii)%dat_no),stat=alloc_stat)
          call alloc_error(alloc_stat, "values", "pd_load_leaf_6", &
               branch%leaves(ii)%dat_no)
      
          read(streams%units(6),pos=(branch%leaves(ii)%lbound-1)*1+1 )&
               values

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_6&
            &(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_6

  !============================================================================
  !> Function which retrieves 2D int8 leaf data from stream files
  !>
  !> The Subroutine retrieves leaf data of type 4 (int8) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is of dimension 2 and is allocated to the correct size determined by
  !> branch%leaves(ii)%dat_no/factor. The read is directly performed
  !> on values. <br>
  !> !! If mod(branch%leaves(ii)%dat_no,factor) is not equal 0 a fatal error is
  !> issued and the program is halted.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.  
  Subroutine pd_load_leaf_4_2D(streams,branch,desc, values, factor)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Integer(kind=8)  , Intent(out), Allocatable, Dimension(:,:) :: values

    Integer(pd_ik)   , intent(in)                         :: factor

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          if (mod(branch%leaves(ii)%dat_no,factor) == 0) then

             Allocate(values(factor,branch%leaves(ii)%dat_no/factor),&
                  stat=alloc_stat)

             call alloc_error(alloc_stat, "values", "pd_load_leaf_4_2D", &
                  branch%leaves(ii)%dat_no)
      
             read(streams%units(4),pos=(branch%leaves(ii)%lbound-1)*8+1 )&
                  values
          Else
             Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
                  &leaf data"
             Write(pd_umon,'(A)')trim(desc)
             Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_4_2D&
                  &(streams,branch,desc,values,factor)"
             Write(pd_umon,'(A)')"The specified specified factor didn't lead to even"
             Write(pd_umon,'(A)')"integer division" 
             Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             stop
          End if

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_4_2D&
            &(streams,branch,desc,values,factor)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_4_2D

  !============================================================================
  !> Function which retrieves 2D real 8 leaf data from stream files 
  !>
  !> The Subroutine retrieves leaf data of type 5 (real8) from the corresponding
  !> stream file. The stream file has to be allready opened. The out parameter
  !> values is of dimension 2 and is allocated to the correct size determined by
  !> branch%leaves(ii)%dat_no/factor. The read is directly performed
  !> on values. <br>
  !> !! If mod(branch%leaves(ii)%dat_no,factor) is not equal 0 a fatal error is
  !> issued and the program is halted.<br>
  !> !! The passed in tBranch structure is NOT cycled recursively.  
  Subroutine pd_load_leaf_5_2D(streams,branch,desc, values, factor)

    Type(tStreams)   , Intent(in)                         :: streams
    Type(tBranch)    , Intent(in)                         :: branch
    Character(Len=*) , Intent(in)                         :: desc

    Real(kind=8)     , Intent(out), Allocatable, Dimension(:,:) :: values

    Integer(pd_ik)   , intent(in)                         :: factor

    Integer              :: ii, alloc_stat
    Logical              :: desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          if (mod(branch%leaves(ii)%dat_no,factor) == 0) then

             Allocate(values(factor,branch%leaves(ii)%dat_no/factor),&
                  stat=alloc_stat)
             call alloc_error(alloc_stat, "values", "pd_load_leaf_5_2D", &
                  branch%leaves(ii)%dat_no)
      
             read(streams%units(5),pos=(branch%leaves(ii)%lbound-1)*8+1 )&
                  values
          Else
             Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
                  &leaf data"
             Write(pd_umon,'(A)')trim(desc)
             Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_5_2D&
                  &(streams,branch,desc,values,factor)"
             Write(pd_umon,'(A)')"The specified specified factor didn't lead to even"
             Write(pd_umon,'(A)')"integer division" 
             Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             stop
          End if

          desc_found = .TRUE.

       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during load of &
            &leaf data"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_load_leaf_5_2D&
            &(streams,branch,desc,values,factor)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_load_leaf_5_2D

  !============================================================================
  !> Function which retrieves integer 1 leaf data from stream files
  Subroutine pd_read_leaf_1(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Integer(kind=1)  , Intent(out), Dimension(*) :: values

    read(streams%units(1),pos=(leaf%lbound-1)*1+1) values(1:leaf%dat_no)

  End Subroutine pd_read_leaf_1

  !============================================================================
  !> Function which retrieves integer 2 leaf data from stream files
  Subroutine pd_read_leaf_2(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Integer(kind=2)  , Intent(out), Dimension(*) :: values

    read(streams%units(2),pos=(leaf%lbound-1)*2+1) values(1:leaf%dat_no)

  End Subroutine pd_read_leaf_2

  !============================================================================
  !> Function which retrieves integer 4 leaf data from stream files
  Subroutine pd_read_leaf_3(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Integer(kind=4)  , Intent(out), Dimension(*) :: values

    read(streams%units(3),pos=(leaf%lbound-1)*4+1) values(1:leaf%dat_no)

  End Subroutine pd_read_leaf_3

  !============================================================================
  !> Function which retrieves integer 8 leaf data from stream files
  Subroutine pd_read_leaf_4(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Integer(kind=8)  , Intent(out), Dimension(:) :: values

    read(streams%units(4),pos=(leaf%lbound-1)*8+1) values

  End Subroutine pd_read_leaf_4

  !============================================================================
  !> Function which retrieves integer 8 leaf data from stream files
  Subroutine pd_read_leaf_4_scal(streams, leaf, value)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Integer(kind=8)  , Intent(out)               :: value

    read(streams%units(4),pos=(leaf%lbound-1)*8+1) value

  End Subroutine pd_read_leaf_4_scal

  !============================================================================
  !> Function which retrieves integer 8 leaf data from stream files
  Subroutine pd_read_leaf_4_2D(streams, leaf, values)

    Type(tStreams)   , Intent(in)                  :: streams
    Type(tLeaf)      , Intent(In)                  :: leaf
    Integer(kind=8)  , Intent(out), Dimension(:,:) :: values

    read(streams%units(4),pos=(leaf%lbound-1)*8+1) values

  End Subroutine pd_read_leaf_4_2D

  !============================================================================
  !> Function which retrieves real 8 leaf data from stream files
  Subroutine pd_read_leaf_5(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Real(kind=8)     , Intent(out), Dimension(:) :: values

    read(streams%units(5),pos=(leaf%lbound-1)*8+1) values

  End Subroutine pd_read_leaf_5

  !============================================================================
  !> Function which retrieves real 8 leaf data from stream files
  Subroutine pd_read_leaf_5_2D(streams, leaf, values)

    Type(tStreams)   , Intent(in)                  :: streams
    Type(tLeaf)      , Intent(In)                  :: leaf
    Real(kind=8)     , Intent(out), Dimension(:,:) :: values

    read(streams%units(5),pos=(leaf%lbound-1)*8+1) values

  End Subroutine pd_read_leaf_5_2D

  !============================================================================
  !> Function which retrieves character leaf data from stream files
  Subroutine pd_read_leaf_6(streams, leaf, values)

    Type(tStreams)   , Intent(in)                :: streams
    Type(tLeaf)      , Intent(In)                :: leaf
    Character        , Intent(out), Dimension(*) :: values

    read(streams%units(6),pos=(leaf%lbound-1)*1+1) values(1:leaf%dat_no)

  End Subroutine pd_read_leaf_6
  !> @} 
  !# End of memeber group "Puredat load routines" #############################

  !############################################################################
  !> \name Puredat storage routines
  !> @{
  !> Module procedures for storing data to puredat stream files along with 
  !> keeping the puredat tree structure up to date

  !============================================================================
  !> Function which writes a complete puredat tree to disk
  Subroutine write_tree(tree)

    Type(tBranch)   , Intent(In)    :: tree 
    Integer                         :: un_head
    !**********************************************************

    un_head = pd_give_new_unit()

    Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.head', status='replace', &
         action='write')

    Call write_branch(tree,un_head)

    Close(un_head)

  End Subroutine write_tree

  !==========================================================================
  !> Function which writes a purdat branch to disk
  Recursive Subroutine write_branch(branch,un_head)

    Type(tBranch)   , Intent(In) :: branch
    Integer         , Intent(In) :: un_head    

    Integer(kind=pd_ik) :: ii
    
    Write(un_head, fmt_bsep)
    Write(un_head, '(4(A))')'<description> ', "'", Trim(branch%desc), "'"
    Write(un_head, '(A,I0)')'<no_of_branches> ',branch%no_branches
    Write(un_head, '(A,I0)')'<no_of_leaves> ',branch%no_leaves
    
    If (allocated(branch%streams)) then
       Write(un_head, '(A,L1)')'<streams_allocated> ',.TRUE.
       Write(un_head, '(A,I0)')'<size_int1_stream> ',branch%streams%dim_st(1)
       Write(un_head, '(A,I0)')'<size_int2_stream> ',branch%streams%dim_st(2)
       Write(un_head, '(A,I0)')'<size_int4_stream> ',branch%streams%dim_st(3)
       Write(un_head, '(A,I0)')'<size_int8_stream> ',branch%streams%dim_st(4)
       Write(un_head, '(A,I0)')'<size_real_stream> ',branch%streams%dim_st(5)
       Write(un_head, '(A,I0)')'<size_char_stream> ',branch%streams%dim_st(6)
       Write(un_head, '(A,I0)')'<size_log_stream>  ',branch%streams%dim_st(7)
    Else
       Write(un_head, '(A,L1)')'<streams_allocated> ',.FALSE.
    End If
    
    Do ii = 1, branch%no_leaves
       Call write_leaf(branch%leaves(ii),un_head)
    End Do
    
    Do ii = 1, branch%no_branches
       Call write_branch(branch%branches(ii),un_head)
    End Do
    
  End Subroutine write_branch
  
  !==========================================================================
  !> Function which writes a purdat leaf to disc
  Subroutine write_leaf(leaf,un_head)
    
    Type(tLeaf)   , Intent(In) :: leaf 
    Integer         , Intent(In) :: un_head    

    Write(un_head, fmt_lsep)
    Write(un_head, '(4(A))')'<description> ' , "'", Trim(leaf%desc), "'"
    Write(un_head, '(A,I0)')'<no_of_data> '  ,leaf%dat_no
    Write(un_head, '(A,I0)')'<type_of_data> ',leaf%dat_ty
    Write(un_head, '(A,I0)')'<lower_bound> ' ,leaf%lbound
    Write(un_head, '(A,I0)')'<upper_bound> ' ,leaf%ubound
    
  End Subroutine write_leaf

  !============================================================================
  !> Function which stores integer 1 leaf data
  Subroutine  pd_store_1(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Integer(Kind=1)  , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 1 ) call cons_error("pd_store_1",&
               "Storing data other than Type 1 is not allowed with this routine")

          If (.not.loc_blind) then
             If ( streams%ifopen(1) ) then
                Write(streams%units(1))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_1",&
                     "Storing data to closed file is not possible with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(1)
          streams%ii_st(1) = streams%ii_st(1) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(1) - 1
          streams%dim_st(1) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
          
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_1(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_1

  !============================================================================
  !> Function which stores integer 2 leaf data
  Subroutine  pd_store_2(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Integer(kind=2)  , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 2 ) call cons_error("pd_store_2",&
               "Storing data other than Type 2 is not allowed with this routine")
          
          If (.not.loc_blind) then
             If ( streams%ifopen(2) ) then
                Write(streams%units(2))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_2",&
                     "Storing data to closed file is not possible with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(2)
          streams%ii_st(2) = streams%ii_st(2) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(2) - 1
          streams%dim_st(2) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_2(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_2

  !============================================================================
  !> Function which stores integer 4 leaf data
  Subroutine  pd_store_3(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Integer(kind=4)  , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 3 ) call cons_error("pd_store_3",&
               "Storing data other than Type 3 is not allowed with this routine")

          If (.not.loc_blind) then
             If ( streams%ifopen(3) ) then
                Write(streams%units(3))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_3",&
                     "Storing data to closed file is not possible with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(3)
          streams%ii_st(3) = streams%ii_st(3) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(3) - 1
          streams%dim_st(3) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_3(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_3

  !============================================================================
  !> Function which stores integer 8 leaf data
  Subroutine  pd_store_4(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Integer(kind=8)  , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 4 ) call cons_error("pd_store_4",&
               "Storing data other than Type 4 is not allowed with this routine")

          If (.not.loc_blind) then
             If ( streams%ifopen(4) ) then
                Write(streams%units(4))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_4",&
                     "Storing data to closed file is not possible with this routine")
             End If
          END If

          branch%leaves(ii)%lbound = streams%ii_st(4)
          streams%ii_st(4) = streams%ii_st(4) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(4) - 1
          streams%dim_st(4) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_4(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_4

  !============================================================================
  !> Function which stores integer 8 leaf data
  Subroutine  pd_store_4_2D(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Integer(kind=pd_ik)  , Intent(in), Dimension(:,:)       :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 4 ) call cons_error("pd_store_4",&
               "Storing data other than Type 4 is not allowed with this routine")

          if ( branch%leaves(ii)%dat_no /= size(values) ) &
               call cons_error("pd_store_4_2D", &
               "In leaf "//trim(desc)//char(10)//&
               &"branch%leaves(ii)%dat_no /= size(values) this should not happen !!")

          If (.not.loc_blind) then
             If ( streams%ifopen(4) ) then
                Write(streams%units(4))values
             else
                call cons_error("pd_store_4",&
                     "Storing data to closed file is not possible with this routine")
             End If
          END If

          branch%leaves(ii)%lbound = streams%ii_st(4)
          streams%ii_st(4) = streams%ii_st(4) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(4) - 1
          streams%dim_st(4) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_4_2D(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_4_2D

  !============================================================================
  !> Function which stores real 8 leaf data
  Subroutine  pd_store_5(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Real(kind=pd_rk) , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 5 ) call cons_error("pd_store_5",&
               "Storing data other than Type 5 is not allowed with this routine")
          
          If (.not.loc_blind) then
             If ( streams%ifopen(5) ) then
                Write(streams%units(5))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_5",&
                     "Storing data to closed file is not possible with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(5)
          streams%ii_st(5) = streams%ii_st(5) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(5) - 1
          streams%dim_st(5) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_5(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch with desc"
       Write(pd_umon,'(A)')trim(branch%desc)
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_5

  !============================================================================
  !> Function which stores real 8 leaf data
  Subroutine  pd_store_5_2D(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Real(kind=pd_rk) , Intent(in), Dimension(:,:)       :: values
    Logical          , intent(in), Optional             :: blind
    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If (branch%leaves(ii)%dat_no /= size(values)) & 
               call cons_error("pd_store_5_2D", &
               "In leaf "//trim(desc)//char(10)//&
               &"branch%leaves(ii)%dat_no /= size(values) this should not happen !!")

          If ( branch%leaves(ii)%dat_ty /= 5 ) call cons_error("pd_store_5",&
               "Storing data other than Type 5 is not allowed with this routine")
          
          If (.not.loc_blind) then
             If ( streams%ifopen(5) ) then
                Write(streams%units(5))values
             else
                call cons_error("pd_store_5",&
                     "Storing data to closed file is not possible with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(5)
          streams%ii_st(5) = streams%ii_st(5) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(5) - 1
          streams%dim_st(5) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_5_2D(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_5_2D

  !============================================================================
  !> Function which stores character leaf data
  Subroutine  pd_store_6(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Character        , Intent(in), Dimension(*)         :: values
    Logical          , intent(in), Optional             :: blind

    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 6 ) call cons_error("pd_store_6",&
               "Storing data other than Type 6 is not allowed with this routine")

          If (.not.loc_blind) then
             If ( streams%ifopen(6) ) then
                Write(streams%units(6))values(1:branch%leaves(ii)%dat_no)
             else
                call cons_error("pd_store_6",&
                     "Storing data to closed file is not possible "//&
                     "with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(6)
          streams%ii_st(6) = streams%ii_st(6) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(6) - 1
          streams%dim_st(6) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_6(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_6

  !============================================================================
  !> Function which stores character leaf data
  Subroutine  pd_store_6_str(streams,branch,desc,values,blind)

    Type(tStreams)   , Intent(inout)                    :: streams
    Type(tBranch)    , Intent(inout)                    :: branch
    Character(Len=*) , Intent(in)                       :: desc
    Character(Len=*) , Intent(in)                       :: values
    Logical          , intent(in), Optional             :: blind

    Integer              :: ii
    Logical              :: desc_found, loc_blind

    if(present(blind))then
       loc_blind=blind
    else
       loc_blind=.FALSE.
    End if

    desc_found=.FALSE.

    Do ii = 1, branch%no_leaves

       If (trim(branch%leaves(ii)%desc) == trim(desc)) then    

          If ( branch%leaves(ii)%dat_ty /= 6 ) call cons_error("pd_store_6",&
               "Storing data other than Type 6 is not allowed with this routine")

          If (.not.loc_blind) then
             If ( streams%ifopen(6) ) then
                Write(streams%units(6))values
             else
                call cons_error("pd_store_6",&
                     "Storing data to closed file is not possible "//&
                     "with this routine")
             End If
          End If

          branch%leaves(ii)%lbound = streams%ii_st(6)
          streams%ii_st(6) = streams%ii_st(6) + branch%leaves(ii)%dat_no
          branch%leaves(ii)%ubound = streams%ii_st(6) - 1
          streams%dim_st(6) = branch%leaves(ii)%ubound
          desc_found = .TRUE.
       End If

    end Do

    If (.NOT. desc_found) then
       Write(pd_umon,'(A)')"Something bad and unexpected happend during storage of data to leaf"
       Write(pd_umon,'(A)')trim(desc)
       Write(pd_umon,'(A)')"In puredat subroutine pd_store_6_str(streams,branch,desc,values)"
       Write(pd_umon,'(A)')"The specified leaf was not found in the passed branch" 
       Write(pd_umon,'(A)')"PROGRAM STOPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       stop
    End If

  End Subroutine pd_store_6_str

  !============================================================================
  !> Function which stores the data of a tBranch structure
  !>
  !> Function which stores the data of a tBranch structure with one direct
  !> mpi write per leaf.
  !> Remark: The routine is to be used carfully if the branch is deeply
  !> structured with only small data chunks in the leaves.
  Recursive Subroutine store_parallel_branch(br, FH_MPI)

    type(tBranch)          , Intent(in)                        :: br
    Integer(kind=pd_mpi_ik), Intent(in), Dimension(no_streams) :: fh_mpi

    Integer(kind=pd_mpi_ik)                                     :: ierr
    Integer(kind=pd_mpi_ik), Dimension(MPI_STATUS_SIZE)         :: status_mpi
    Integer(kind=pd_ik)                                         :: ii
    
    Do ii = 1, br%no_leaves

       if (br%leaves(ii)%pstat >= 0) then

          Select Case (br%leaves(ii)%dat_ty)
             
          Case (1)
             Call MPI_FILE_WRITE_AT(FH_MPI(1), &
                  Int(br%leaves(ii)%lbound-1, MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_int1, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_INTEGER1, &
                  status_mpi, ierr)
          Case (2)
             Call MPI_FILE_WRITE_AT(FH_MPI(2), &
                  Int((br%leaves(ii)%lbound-1), MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_int2, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_INTEGER2, &
                  status_mpi, ierr)
          Case (3)
             Call MPI_FILE_WRITE_AT(FH_MPI(3), &
                  Int((br%leaves(ii)%lbound-1), MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_int4, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_INTEGER4, &
                  status_mpi, ierr)
          Case (4)
             Call MPI_FILE_WRITE_AT(FH_MPI(4), &
                  Int((br%leaves(ii)%lbound-1), MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_int8, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_INTEGER8, &
                  status_mpi, ierr)
          Case (5)
             Call MPI_FILE_WRITE_AT(FH_MPI(5), &
                  Int((br%leaves(ii)%lbound-1), MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_real8, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_REAL8, &
                  status_mpi, ierr)
          Case (6)
             Call MPI_FILE_WRITE_AT(FH_MPI(6), &
                  Int(br%leaves(ii)%lbound-1, MPI_OFFSET_KIND), &
                  br%leaves(ii)%p_char, &
                  Int(br%leaves(ii)%dat_no,pd_mpi_ik), MPI_character, &
                  status_mpi, ierr)

          Case default
             Write(*,*)"Data type ",br%leaves(ii)%dat_ty," is not yet implemented"

          End Select

       End if

    End Do
    
    Do ii = 1, br%no_branches
       Call store_parallel_branch(br%branches(ii), FH_MPI)
    End Do
    
  End Subroutine store_parallel_branch

  !============================================================================
  !> Function which dumps a tbranch structure recursively to disk
  !>
  !> The data contained in the branch are dumped sequentially with a leaves
  !> first rule applied. A new tstreams structure is generated from which the
  !> resulting sequential leaf bounds are derived
  Subroutine dump_branch(branch)

    Type(tBranch)   , Intent(In)    :: branch
    Integer                         :: un_head
    Type(tBranch)                   :: tmp_tree
    Integer(kind=pd_ik)             :: ii, wskip
    Character(len=pd_mcl)           :: tmp_line
    !**********************************************************

    !** Raise temporary tree structure ***
    call raise_tree('',tmp_tree)
    call open_stream_files(tmp_tree%streams, "write", "new")
    
    un_head = pd_give_new_unit()

    Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.head', status='new', &
         action='write')

    Write(un_head, fmt_bsep)
    !------------------------12345678901234567890123456789012345678901234567890
    Write(un_head, '(4(A))')'<description> ', "'", Trim(branch%desc), "'"
    Write(un_head, '(A,I15)')'<no_of_branches> ',branch%no_branches
    Write(un_head, '(A,I15)')'<no_of_leaves> ',branch%no_leaves
    
    !** Initially dump the tstream component of tmp_tree ***
    !------------------------12345678901234567890123456789012345678901234567890
    Write(un_head, '(A,L1)')'<streams_allocated> ',.TRUE.
    Write(un_head, '(A,I15)')'<size_int1_stream> ',tmp_tree%streams%dim_st(1)
    Write(un_head, '(A,I15)')'<size_int2_stream> ',tmp_tree%streams%dim_st(2)
    Write(un_head, '(A,I15)')'<size_int4_stream> ',tmp_tree%streams%dim_st(3)
    Write(un_head, '(A,I15)')'<size_int8_stream> ',tmp_tree%streams%dim_st(4)
    Write(un_head, '(A,I15)')'<size_real_stream> ',tmp_tree%streams%dim_st(5)
    Write(un_head, '(A,I15)')'<size_char_stream> ',tmp_tree%streams%dim_st(6)
    Write(un_head, '(A,I15)')'<size_log_stream>  ',tmp_tree%streams%dim_st(7)
    
    Do ii = 1, branch%no_leaves
       Call dump_leaf(branch%leaves(ii),un_head, tmp_tree%streams)
    End Do
    
    Do ii = 1, branch%no_branches
       Call dump_branch_rec(branch%branches(ii),un_head, tmp_tree%streams)
    End Do

    Close(un_head)

    Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.head', status='old', &
         action='write',access="stream")

    !** Calc skip from begin of file :
    !** fmt_bsep               = 12 + 1
    !** branch desc. line      = 16 + len_trim(branch%desc) + 1
    !** no of branches line    = 17 + 15 + 1
    !** no of leaves line      = 15 + 15 + 1
    !** streams allocated line = 20 + 1 + 1
    wskip = 13 + 16+len_trim(branch%desc)+1 + 17+15+1 + 15+15+1 + 20+1+1 + 1
    
    !** Update the sizes of the tstream components of tmp_tree ***
    Write(tmp_line,'(A,I15,A)')'<size_int1_stream> ',tmp_tree%streams%dim_st(1), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_int2_stream> ',tmp_tree%streams%dim_st(2), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_int4_stream> ',tmp_tree%streams%dim_st(3), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_int8_stream> ',tmp_tree%streams%dim_st(4), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_real_stream> ',tmp_tree%streams%dim_st(5), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_char_stream> ',tmp_tree%streams%dim_st(6), char(10)
    Write(un_head, pos=wskip)trim(tmp_line)
    wskip = wskip + 35
    Write(tmp_line,'(A,I15,A)')'<size_log_stream>  ',tmp_tree%streams%dim_st(7)
    Write(un_head, pos=wskip)trim(tmp_line)

    close(un_head)
    
  End Subroutine dump_branch

  !==========================================================================
  !> Function dumps writes a purdat branch recursively to disk
  Recursive Subroutine dump_branch_rec(branch,un_head, streams)

    Type(tBranch)    , Intent(In)    :: branch
    Integer          , Intent(In)    :: un_head

    Type(tStreams)   , Intent(InOut) :: streams
    Integer(kind=pd_ik)              :: ii
    
    Write(un_head, fmt_bsep)
    Write(un_head, '(4(A))')'<description> ', "'", Trim(branch%desc), "'"
    Write(un_head, '(A,I0)')'<no_of_branches> ',branch%no_branches
    Write(un_head, '(A,I0)')'<no_of_leaves> ',branch%no_leaves

    !** We do not dump any contained tstream comonents ***
    Write(un_head, '(A,L1)')'<streams_allocated> ',.FALSE.
    
    Do ii = 1, branch%no_leaves
       Call dump_leaf(branch%leaves(ii),un_head, streams)
    End Do
    
    Do ii = 1, branch%no_branches
       Call dump_branch_rec(branch%branches(ii),un_head, streams)
    End Do
    
  End Subroutine dump_branch_rec
  
  !==========================================================================
  !> Function which dumpes a purdat leaf to disc
  Subroutine dump_leaf(leaf, un_head, streams)
    
    Type(tLeaf)      , Intent(In)    :: leaf 
    Integer          , Intent(In)    :: un_head    
    Type(tStreams)   , Intent(InOut) :: streams
    
    Write(un_head, fmt_lsep)
    Write(un_head, '(4(A))')'<description> ' , "'", Trim(leaf%desc), "'"
    Write(un_head, '(A,I0)')'<no_of_data> '  ,leaf%dat_no
    Write(un_head, '(A,I0)')'<type_of_data> ',leaf%dat_ty
    Write(un_head, '(A,I0)')'<lower_bound> ' ,streams%ii_st(leaf%dat_ty)
    Write(un_head, '(A,I0)')'<upper_bound> ' ,streams%ii_st(leaf%dat_ty) + &
                                              leaf%dat_no - 1
    
    streams%ii_st(leaf%dat_ty)  = streams%ii_st(leaf%dat_ty) + leaf%dat_no
    streams%dim_st(leaf%dat_ty) = streams%ii_st(leaf%dat_ty) - 1

    If ( leaf%pstat > 0 ) then
       
       Select case (leaf%dat_ty)

       Case(1)
          If (associated(leaf%p_int1)) then
             Write(streams%units(1))leaf%p_int1
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(2)
          If (associated(leaf%p_int2)) then
             Write(streams%units(2))leaf%p_int2
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(3)
          If (associated(leaf%p_int4)) then
             Write(streams%units(3))leaf%p_int4
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(4)
          If (associated(leaf%p_int8)) then
             Write(streams%units(4))leaf%p_int8
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(5)
          If (associated(leaf%p_real8)) then
             Write(streams%units(5))leaf%p_real8
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(6)
          If (associated(leaf%p_char)) then
             Write(streams%units(6))leaf%p_char
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       Case(7)
          If (associated(leaf%p_log)) then
             Write(streams%units(7))leaf%p_log
          Else
             write(pd_umon,PDF_W_A)"p_int1 not associated in leaf "//trim(leaf%desc)
          End If
          
       End Select
    End If
  End Subroutine dump_leaf
  !> @} 
  !# End of memeber group "Puredat storage routines" ##########################


  !############################################################################
  !> \name Puredat serialization routines
  !> @{
  !> Module procedures for serailizing puredat header structures
  
  !============================================================================
  !> Subroutine for serializing a puredat tree to arrays
  !>
  !> Subroutine which recursively  serializes a puredat tree to arrays 
  !> regardless of the branch structures.
  !> Rule: Leaves first !!
  Subroutine serialize_tree(tree, desc, dat_no, dat_ty, offsets)
    
    type(tBranch), Intent(in) :: tree
    
    Character(len=pd_mcl), Dimension(:), Allocatable, Intent(Out) :: desc
    Integer(kind=pd_ik)  , Dimension(:), Allocatable, Intent(Out) :: dat_no
    Integer(kind=1)      , Dimension(:), Allocatable, Intent(Out) :: dat_ty
    Integer(kind=pd_ik)  , Dimension(:), Allocatable, Intent(Out) :: offsets

    Integer(kind=pd_ik)              :: c
    Integer(kind=pd_ik)              :: alloc_stat
    Integer(kind=pd_ik)              :: ii

    c=0

    Allocate(desc(tree%streams%no_leaves), stat=alloc_stat)
    Call alloc_error(alloc_stat, "desc", &
                     "serialize_tree", tree%streams%no_leaves)
    desc = ""

    Allocate(dat_no(tree%streams%no_leaves), stat=alloc_stat)
    Call alloc_error(alloc_stat, "dat_no", &
                     "serialize_tree", tree%streams%no_leaves)
    dat_no = 0

    Allocate(dat_ty(tree%streams%no_leaves), stat=alloc_stat)
    Call alloc_error(alloc_stat, "dat_ty", &
                     "serialize_tree", tree%streams%no_leaves)
    dat_ty = 0

    Allocate(offsets(tree%streams%no_leaves), stat=alloc_stat)
    Call alloc_error(alloc_stat, "offsets", &
                     "serialize_tree", tree%streams%no_leaves)
    offsets = 0

    Do ii = 1, tree%no_leaves

       c          = c + 1
       desc(c)    = trim(tree%desc)//"#"//trim(tree%leaves(ii)%desc)
       dat_no(c)  = tree%leaves(ii)%dat_no
       dat_ty(c)  = tree%leaves(ii)%dat_ty
       offsets(c) = tree%leaves(ii)%lbound

    End Do

    Do ii = 1, tree%no_branches

       call serialize_tree_rec(tree%branches(ii), &
                               desc, dat_no, dat_ty, offsets ,c)

    End Do
   
  End Subroutine serialize_tree

  !============================================================================
  !> Subroutine for serializing a puredat tree to arrays
  !>
  !> Subroutine which recursively  serializes a puredat tree to arrays 
  !> regardless of the branch structures.
  !> Rule: Leaves first !!
  Recursive Subroutine serialize_tree_rec(tree, desc, dat_no, dat_ty, offsets, c)
    
    type(tBranch), Intent(in) :: tree
    
    Character(len=pd_mcl), Dimension(:), Intent(InOut) :: desc
    Integer(kind=pd_ik)  , Dimension(:), Intent(InOut) :: dat_no
    Integer(kind=1)      , Dimension(:), Intent(InOut) :: dat_ty
    Integer(kind=pd_ik)  , Dimension(:), Intent(InOut) :: offsets

    Integer(kind=pd_ik)                , Intent(inout) :: c

    Integer(kind=pd_ik)                                :: ii

    Do ii = 1, tree%no_leaves

       c          = c + 1  
       desc(c)    = trim(tree%desc)//"#"//trim(tree%leaves(ii)%desc)
       dat_no(c)  = tree%leaves(ii)%dat_no
       dat_ty(c)  = tree%leaves(ii)%dat_ty
       offsets(c) = tree%leaves(ii)%lbound

    End Do

    Do ii = 1, tree%no_branches

       call serialize_tree_rec(tree%branches(ii), &
                               desc, dat_no, dat_ty, offsets, c)

    End Do
   
  End Subroutine serialize_tree_rec

  !============================================================================
  !> Subroutine which serializes a complete tBranch structure
  Subroutine serialize_branch(branch,head,size,sdat)

    Type(tBranch)      , Intent(In)                               :: branch
    Integer(kind=pd_ik), Intent(inout), Dimension(:), Allocatable :: head
    Integer(kind=pd_ik), Intent(Out)                              :: size
    Logical, optional                                             :: sdat

    Integer(kind=pd_ik)              :: alloc_stat
    Logical                          :: loc_sdat
    
    !**********************************************************

    If ( present(sdat) ) then
       loc_sdat = sdat
    Else
       loc_sdat = .FALSE.
    End If 
    
    size = 0

    If (loc_sdat) then
       Call get_serial_branch_size_with_data(branch,size)
    Else
       Call get_serial_branch_size(branch,size)
    End If

    !Write(pd_umon,*)"Determined ",size," elements to serialize header"

    Allocate(head(size), stat=alloc_stat)
    Call alloc_error(alloc_stat, "head", "serialize_branch", size) 

    head = 0
    size = 1

    If (loc_sdat) then
       call serialize_branch_with_data_rec(branch,head,size)
    Else
       call serialize_branch_rec(branch,head,size)
    End If

    size = size - 1

!!$    un_head = pd_give_new_unit()
!!$    Open(unit=un_head, file=Trim(pro_path)//Trim(pro_name)//'.serialhead', status='replace', &
!!$         action='write',access="stream")
!!$    write(un_head)head
!!$    close(un_head)

  End Subroutine serialize_branch

  !============================================================================
  !> Subroutine which returns the size of a field of integer(kind=pd_ik)
  !> necessary to serializes a complete tBranch structure
  Recursive Subroutine get_serial_branch_size(branch, size)
    
    Type(tBranch)      , Intent(In)    :: branch
    Integer(kind=pd_ik), Intent(InOut) :: size

    Integer(kind=pd_ik)                :: ii

    !** Account for fixed components ***
    size = size + pd_ce + 2 

    !** Account for streams component if allocated ***
    If (allocated(branch%streams)) then
       size = size + 1 + 4*no_streams + no_streams*pd_ce
    Else
       size = size + 1
    End If

    !** Account for leaves if any *************
    size = size + branch%no_leaves * (no_streams + 4 + pd_ce)

    !** Account for branches if any ***
    Do ii = 1, branch%no_branches
       Call get_serial_branch_size(branch%branches(ii), size)
    End Do
       
  End Subroutine get_serial_branch_size
  
  !============================================================================
  !> Subroutine which returns the size of a field of integer(kind=pd_ik)
  !> necessary to serializes a complete tBranch structure with data.
  !>
  !> Subroutine which returns the size of a field of integer(kind=pd_ik)
  !> necessary to serializes a complete tBranch structure with data.<br>
  !> The routine only accounts for leafs with pstat >= 0 . This means only leafs
  !> with directly allocated data or data residing in serial streams are
  !> serialized. Leaves with parallel data marked with pstat = -1 are not
  !> regarded !
  Recursive Subroutine get_serial_branch_size_with_data(branch, size)
    
    Type(tBranch)      , Intent(In)    :: branch
    Integer(kind=pd_ik), Intent(InOut) :: size

    Integer(kind=pd_ik)                :: ii

    !** Account for fixed components ***
    size = size + pd_ce + 2 

    !** Account for streams component if allocated ***
    If (allocated(branch%streams)) then
       size = size + 1 + 4*no_streams + no_streams*pd_ce
    Else
       size = size + 1
    End If

    !** Account for leaves if any ****************
    Do ii = 1, branch%no_leaves

       size = size + (no_streams + 4 + pd_ce)

       !** Account for leaf data *****************
       if ( (branch%leaves(ii)%dat_no > 0) .AND. &
            (branch%leaves(ii)%pstat >= 0)        ) then
          
          Select Case (branch%leaves(ii)%dat_ty)
          
          Case (1)
             size = size + Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
          Case (2)
             size = size + Int(branch%leaves(ii)%dat_no/4,pd_ik)+1_pd_ik
          Case (3)
             size = size + Int(branch%leaves(ii)%dat_no/2,pd_ik)+1_pd_ik
          Case (4)
             size = size +     branch%leaves(ii)%dat_no
          Case (5)
             size = size +     branch%leaves(ii)%dat_no
          Case (6)
             size = size + Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
          Case (7)
             size = size + Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
          Case default
             Write(pd_umon,*)"Serialization of data type ",&
                  branch%leaves(ii)%dat_ty," is not jet implemented"
          End Select

       End if
       
    End Do
    
    !** Account for branches if any ***
    Do ii = 1, branch%no_branches
       Call get_serial_branch_size_with_data(branch%branches(ii), size)
    End Do
       
  End Subroutine get_serial_branch_size_with_data

  !============================================================================
  !> Subroutine which serializes a tBranch structure recursively
  Recursive Subroutine serialize_branch_rec(branch,head,pos)

    Type(tBranch)                    , Intent(In)    :: branch
    Integer(kind=pd_ik), Dimension(:), Intent(out)   :: head
    Integer(kind=pd_ik)              , Intent(InOut) :: pos

    Integer(kind=pd_ik),Dimension(pd_ce)             :: char_mold
    Integer(kind=pd_ik),Dimension(pd_ce*no_streams)  :: char_mold7

    Integer(kind=pd_ik)                              :: ii

    !** Fixed Components ******************************************************
    head(pos:pos+pd_ce-1) = Transfer(branch%desc,char_mold)
    pos = pos+pd_ce

    head(pos) = branch%no_branches
    pos = pos+1
    head(pos) = branch%no_leaves
    pos = pos+1

    !** Streams ***************************************************************
    If (allocated(branch%streams)) then
       
       head(pos) = 1
       pos = pos+1

       head(pos:pos+no_streams-1) = branch%streams%dim_st
       pos = pos+no_streams
       
       head(pos:pos+no_streams-1) = branch%streams%ii_st
       pos = pos+no_streams

       head(pos:pos+no_streams*pd_ce-1) = Transfer(branch%streams%stream_files,char_mold7)
       pos = pos+no_streams*pd_ce
      
       Do ii = 1, no_streams
          if (branch%streams%ifopen(ii)) then
             head(pos) = 1
          Else 
             head(pos) = 0
          End if
          pos = pos+1
       End Do

       head(pos:pos+no_streams-1) = branch%streams%units
       pos = pos+no_streams

    Else

       head(pos) = 0
       pos = pos+1

    End If

    !** Leaves ****************************************************************
    Do ii = 1, branch%no_leaves
       
       head(pos:pos+pd_ce-1) = Transfer(branch%leaves(ii)%desc,char_mold)
       pos = pos+pd_ce

       head(pos) = branch%leaves(ii)%dat_no
       pos = pos+1
       head(pos) = branch%leaves(ii)%dat_ty
       pos = pos+1
       head(pos) = branch%leaves(ii)%lbound
       pos = pos+1
       head(pos) = branch%leaves(ii)%ubound
       pos = pos+1

       head(pos) = branch%leaves(ii)%pstat
       pos = pos + 1

    End Do

    !** Branches **************************************************************
    Do ii = 1, branch%no_branches
       Call serialize_branch_rec(branch%branches(ii),head,pos)
    End Do

  End Subroutine serialize_branch_rec

  !============================================================================
  !> Subroutine which serializes a tBranch structure recursively
  !>
  !> Subroutine which serializes a tBranch structure recursively<br>
  !> The routine only accounts for leafs with pstat >= 0 . This means only leafs
  !> with directly allocated data or data residing in serial streams are
  !> serialized. Leaves with parallel ata marked with pstat = -1 are not
  !> regarded !
  Recursive Subroutine serialize_branch_with_data_rec(branch,head,pos)

    Type(tBranch)                    , Intent(In)    :: branch
    Integer(kind=pd_ik), Dimension(:), Intent(out)   :: head
    Integer(kind=pd_ik)              , Intent(InOut) :: pos

    Integer(kind=pd_ik),Dimension(pd_ce)             :: char_mold
    Integer(kind=pd_ik),Dimension(pd_ce*no_streams)  :: char_mold7

    Integer(kind=pd_ik)                              :: ii, no_int8_elems

    !** Fixed Components ******************************************************
    head(pos:pos+pd_ce-1) = Transfer(branch%desc,char_mold)
    pos = pos+pd_ce

    head(pos) = branch%no_branches
    pos = pos+1
    head(pos) = branch%no_leaves
    pos = pos+1

    !** Streams ***************************************************************
    If (allocated(branch%streams)) then
       
       head(pos) = 1
       pos = pos+1

       head(pos:pos+no_streams-1) = branch%streams%dim_st
       pos = pos+no_streams
       
       head(pos:pos+no_streams-1) = branch%streams%ii_st
       pos = pos+no_streams

       head(pos:pos+no_streams*pd_ce-1) = Transfer(branch%streams%stream_files,char_mold7)
       pos = pos+no_streams*pd_ce
      
       Do ii = 1, no_streams
          if (branch%streams%ifopen(ii)) then
             head(pos) = 1
          Else 
             head(pos) = 0
          End if
          pos = pos+1
       End Do

       head(pos:pos+no_streams-1) = branch%streams%units
       pos = pos+no_streams

    Else

       head(pos) = 0
       pos = pos+1

    End If

    !** Leaves ****************************************************************
    Do ii = 1, branch%no_leaves
       
       head(pos:pos+pd_ce-1) = Transfer(branch%leaves(ii)%desc,char_mold)
       pos = pos+pd_ce

       head(pos) = branch%leaves(ii)%dat_no
       pos = pos+1
       head(pos) = branch%leaves(ii)%dat_ty
       pos = pos+1
       head(pos) = branch%leaves(ii)%lbound
       pos = pos+1
       head(pos) = branch%leaves(ii)%ubound
       pos = pos+1

       head(pos) = branch%leaves(ii)%pstat
       pos = pos + 1

       !** Serialize leaf data *****************
       if ( (branch%leaves(ii)%dat_no > 0) .AND. &
            (branch%leaves(ii)%pstat >= 0)        ) then
          
          Select Case (branch%leaves(ii)%dat_ty)
          
          Case (1)

             no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_int1,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case (2)
             
             no_int8_elems = Int(branch%leaves(ii)%dat_no/4,pd_ik)+1_pd_ik
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_int2,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case (3)
             
             no_int8_elems = Int(branch%leaves(ii)%dat_no/2,pd_ik)+1_pd_ik
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_int4,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case (4)
             
             no_int8_elems =     branch%leaves(ii)%dat_no
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_int8,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case (5)
             
             no_int8_elems =     branch%leaves(ii)%dat_no
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_real8,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case (6)
             
             no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_char,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems

          Case (7)
             
             no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
             head(pos:pos+no_int8_elems-1) = Transfer(branch%leaves(ii)%p_log,head(pos:pos+no_int8_elems-1))
             pos = pos + no_int8_elems
             
          Case default
             Write(pd_umon,*)"Serialization of data type ",branch%leaves(ii)%dat_ty," is not jet implemented"
          End Select

       End if
       
    End Do

    !** Branches **************************************************************
    Do ii = 1, branch%no_branches
       Call serialize_branch_with_data_rec(branch%branches(ii),head,pos)
    End Do

  End Subroutine serialize_branch_with_data_rec

  !============================================================================
  !> Subroutine which deserializes a tBranch structure
  Subroutine deserialize_branch(branch, head, sdat)

    Type(tBranch)                    , Intent(Out) :: branch
    Integer(kind=pd_ik), Dimension(:), Intent(in)  :: head
    Logical, optional                              :: sdat

    Logical                                        :: loc_sdat
    Integer(kind=pd_ik)                            :: pos

    If ( present(sdat) ) then
       loc_sdat = sdat
    Else
       loc_sdat = .FALSE.
    End If
    
    pos = 1

    If (loc_sdat) then
       
       call deserialize_branch_with_data_rec(branch,head,pos)
       
    Else
       
       call deserialize_branch_rec(branch,head,pos)

    End If

  End Subroutine deserialize_branch


  !============================================================================
  !> Subroutine which deserializes a tBranch structure recursively
  Recursive Subroutine deserialize_branch_rec(branch,head,pos)!,no_l,no_b)

    Type(tBranch)                    , Intent(InOut) :: branch
    Integer(kind=pd_ik), Dimension(:), Intent(in)    :: head
    Integer(kind=pd_ik)              , Intent(InOut) :: pos!, no_l,no_b

    CHARACTER(len=pd_mcl)                            :: char_mold
    Integer(kind=pd_ik)                              :: ii

    !** Fixed Components ******************************************************
    branch%desc = Transfer(head(pos:pos+pd_ce-1),branch%desc)
    pos = pos+pd_ce

    branch%no_branches = head(pos)
    pos = pos+1
    branch%no_leaves = head(pos)
    pos = pos+1

    !** Streams ***************************************************************
    If (head(pos) == 1) then

       pos = pos+1       
       allocate(branch%streams)

       branch%streams%dim_st = head(pos:pos+no_streams-1)
       pos = pos+no_streams
       
       branch%streams%ii_st = head(pos:pos+no_streams-1)
       pos = pos+no_streams

       branch%streams%stream_files = ""
       branch%streams%stream_files = Transfer(head(pos:pos+no_streams*pd_ce-1),&
            char_mold,no_streams)
       pos = pos+no_streams*pd_ce
      
       Do ii = 1, no_streams
          if (head(pos) == 1) then
             branch%streams%ifopen(ii) = .TRUE.
          Else 
             branch%streams%ifopen(ii) = .FALSE.
          End if
          pos = pos+1
       End Do

       branch%streams%units = head(pos:pos+no_streams-1)
       pos = pos+no_streams

    Else

       pos = pos+1

    End If

    !** Leaves ****************************************************************
    If (branch%no_leaves > 0) then

       Allocate(branch%leaves(branch%no_leaves))
       !no_l = no_l + branch%no_leaves
       !write(pd_umon,'(A,I0)')"Allocated leaves   : ",branch%no_leaves
       Do ii = 1, branch%no_leaves
       
          branch%leaves(ii)%desc = Transfer(head(pos:pos+pd_ce-1),branch%leaves(ii)%desc)
          pos = pos+pd_ce
          
          branch%leaves(ii)%dat_no = head(pos)
          pos = pos+1
          branch%leaves(ii)%dat_ty = Int(head(pos),1)
          pos = pos+1
          branch%leaves(ii)%lbound = head(pos)
          pos = pos+1
          branch%leaves(ii)%ubound = head(pos)
          pos = pos+1
          
          branch%leaves(ii)%pstat = head(pos)
          pos = pos + 1
             
       End Do

    End If

    !** Branches **************************************************************
    If ( branch%no_branches > 0 ) then
       
       Allocate(branch%branches(branch%no_branches))
       !no_b = no_b + branch%no_branches
       !write(pd_umon,'(A,I0)')"Allocated branches : ",branch%no_branches
       Do ii = 1, branch%no_branches
          Call deserialize_branch_rec(branch%branches(ii),head,pos)!,no_l,no_b)
       End Do

    End If

  End Subroutine deserialize_branch_rec

  !============================================================================
  !> Subroutine which deserializes a tBranch structure recursively
  !>
  !> !! The deserialization is done directly to the leaf pointers
  !> \TODO Implement deserialisation of the leaf data directly to stream arrays
  Recursive Subroutine deserialize_branch_with_data_rec(branch,head,pos)!,no_l,no_b)

    Type(tBranch)                    , Intent(InOut) :: branch
    Integer(kind=pd_ik), Dimension(:), Intent(in)    :: head
    Integer(kind=pd_ik)              , Intent(InOut) :: pos!, no_l,no_b

    CHARACTER(len=pd_mcl)                            :: char_mold
    Integer(kind=pd_ik)                              :: ii, no_int8_elems

    !** Fixed Components ******************************************************
    branch%desc = Transfer(head(pos:pos+pd_ce-1),branch%desc)
    pos = pos+pd_ce

    branch%no_branches = head(pos)
    pos = pos+1
    branch%no_leaves = head(pos)
    pos = pos+1

    !** Streams ***************************************************************
    If (head(pos) == 1) then

       pos = pos+1       
       allocate(branch%streams)

       branch%streams%dim_st = head(pos:pos+no_streams-1)
       pos = pos+no_streams
       
       branch%streams%ii_st = head(pos:pos+no_streams-1)
       pos = pos+no_streams

       branch%streams%stream_files = ""
       branch%streams%stream_files = Transfer(head(pos:pos+no_streams*pd_ce-1),&
            char_mold,no_streams)
       pos = pos+no_streams*pd_ce
      
       Do ii = 1, no_streams
          if (head(pos) == 1) then
             branch%streams%ifopen(ii) = .TRUE.
          Else 
             branch%streams%ifopen(ii) = .FALSE.
          End if
          pos = pos+1
       End Do

       branch%streams%units = head(pos:pos+no_streams-1)
       pos = pos+no_streams

    Else

       pos = pos+1

    End If

    !** Leaves ****************************************************************
    If (branch%no_leaves > 0) then

       Allocate(branch%leaves(branch%no_leaves))
       !no_l = no_l + branch%no_leaves
       !write(pd_umon,'(A,I0)')"Allocated leaves   : ",branch%no_leaves
       Do ii = 1, branch%no_leaves
       
          branch%leaves(ii)%desc = Transfer(head(pos:pos+pd_ce-1),branch%leaves(ii)%desc)
          pos = pos+pd_ce
          
          branch%leaves(ii)%dat_no = head(pos)
          pos = pos+1
          branch%leaves(ii)%dat_ty = Int(head(pos),1)
          pos = pos+1
          branch%leaves(ii)%lbound = head(pos)
          pos = pos+1
          branch%leaves(ii)%ubound = head(pos)
          pos = pos+1
          
          branch%leaves(ii)%pstat = head(pos)
          pos = pos + 1

          !** DeSerialize leaf data *****************
          if ( (branch%leaves(ii)%dat_no > 0) .AND. &
               (branch%leaves(ii)%pstat >= 0) ) then
                    
             Select Case (branch%leaves(ii)%dat_ty)
          
             Case (1)

                Allocate(branch%leaves(ii)%p_int1(branch%leaves(ii)%dat_no))
                no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
                branch%leaves(ii)%p_int1 = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_int1)
                pos = pos + no_int8_elems
                
             Case (2)

                Allocate(branch%leaves(ii)%p_int2(branch%leaves(ii)%dat_no))
                no_int8_elems = Int(branch%leaves(ii)%dat_no/4,pd_ik)+1_pd_ik
                branch%leaves(ii)%p_int2 = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_int2)
                pos = pos + no_int8_elems
                
             Case (3)

                Allocate(branch%leaves(ii)%p_int4(branch%leaves(ii)%dat_no))
                no_int8_elems = Int(branch%leaves(ii)%dat_no/2,pd_ik)+1_pd_ik
                branch%leaves(ii)%p_int4 = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_int4)
                pos = pos + no_int8_elems
                
             Case (4)

                Allocate(branch%leaves(ii)%p_int8(branch%leaves(ii)%dat_no))
                no_int8_elems =     branch%leaves(ii)%dat_no
                branch%leaves(ii)%p_int8 = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_int8)
                pos = pos + no_int8_elems
                
             Case (5)

                Allocate(branch%leaves(ii)%p_real8(branch%leaves(ii)%dat_no))
                no_int8_elems =     branch%leaves(ii)%dat_no
                branch%leaves(ii)%p_real8 = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_real8)
                pos = pos + no_int8_elems
                
             Case (6)

                Allocate(branch%leaves(ii)%p_char(branch%leaves(ii)%dat_no))
                no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
                branch%leaves(ii)%p_char = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_char)
                pos = pos + no_int8_elems

             Case (7)

                Allocate(branch%leaves(ii)%p_log(branch%leaves(ii)%dat_no))
                no_int8_elems = Int(branch%leaves(ii)%dat_no/8,pd_ik)+1_pd_ik
                branch%leaves(ii)%p_log = Transfer(head(pos:pos+no_int8_elems-1),branch%leaves(ii)%p_log)
                pos = pos + no_int8_elems
                
             Case default
                Write(pd_umon,*)"DeSerialization of data type ",branch%leaves(ii)%dat_ty," is not jet implemented"
             End Select

          End if
             
       End Do

    End If

    !** Branches **************************************************************
    If ( branch%no_branches > 0 ) then
       
       Allocate(branch%branches(branch%no_branches))
       !no_b = no_b + branch%no_branches
       !write(pd_umon,'(A,I0)')"Allocated branches : ",branch%no_branches
       Do ii = 1, branch%no_branches
          Call deserialize_branch_with_data_rec(branch%branches(ii),head,pos)!,no_l,no_b)
       End Do

    End If

  End Subroutine deserialize_branch_with_data_rec
  !> @} 
  !# End of memeber group "Puredat  serialization routines" ###################


  !############################################################################
  !> \name Puredat search routines
  !> @{
  !> Module procedures for searching elements in a puredat tree structure

  !============================================================================
  !> Subroutine to search a branch in a tTree structure with success warning
  !>
  !> The subroutine calls search_branch rec and warns to std out in case a
  !> branch withdescr was not found in the tBranch structure.
  Subroutine Search_branch_wrn(descr, branch, out_branch, success, wrn)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In), Target :: branch
    Type(tBranch), Pointer , Intent(out)        :: out_branch
    Logical                , Intent(inout)      :: success
    Logical                , Intent(in)         :: wrn

    Integer                                     :: ii

    If (Trim(descr) == Trim(branch%desc)) Then

       out_branch => branch
       success    = .True.

    Else 

       success = .FALSE.

       Do ii = 1, branch%no_branches

          Call Search_branch(descr, branch%branches(ii), out_branch, success)
          if (success) exit

       End Do
       
    End If

    If (wrn .AND. (.NOT.success)) then
       write(*,PDF_W_A)"Branch with descr"
       write(*,PDF_W_A)trim(descr)
       write(*,PDF_W_A)"was not found in branch with descr"
       write(*,PDF_W_A)trim(branch%desc)
    End If
    
  End Subroutine Search_branch_wrn
  
  !============================================================================
  !> Subroutine to search a branch in a tTree structure
  !>
  !> The subroutine recursively cycles the tBranch members of a tTree 
  !> structure searching the specified branch by full comarison of 
  !> tBranch::desc and returning the first element found
  Recursive Subroutine Search_branch_rec(descr, branch, out_branch, success)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In), Target :: branch
    Type(tBranch), Pointer , Intent(out)        :: out_branch
    Logical                , Intent(inout)      :: success

    Integer                                     :: ii

    If (Trim(descr) == Trim(branch%desc)) Then

       out_branch => branch
       success    = .True.

    Else 

       success = .FALSE.

       Do ii = 1, branch%no_branches

          Call Search_branch(descr, branch%branches(ii), out_branch, success)
          if (success) exit

       End Do
       
    End If

  End Subroutine Search_branch_rec

  !============================================================================
  !> Subroutine to count the multiplicity of branches in a tBranch structure
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_branch_num_rec
  Subroutine get_branch_num(descr, branch, num)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In)         :: branch

    Integer(kind=pd_ik)    , Intent(inout)      :: num

    num = 0

    call get_branch_num_rec(descr, branch, num)

  End Subroutine get_branch_num

  !============================================================================
  !> Subroutine to count the multiplicity of branches in a tBranch structure
  !>
  !> The function recursively cycles the members of a tBranch  
  !> structure searching the specified branch by full comarison of 
  !> tbranch::desc and counting its multiplicity
  Recursive Subroutine get_branch_num_rec(descr, branch, num)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In)         :: branch

    Integer(kind=pd_ik)    , Intent(inout)      :: num

    Integer                                     :: ii

    If (Trim(descr) == Trim(branch%desc)) Then
       num = num + 1
    End If

    Branch_Cycle : Do ii = 1, branch%no_branches   
       call get_branch_num_rec(descr, branch%branches(ii),num)
    End Do Branch_Cycle

  End Subroutine get_branch_num_rec

  !============================================================================
  !> Subroutine to return a branch specified by its number
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_branch_with_num_rec
  Recursive Subroutine get_branch_with_num(descr, branch, num, out_branch, success)
    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In), Target :: branch
    Type(tBranch), Pointer , Intent(out)        :: out_branch
    Integer(kind=pd_ik)    , Intent(in)         :: num
    Logical                , Intent(inout)      :: success

    Integer(kind=pd_ik)                         :: count    

    count = 0
    call get_branch_with_num_rec(descr, branch, num, out_branch, count, success)

  End Subroutine get_branch_with_num

  !============================================================================
  !> Subroutine to return a branch specified by its number 
  !>
  !> The function recursively cycles the members of a tBranch counting each
  !> member and returning the branch with the given number and full comarison of 
  !> tbranch::desc
  Recursive Subroutine get_branch_with_num_rec(descr, branch, num, out_branch, count, success)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In), Target :: branch
    Type(tBranch), Pointer , Intent(out)        :: out_branch
    Integer(kind=pd_ik)    , Intent(in)         :: num
    Integer(kind=pd_ik)    , Intent(inout)      :: count
    Logical                , Intent(inout)      :: success
    Integer                                     :: ii

    success = .FALSE.
    count = count + 1

    if (count /= num) then
       
       count = count + branch%no_leaves

       Branch_Cycle : Do ii = 1, branch%no_branches
          call get_branch_with_num_rec(descr, branch%branches(ii), num, out_branch, count, success)
          if (success) exit
       End Do Branch_Cycle

    Else

       If ( trim(descr) /= trim(branch%desc) ) then

          WRITE(pd_umon,*)
          WRITE(pd_umon,PDF_SEP)
          WRITE(pd_umon,PDF_E_A)   'An error occured in routine : get_branch_with_num_rec'
          WRITE(pd_umon,PDF_E_AI0) 'The branch with number      : ',num
          WRITE(pd_umon,PDF_E_A)   'has description :'
          WRITE(pd_umon,PDF_E_A)   trim(branch%desc)
          WRITE(pd_umon,PDF_E_A)   'but it was searched for description :'
          WRITE(pd_umon,PDF_E_A)   descr
          WRITE(pd_umon,PDF_E_STOP)
          STOP

       end If

       success = .TRUE.
       out_branch => branch

    End if  
    
  End Subroutine get_branch_with_num_rec

  !============================================================================
  !> Subroutine to count the multiplicity of leaf in a tBranch structure
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_leaf_num_rec
  Subroutine get_leaf_num(descr, branch, num)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In)         :: branch

    Integer(kind=pd_ik)    , Intent(inout)      :: num

    num = 0

    call get_leaf_num_rec(descr, branch, num)

  End Subroutine get_leaf_num

  !============================================================================
  !> Subroutine to count the multiplicity of leaf in a tBranch structure
  !>
  !> The function recursively cycles the members of a tBranch  
  !> structure searching the specified leaf by full comarison of 
  !> tleaf::desc and counting its multiplicity
  Recursive Subroutine get_leaf_num_rec(descr, branch, num)

    Character(len=*)       , Intent(In)         :: descr
    Type(tBranch)          , Intent(In)         :: branch

    Integer(kind=pd_ik)    , Intent(inout)      :: num

    Integer                                     :: ii

    Leaf_cycle : Do ii = 1, branch%no_leaves

       If (Trim(descr) == Trim(branch%leaves(ii)%desc)) Then
          num = num + 1
       End If

    End Do Leaf_cycle

    Branch_Cycle : Do ii = 1, branch%no_branches
       call get_leaf_num_rec(descr, branch%branches(ii),num)
    End Do Branch_Cycle

  End Subroutine get_leaf_num_rec

  !============================================================================
  !> Subroutine to get a list of leafs out of a tBranch structure
  !>
  !> The function recursively cycles the members of a tBranch  
  !> structure searching the specified leafs by full comarison of 
  !> tleaf::desc and retrieving them to a list. This Routine is a wrapper 
  !> for the recursive subroutine which has to carry the aditional parameter 
  !> leaf_count
  Subroutine get_leaf_list(descr, branch, leaf_num, leaf_list)
    
    Character(len=*)       , Intent(In)                 :: descr
    Type(tBranch)          , Intent(In)                 :: branch
    Integer(kind=pd_ik)    , Intent(out)                :: leaf_num
        
    Type(tLeaf), Dimension(:), Allocatable, intent(out) :: leaf_list
    
    Integer(kind=pd_ik)                                 :: leaf_count

    integer                                             :: alloc_stat

    leaf_num = 0

    call  get_leaf_num(descr, branch, leaf_num)

    leaf_count = 0

    Allocate(leaf_list(leaf_num), stat=alloc_stat)
    call alloc_error(alloc_stat, "leaf_list", "get_leaf_list", leaf_num)

    call  get_leaf_list_rec(descr, branch, leaf_num, leaf_count, leaf_list)

  End Subroutine get_leaf_list

  !============================================================================
  !> Subroutine to get a list of leafs out of a tBranch structure
  !>
  !> The function recursively cycles the members of a tBranch  
  !> structure searching the specified leafs by full comarison of 
  !> tleaf::desc and retrieving them to a list. To be used in conjunction with
  !> get_leaf_num
  Recursive Subroutine get_leaf_list_rec(descr, branch, leaf_num, &
                                         leaf_count, leaf_list)

    Character(len=*)       , Intent(In)             :: descr
    Type(tBranch)          , Intent(In)             :: branch
    Integer(kind=pd_ik)    , Intent(in)             :: leaf_num
    Integer(kind=pd_ik)    , Intent(inOut)          :: leaf_count

    Type(tLeaf), Dimension(leaf_num), intent(inout) :: leaf_list

    Integer                                         :: ii

    Leaf_cycle : Do ii = 1, branch%no_leaves

       If (Trim(descr) == Trim(branch%leaves(ii)%desc)) Then
          leaf_count = leaf_count + 1
          leaf_list(leaf_count) = branch%leaves(ii)
       End If

    End Do Leaf_cycle

    Branch_Cycle : Do ii = 1, branch%no_branches
       call get_leaf_list_rec(descr, branch%branches(ii),leaf_num, &
            leaf_count, leaf_list)
    End Do Branch_Cycle

  End Subroutine get_leaf_list_rec

  !============================================================================
  !> Subroutine to get a tLeaf Structure by number from a tBranch structure
  !>
  !> The subroutine is a wrapper for the recursive subroutine 
  !> get_leaf_with_num_rec
  Subroutine get_leaf_with_num(tree, num, out_leaf, success)

    Type(tBranch)          , Intent(In), Target :: tree
    Type(tLeaf),   Pointer , Intent(out)        :: out_leaf
    Integer(kind=pd_ik)    , Intent(in)         :: num
    Logical                , Intent(inout)      :: success

    Integer(kind=pd_ik)                         :: count    

    count = 0
    call get_leaf_with_num_rec(tree, num, out_leaf, count, success)

  End Subroutine get_leaf_with_num

  !============================================================================
  !> Subroutine to get a tLeaf Structure by number from a tBranch structure
  !>
  !> The function recursively cycles the members of a tBranch counting each
  !> member and returning the leaf with the given number and full comarison of 
  !> tbranch::desc
  Recursive Subroutine get_leaf_with_num_rec(tree, num, out_leaf, count, success)

    Type(tBranch)          , Intent(In), Target :: tree
    Type(tLeaf),   Pointer , Intent(out)        :: out_leaf
    Integer(kind=pd_ik)    , Intent(in)         :: num
    Integer(kind=pd_ik)    , Intent(inout)      :: count
    Logical                , Intent(inout)      :: success
    Integer                                     :: ii

    success = .FALSE.
    count = count + 1

    if (count /= num) then
       
       Do ii = 1, tree%no_leaves

          count = count + 1

          if (count == num) then

             success = .TRUE.
             out_leaf => tree%leaves(ii)
             exit
          End if

       End Do

       if  (count /= num) then
          Branch_Cycle : Do ii = 1, tree%no_branches
             call get_leaf_with_num_rec(tree%branches(ii), num, out_leaf, count, success)
             if (success) exit
          End Do Branch_Cycle
       End if

    Else

       WRITE(pd_umon,*)
       WRITE(pd_umon,PDF_SEP)
       WRITE(pd_umon,PDF_E_A)   'An error occured in routine : get_branch_with_num_rec'
       WRITE(pd_umon,PDF_E_AI0) 'The element with number     : ',num
       WRITE(pd_umon,PDF_E_A)   'was found but is a branch with description :'
       WRITE(pd_umon,PDF_E_A)   trim(tree%desc)
       WRITE(pd_umon,PDF_E_A)   'To get a tBranch structure by number please use the'
       WRITE(pd_umon,PDF_E_A)   'PureDat subroutine get_branch_with_num.'
       WRITE(pd_umon,PDF_E_STOP)
       STOP

    end If
    
  End Subroutine get_leaf_with_num_rec

  !> @} 
  !# End of memeber group "Puredat search routines" ###########################

  !############################################################################
  !> \name Puredat logging routines
  !> @{
  !> These routines are intended for debug output of puredat data structures  
  !> NOT for production output !!             

  !============================================================================
  !> Subroutine for logging a complete puredat tree with all content
  Subroutine log_tree(tree, unit_lf, data, commands)

    Type(tBranch)   , Intent(In)                 :: tree
    Integer         , Intent(in)                 :: unit_lf
    Logical         , Intent(In), Optional       :: data, commands

    Integer                                      :: spacer

    Character(Len=pd_mcl)                        :: sep_dash  = ''
    Logical                                      :: loc_data, loc_commands

    If (present(data)) Then
       loc_data = data
    Else
       loc_data = .FALSE.
    End If
    If (present(commands)) Then
       loc_commands = commands
    Else
       loc_commands = .FALSE.
    End If

    spacer=Len_Trim(tree%desc)+23
    Write(sep_dash ,'(A,I0,A)')"('+-',",spacer-3,"('-'),'-+')"

    Write(unit_lf,*)
    Write(unit_lf, sep_dash)
    Write(unit_lf, '(A,A,A)')"| Structure of tree '",Trim(tree%desc),"' |"
    Write(unit_lf, sep_dash)

    Call log_branch(tree ,unit_lf, " ", loc_data, loc_commands)

  End Subroutine log_tree
  
  !============================================================================
  !> Subroutine for logging a tBranch structure 
  Recursive Subroutine log_branch(branch, unit_lf, fmt_str_in, data, commands)

    Type(tBranch)   , Intent(In)            :: branch
    Integer         , Intent(in)            :: unit_lf
    Character(len=*), Intent(In) , optional :: fmt_str_in
    Logical         , Intent(In) , optional :: data, commands

    Integer(kind=pd_ik)                     :: ii, lt_desc, len_fmt_str
    Integer                                 :: spacer

    Character(Len=pd_mcl)                   :: b_sep

    Character(len=pd_mcl)                   :: fmt_str   = ''

    Logical                                 :: loc_data, loc_commands

    spacer=Len_Trim(branch%desc)+23

    if (present(fmt_str_in)) then
       fmt_str     = fmt_str_in
       len_fmt_str = Len(fmt_str_in)
    Else
       fmt_str     = ' '
       len_fmt_str = 1
    end if
    
    If (present(data)) Then
       loc_data = data
    Else
       loc_data = .FALSE.
    End If
    If (present(commands)) Then
       loc_commands = commands
    Else
       loc_commands = .FALSE.
    End If

    Write(unit_lf, "('"//fmt_str(1:len_fmt_str-1)//"|')")

    If ( Len_Trim(branch%desc) < 1 ) then
       lt_desc = 1 
    Else
       lt_desc = Len_Trim(branch%desc)
    end If

    Write(b_sep,'(A,I0,A)')"('"//fmt_str(1:len_fmt_str-1)//"|   +----------',",LT_desc,"('-'),'-+')"

    Write(unit_lf, b_sep)
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str-1)//"+---|',A,A,A)") ' Branch : ',Trim(branch%desc),' |'

    Write(b_sep,'(A,I0,A)')"('"//fmt_str(1:len_fmt_str)//"   +----------',",LT_desc,"('-'),'-+')"
    Write(unit_lf, b_sep)

    If (allocated(branch%streams)) then

       Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
       Do ii = 1, no_streams
          Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Dimension of stream        ',I0,' = ',I0)")&
               ii,branch%streams%dim_st(ii)
       End Do
       Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
 
       Do ii = 1, no_streams
          Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Position of stream counter ',I0,' = ',I0)")&
               ii, branch%streams%ii_st(ii)
       End Do
       Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
       Do ii = 1, no_streams
          Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Stream-file                ',I0,' = ',A)")&
               ii, trim(branch%streams%Stream_files(ii))
       End Do
       Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
       Do ii = 1, no_streams
          Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Units for stream files     ',I0,' = ',I0)")&
               ii, branch%streams%units(ii)
       End Do
       Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
       Do ii = 1, no_streams
          Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Open state of stream file  ',I0,' = ',L1)")&
               ii, branch%streams%ifopen(ii)
       End Do

    End If

    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
 
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- No of branches = ',I0)")branch%no_branches
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- No of leaves   = ',I0)")branch%no_leaves

    Do ii = 1, branch%no_branches

       If ((ii < branch%no_branches) .Or. (branch%no_leaves > 0)) Then
          Call log_branch(branch%branches(ii) ,unit_lf, fmt_str(1:len_fmt_str)//"   |", loc_data, loc_commands)
       Else
          Call log_branch(branch%branches(ii) ,unit_lf, fmt_str(1:len_fmt_str)//"    ", loc_data, loc_commands)
       End If

    End Do

    Do ii = 1, branch%no_leaves

       If (ii < branch%no_leaves) Then
          Call log_leaf(branch%leaves(ii), unit_lf, fmt_str(1:len_fmt_str)//"   |", data, commands, branch)
       Else
          Call log_leaf(branch%leaves(ii), unit_lf, fmt_str(1:len_fmt_str)//"    ", data, commands, branch) 
       End If
       
    End Do

  End Subroutine log_branch

  !============================================================================
  !> Subroutine for logging a tBranch structure
  Subroutine log_leaf(tree, unit_lf, fmt_str_in, data, commands, parent)

    Type(tLeaf)     , Intent(In)            :: tree
    Integer         , Intent(in)            :: unit_lf
    Character(len=*), Intent(In) , optional :: fmt_str_in
    Logical         , Intent(In) , optional :: data, commands
    Type(tBranch)   , Intent(in) , optional :: parent

    Character(Len=pd_mcl)           :: b_sep, fmt_str
    Integer                         :: lt_desc, len_fmt_str

    Logical                         :: loc_data, loc_commands

    if (present(fmt_str_in)) then
       fmt_str     = fmt_str_in
       len_fmt_str = Len(fmt_str_in)
    Else
       fmt_str     = ' '
       len_fmt_str = 1
    end if

    If (present(data)) Then
       loc_data = data
    Else
       loc_data = .FALSE.
    End If
    If (present(commands) .AND. present(parent) ) Then
       loc_commands = commands
    Else
       loc_commands = .FALSE.
    End If

    if(present(commands) .AND. (.NOT.present(parent))) then
       Write(unit_lf,*)"Optional parameter 'commands' is present &
            &without optional parameter 'parent'"
       Write(unit_lf,*)"This parameter combination is not implemented. &
            &'commands' is reseted to .FALSE."
       loc_commands = .FALSE.
    End if

    Write(unit_lf, "('"//fmt_str(1:Len_fmt_str-1)//"|')")

    If ( Len_Trim(tree%desc) < 1 ) then
       lt_desc = 1 
    Else
       lt_desc = Len_Trim(tree%desc)
    end If

    Write(b_sep,'(A,I0,A)')"('"//fmt_str(1:Len_fmt_str-1)//"|   +--------',",LT_desc,"('-'),'-+')"
    Write(unit_lf, b_sep)
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str-1)//"+---| ',A,A)")'Leaf : ',Trim(tree%desc)//' |'
    Write(b_sep,'(A,I0,A)')"('"//fmt_str(1:len_fmt_str)//"   +--------',",LT_desc,"('-'),'-+')"
    Write(unit_lf, b_sep)

    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   |')")
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- No of data           : ',I0)")tree%dat_no
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Type of data         : ',I0)")tree%dat_ty
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Lower bound in stream: ',I0)")tree%lbound
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- Upper bound in stream: ',I0)")tree%ubound
    Write(unit_lf, "('"//fmt_str(1:len_fmt_str)//"   +-- PStat                : ',I0)")tree%pstat

    if (loc_commands) then
       Write(unit_lf,"('"//fmt_str(1:len_fmt_str)//"   |')")
       Write(unit_lf,&
            "('"//fmt_str(1:len_fmt_str)//"   +-- Command for leaf extraction : ',3(A,1X),2(A,A,A,1X),A,1X,A)") &
            "pd_leaf_to_file_x86_64", &
            trim(pro_path),  trim(pro_name), &
            '"',trim(tree%desc),'"', '"',trim(parent%desc),'"', &
            trim(pro_path),  trim(pro_name)//"--leaf_data.dat"
    End if

    If (loc_data) then

       Call log_Leaf_data(unit_lf, tree, fmt_str(1:len_fmt_str))
      
    End If

  contains

    !**************************************************************************
    !> This subroutine prints out data contained in a tLeaf structure.
    !> It is meant to be called only by log_leaf
    Subroutine Log_Leaf_data(un_lf, leaf, fmt_str)

      Integer         , Intent(in) :: un_lf
      Type(tLeaf)     , Intent(In) :: leaf
      Character(Len=*), Intent(In) :: fmt_str

      Integer :: ii, pos, cpos

      !** We have to have data residing in the leaf and we have to have a
      !** leaf pointer status > 0 which marks process owned data residing
      !** either directly below the leaf pointer or in a tStreams structure
      !** pointed to by the leaf pointer
      If ( (leaf%dat_no > 0) .AND. leaf%pstat >= 0) then
         
         Write(un_lf, "('"//fmt_str//"   |')")
         Write(un_lf, "('"//fmt_str//"   +',170('-'),('+'))")
         Write(un_lf, "('"//fmt_str//"   |',' Leaf Data ',159(' '),('|'))")

         !***********************************************************
         !** dat_no in [0,30] ***************************************
         If (leaf%dat_no <= 30) then
            
            Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
            
            pos = 1
            Do ii = 1, leaf%dat_no
               
               Select Case (leaf%dat_ty)
                  
               Case (1)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int1(ii)
               Case (2)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int2(ii)
               Case (3)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int4(ii)
               Case (4)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int8(ii)
               Case (5)
                  Write(un_lf, "(E16.6,' ')", ADVANCE='NO')leaf%p_real8(ii)
               Case (6)
                  if (Ichar(leaf%p_char(ii)) /= 10) then
                     Write(un_lf, "(A1)", ADVANCE='NO')leaf%p_char(ii)
                  Else
                     Write(un_lf, "(A1)", ADVANCE='NO')'_'
                  End if
               Case (7)
                  Write(un_lf, "(L1,' ')", ADVANCE='NO')leaf%p_log(ii)
                  
               End Select
               
               If ( mod(pos,10) == 0) then
                  Write(un_lf, "('|')")
                  Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
                  pos = 0
               End If
               pos = pos + 1
               
            End Do
            
            Do ii = 1, 10-mod(leaf%dat_no,10_pd_ik)
               Write(un_lf, "(17(' '))", ADVANCE='NO')
            End Do
            Write(un_lf, "('|')")

         !***********************************************************
         !** dat_no in (30,1000] ************************************
         Else if ((leaf%dat_no > 30) .AND. (leaf%dat_no <= 1000)) then

            Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
            
            If ( leaf%dat_ty /= 6 ) then
               
               pos = 1
               Do ii = 1, leaf%dat_no
                  
                  Select Case (leaf%dat_ty)
                     
                  Case (1)
                     Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int1(ii)
                  Case (2)
                     Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int2(ii)
                  Case (3)
                     Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int4(ii)
                  Case (4)
                     Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int8(ii)
                  Case (5)
                     Write(un_lf, "(E16.6,' ')", ADVANCE='NO')leaf%p_real8(ii)
                  Case (7)
                     Write(un_lf, "(L16,' ')", ADVANCE='NO')leaf%p_log(ii)
                     
                  End Select
                  
                  If ( mod(pos,10) == 0) then
                     Write(un_lf, "('|')")
                     Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
                     pos = 0
                  End If
                  pos = pos + 1
                  
               End Do
               
               Do ii = 1, 10-mod(leaf%dat_no,10_pd_ik)
                  Write(un_lf, "(17(' '))", ADVANCE='NO')
               End Do
               Write(un_lf, "('|')")
         
            Else
               
               pos  = 1
               cpos = 1
               Do ii = 1, leaf%dat_no, 17
                  
                  if (Ichar(leaf%p_char(ii)) /= 10) then
                     Write(un_lf, "(*(A1))", ADVANCE='NO')leaf%p_char(ii:min(ii+16,leaf%dat_no))
                  Else
                     Write(un_lf, "(A1)", ADVANCE='NO')'_'
                  End if
                  
                  If ( mod(pos,10) == 0) then
                     Write(un_lf, "('|')")
                     Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
                     pos = 0
                  End If
                  pos = pos + 1
                  cpos = ii   
               End Do
               
               If ( mod(leaf%dat_no,17) /= 0 ) then
                  
                  Do ii = 1, 17-mod(leaf%dat_no,17)
                     Write(un_lf, "(' ')", ADVANCE='NO')
                  End Do
                  
                  Do ii = 1, 9-mod(leaf%dat_no/17,10_pd_ik)
                     Write(un_lf, "(17(' '))", ADVANCE='NO')
                  End Do
                  Write(un_lf, "('|')")
                  
               else
                  
                  Do ii = 1, 10-mod(leaf%dat_no/17,10_pd_ik)
                     Write(un_lf, "(17(' '))", ADVANCE='NO')
                  End Do
                  Write(un_lf, "('|')")
                  
               end If
               
            End If

         !***********************************************************
         !** dat_no > 1000 ******************************************
         Else

            Write(un_lf, "('"//fmt_str//"   |')",ADVANCE='NO')
            
            pos = 1
            Do ii = 1, 9
               Select Case (leaf%dat_ty)
                  
               Case (1)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int1(ii)
               Case (2)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int2(ii)
               Case (3)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int4(ii)
               Case (4)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int8(ii)
               Case (5)
                  Write(un_lf, "(E16.6,' ')", ADVANCE='NO')leaf%p_real8(ii)
               Case (6)
                  if (Ichar(leaf%p_char(ii)) /= 10) then
                     Write(un_lf, "(A1)", ADVANCE='NO')leaf%p_char(ii)
                  Else
                     Write(un_lf, "(A1)", ADVANCE='NO')'_'
                  End if
               Case (7)
                  Write(un_lf, "(L1,' ')", ADVANCE='NO')leaf%p_log(ii)
                  
               End Select
            End Do
            
            Write(un_lf, "(10(' '),'...... |')")
            Write(un_lf, "('"//fmt_str//"   |')"   ,ADVANCE='NO')
            Write(un_lf, "(10(' '),'...... ')", ADVANCE='NO')
            
            Do ii = ( leaf%dat_no-1)/2-3, ( leaf%dat_no-1)/2+4
               Select Case (leaf%dat_ty)
                  
               Case (1)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int1(ii)
               Case (2)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int2(ii)
               Case (3)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int4(ii)
               Case (4)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int8(ii)
               Case (5)
                  Write(un_lf, "(E16.6,' ')", ADVANCE='NO')leaf%p_real8(ii)
               Case (6)
                  if (Ichar(leaf%p_char(ii)) /= 10) then
                     Write(un_lf, "(A1)", ADVANCE='NO')leaf%p_char(ii)
                  Else
                     Write(un_lf, "(A1)", ADVANCE='NO')'_'
                  End if
               Case (7)
                  Write(un_lf, "(L1,' ')", ADVANCE='NO')leaf%p_log(ii)
                  
               End Select
            End Do
            
            Write(un_lf, "(10(' '),'...... |')")
            Write(un_lf, "('"//fmt_str//"   |')"  ,ADVANCE='NO')
            Write(un_lf, "(10(' '),'...... ')", ADVANCE='NO')
            
            Do ii = leaf%dat_no - 8,  leaf%dat_no
               Select Case (leaf%dat_ty)
                  
               Case (1)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int1(ii)
               Case (2)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int2(ii)
               Case (3)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int4(ii)
               Case (4)
                  Write(un_lf, "(I16,' ')", ADVANCE='NO')leaf%p_int8(ii)
               Case (5)
                  Write(un_lf, "(E16.6,' ')", ADVANCE='NO')leaf%p_real8(ii)
               Case (6)
                  if (Ichar(leaf%p_char(ii)) /= 10) then
                     Write(un_lf, "(A1)", ADVANCE='NO')leaf%p_char(ii)
                  Else
                     Write(un_lf, "(A1)", ADVANCE='NO')'_'
                  End if
               Case (7)
                  Write(un_lf, "(L1,' ')", ADVANCE='NO')leaf%p_log(ii)
                  
               End Select
            End Do
            
            Write(un_lf, "('|')")
            
         End If
         
         Write(un_lf, "('"//fmt_str//"   +',170('-'),('+'))")

      End If
      
    End Subroutine Log_Leaf_data

  End Subroutine log_leaf
  !> @} 
  !# End of memeber group "Puredat logging routines ###########################

  !############################################################################
  !############################################################################
  !> \name Puredat private auxilliary routines
  !> @{
  !> These routines exist in most codes. So they are declared private to aviod
  !> interference with other implementations
   
  !============================================================================
  !> Function which returns new free unit
  function pd_give_new_unit() result(new_unit)

    Integer :: new_unit

    Integer :: ii

    Logical :: unit_is_open

    Do ii = 3000, huge(new_unit)-1

       inquire(unit=ii, opened=unit_is_open)

       if( .not.unit_is_open ) then
          new_unit = ii
          Exit
       end if

    End Do

    if ( unit_is_open ) then

       WRITE(pd_umon,*)
       WRITE(pd_umon,*)'Something bad and unexpected happened during search ',&
            'for free unit'
       WRITE(pd_umon,*)'Could not find a new unit between 3000 and huge(Int(kind=4))'
       WRITE(pd_umon,*)' '
       WRITE(pd_umon,*)'PROGRAM STOPPED'
       STOP
    END IF

  End function pd_give_new_unit

  !============================================================================
  !> Subroutine for I/O error handling while operating on files
  SUBROUTINE file_err(in_file,io_stat, called, routine)

    INTEGER(kind=pd_ik)        , Intent(in)   :: io_stat
    CHARACTER (LEN=*)          , Intent(in)   :: in_file
    CHARACTER (LEN=*), optional, Intent(in)   :: called, routine

    IF (io_stat /= 0) Then
    
       WRITE(pd_umon,*)
       WRITE(pd_umon,"(80('='))")
       WRITE(pd_umon,"('EE ',A,T77,' EE')")   'Operation on file: '       
       WRITE(pd_umon,"('EE ',A          )")   in_file
       WRITE(pd_umon,"('EE ',A,T77,' EE')",Advance="NO") 'faild'

       If (present(called)) then
          Write(pd_umon,"('EE ',A,T77,' EE')")'during call to'
          Write(pd_umon,"('EE ',A,T77,' EE')")called
       Else
          Write(pd_umon,*)
          Write(pd_umon,"('EE ',A,T77,' EE')")'!!'
       End If
       If (present(routine)) then
          Write(pd_umon,"('EE ',A,T77,' EE')")'in '
          Write(pd_umon,"('EE ',A,T77,' EE')")routine
       End If
       
       WRITE(pd_umon,"('EE ',A,I0,T77,' EE')")'With I/O Status ',io_stat
       WRITE(pd_umon,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
       STOP
       
    End IF

  END SUBROUTINE file_err
  
  !============================================================================
  !> Subroutine for allocation error handling
  Subroutine alloc_error(alloc_stat, field, routine, dim)

   Integer(kind=pd_ik), Intent(in)           :: alloc_stat
   Character(Len=*)   , Intent(in)           :: field, routine
   Integer(kind=pd_ik), Intent(in) ,optional :: dim

   IF (alloc_stat /= 0) Then
      WRITE(mssg, '(3A,I2,3A)') "Allocation of the field/structure ", field, &
      " of dimension ", dim," in routine ", TRIM(routine), " failed."
      CALL handle_err(pd_umon, mssg, io_stat)
   End IF

  End Subroutine alloc_error

  !============================================================================
  !> Subroutine for deallocation error handling
  Subroutine dealloc_error(alloc_stat, field, routine)

    Integer(kind=pd_ik), Intent(in)           :: alloc_stat
    Character(Len=*)   , Intent(in)           :: field, routine

   IF (alloc_stat /= 0) Then
      WRITE(mssg, '(5A)') "Dellocation of the field/structure ", field, &
      " in routine ", TRIM(routine), " failed."
      CALL handle_err(pd_umon, mssg, io_stat)
   End IF

  End Subroutine dealloc_error

  !============================================================================
  !> Subroutine for consistency error handling
  Subroutine cons_error(routine, msg)

   Character(Len=*)   , Intent(in)           :: routine, msg

   WRITE(mssg, '(4A)') "A consistency error occoured in routine: ", TRIM(routine), " ", TRIM(msg)
   CALL handle_err(pd_umon, mssg, 1)
   
  End Subroutine cons_error

  !----------------------------------------------------------------------------
  Function char_to_str(char_arr) result(str)

    character, dimension(:) , intent(in) :: char_arr

    character(len=size(char_arr))        :: str

    Integer :: ii
    
    str = ''

    Do ii = 1, size(char_arr)
       if (char_arr(ii) == Char(10)) exit
       str(ii:ii) = char_arr(ii)
    End Do

  End Function char_to_str

  !----------------------------------------------------------------------------
  Function str_to_char(str) result(char_arr)

    character(len=*)  , intent(in) :: str

    character, dimension(len(str)) :: char_arr

    Integer :: ii
    
    char_arr = ''

    Do ii = 1, size(char_arr)
       char_arr(ii) = str(ii:ii)
    End Do

  End Function str_to_char
 
  !> @}
  !# End of memeber group "Puredat private auxilliary routines" ###############

End Module puredat
