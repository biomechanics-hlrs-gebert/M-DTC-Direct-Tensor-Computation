!******************************************************************************
!>  Program for                                                              **
!**                                                                          **
!** ------------------------------------------------------------------------ **
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Ralf Schneider \n
!>  on : 15.02.2012
!** ------------------------------------------------------------------------ **
Program pd_merge_branch_to_tree

  Use puredat         ! From libpuredat

  implicit none

  !****************************************************************************
  !** Declarations ************************************************************

  Real(Kind=pd_rk)            :: gstart_time, gend_time

  Character(len=pd_mcl)       :: base_path, base_name, source_path, source_name
  Character(len=pd_mcl)       :: base_desc, source_desc, arg

  Integer                     :: num_args, base_no, source_no,base_num_b,source_num_b

  Logical                     :: success

  Type(tBranch)               :: base_tree, source_tree
  Type(tBranch), Pointer      :: base_branch, source_branch

  !== Code ====================================================================

  Call Cpu_time(gstart_time)

  !============================================================================

  num_args = command_argument_count()

  If (num_args < 6) then
     Write(*,'(80("="))')
     Write(*,'(A)')"== Usage:"
     Write(*,'(A)')"== arg 1    : Puredat project path of base tree to merge to"
     Write(*,'(A)')"== arg 2    : Puredat project name of base tree to merge to"
     Write(*,'(A)')"== arg 3    : Puredat project path of source tree to merge from"
     Write(*,'(A)')"== arg 4    : Puredat project name of source tree to merge from"
     Write(*,'(A)')"== arg 5    : Description of base branch to be merged to"
     Write(*,'(A)')"== arg 6    : Description of source branch to be merged from"
     Write(*,'(A)')"==[arg 7]   : Number of base branch to extract if there are equal subbranches in the base tree"
     Write(*,'(A)')"==[arg 8]   : Number of source branch to extract if there are equal subbranches in the source tree"
     Write(*,'(80("="))')
     stop
  End If

  call get_command_argument(1, base_path)
  call get_command_argument(2, base_name)

  call get_command_argument(3, source_path)
  call get_command_argument(4, source_name)

  call get_command_argument(5, base_desc)
  call get_command_argument(6, source_desc)

  If (num_args > 6) then
     call get_command_argument(7, arg)
     read(arg,*)base_no
  End If
  If (num_args > 7) then
     call get_command_argument(8, arg)
     read(arg,*)source_no
  End If

  write(*,*)"=="
  write(*,*)"pd-Base path       : ",trim(base_path)
  write(*,*)"pd-Base name       : ",trim(base_name)
  write(*,*)"pd-Source path     : ",trim(source_path)
  write(*,*)"pd-Source name     : ",trim(source_name)
  write(*,*)"Base description   : ",trim(base_desc)
  write(*,*)"Source description : ",trim(source_desc)
  If  (num_args > 6) write(*,*)"Base   branch no   : ",base_no
  If  (num_args > 7) write(*,*)"Source branch no   : ",source_no

  !-- Read trees --------------------------------------------------------------
  pro_path  = base_path
  pro_name  = base_name
  base_tree = read_tree()

  pro_path  = source_path
  pro_name  = source_name
  source_tree = read_tree()

  !----------------------------------------------------------------------------
  !-- Search for base branch to include source branch into --------------------

  !-- Check whether there are multiple branches with the same description in --
  !-- the base tree                                                          --
  call get_branch_num(trim(base_desc),base_tree,base_num_b)

  !-- Case 1: Only one branch is found ----------------------------------------
  if (base_num_b == 1) then
     call search_branch(trim(base_desc),base_tree,base_branch,success)
     If  (num_args > 6) then
        Write(*,*)"The given global number for the base branch was overwritten by"
        Write(*,*)"the search for the given description of the base branch since"
        Write(*,*)"it was only found one branch with the given description."
     End If

  !-- Case 2: Multiple branches with the same description are found -----------
  else
     If  (num_args == 6) then
        Write(*,*)"Multiple branches with the given description were found."
        Write(*,*)"Since no global number for the branch extraction was given"
        Write(*,*)"the first branch found will be returned."
        call search_branch(trim(base_desc),base_tree,base_branch,success)

     else
        write(*,*)"Getting branch with no ",base_no
        call get_branch_with_num(trim(base_desc), base_tree, base_no, base_branch, success)

     End If

  End if
  If (.NOT. success) then
     write(*,*)"Couldn't find the branch ",trim(base_desc)
     Write(*,*)"given as argument 5      "
     write(*,*)"in base tree             ",trim(base_name)
     write(*,*)"with pro_path            ",trim(base_path)
     write(*,*)"Program halted !!!!!!!!!!!"
     STOP
  End If

  !----------------------------------------------------------------------------
  !-- Search for source branch to include into base branch --------------------

  !-- Check whether there are multiple branches with the same description in --
  !-- the source tree  
  call get_branch_num(trim(source_desc),source_tree,source_num_b)

  !-- Case 1: Only one branch is found ----------------------------------------
  if (source_num_b == 1) then
     call search_branch(trim(source_desc),source_tree,source_branch,success)
     If  (num_args > 7) then
        Write(*,*)"The given global number for the source branch was overwritten by"
        Write(*,*)"the search for the given description of the source branch since"
        Write(*,*)"it was only found one branch with the given description."
     End If

  !-- Case 2: Multiple branches with the same description are found -----------
  else
     If  (num_args == 7) then
        Write(*,*)"Multiple branches with the given description were found."
        Write(*,*)"Since no global number for the branch extraction was given"
        Write(*,*)"the first branch found will be returned."
        call search_branch(trim(source_desc),source_tree,source_branch,success)

     else
        write(*,*)"Getting branch with no ",source_no
        call get_branch_with_num(trim(source_desc), source_tree, source_no, source_branch, success)

     End If

  End if
  If (.NOT. success) then
     write(*,*)"Couldn't find the branch ",trim(source_desc)
     Write(*,*)"given as argument 6      "
     write(*,*)"in base tree             ",trim(source_name)
     write(*,*)"with pro_path            ",trim(source_path)
     write(*,*)"Program halted !!!!!!!!!!!"
     STOP
  End If

  call include_branch_into_branch(s_b=source_branch, t_b=base_branch, &
                                  s_streams = source_tree%streams   , &
                                  t_streams = base_tree%streams     , &
                                  clean_target_files = .TRUE.          )
  pro_path  = base_path
  pro_name  = base_name
  call write_tree(base_tree)

  !============================================================================

  Call Cpu_time(gend_time)
  Write(*,*)'=='
  Write(*,'(A,F9.3)')" Used CPU time : ",gend_time-gstart_time
  
End Program pd_merge_branch_to_tree
