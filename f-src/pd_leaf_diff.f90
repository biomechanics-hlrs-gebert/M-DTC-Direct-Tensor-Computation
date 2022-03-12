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
Program pd_leaf_diff

  Use puredat         ! From libpuredat

  Implicit None

  !****************************************************************************
  !** Declarations ************************************************************

  !-- Chain Variables ---------------------------------------------------------
  Real(Kind=pd_rk)      :: gstart_time, gend_time
  Integer               :: num_args
  Character(len=pd_mcl) :: pro_path_A,pro_path_B,pro_name_A,pro_name_B
  Type(tBranch)         :: tree_A, tree_B

  !== Code ====================================================================
  Call Cpu_time(gstart_time)
  !============================================================================

  num_args = command_argument_count()

  If (num_args < 4) then
     Write(*,'(80("="))')
     Write(*,'(A)')"== Usage:"
     Write(*,'(A)')"== arg 1: Puredat project path of tree A"
     Write(*,'(A)')"== arg 2: Puredat project name of tree A"
     Write(*,'(A)')"== arg 3: Puredat project path of tree B"
     Write(*,'(A)')"== arg 4: Puredat project name of tree B"
     Write(*,'(80("="))')
     Stop
  End If

  call get_command_argument(1, pro_path_A)
  call get_command_argument(2, pro_name_A)
  call get_command_argument(3, pro_path_B)
  call get_command_argument(4, pro_name_B)

  write(*,*)"=="
  write(*,*)"Pro path A: ",trim(pro_path_A)
  write(*,*)"Pro name A: ",trim(pro_name_A)
  write(*,*)"Pro path B: ",trim(pro_path_B)
  write(*,*)"Pro name B: ",trim(pro_name_B)
  write(*,*)"=="

  pro_name = pro_name_A
  pro_path = pro_path_A

  tree_A = read_tree()

  call read_streams(tree_A)
  call connect_pointers(tree_A%streams,tree_A)

  pro_name = pro_name_B
  pro_path = pro_path_B

  tree_B = read_tree()

  call read_streams(tree_B)
  call connect_pointers(tree_B%streams,tree_B)

  write(*,*)"!!! Fixed code !!! Implement correctly !!!"
  write(*,*)"MINVAL of diff =",minval(tree_B%branches(2)%branches(1)%leaves(1)%p_real8 - &
                                      tree_A%branches(2)%branches(1)%leaves(1)%p_real8)
  write(*,*)"MAXVAL of diff =",maxval(tree_B%branches(2)%branches(1)%leaves(1)%p_real8 - &
                                      tree_A%branches(2)%branches(1)%leaves(1)%p_real8)

  !============================================================================

  Call Cpu_time(gend_time)
  Write(*,*)'=='
  Write(*,'(A,F9.3)')" Used CPU time : ",gend_time-gstart_time
 
End Program pd_leaf_diff
