!******************************************************************************
!>  Program for                                                              **
!**                                                                          **
!** ------------------------------------------------------------------------ **
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on : 12.03.2022
!** ------------------------------------------------------------------------ **
Program pd_dump_tree

  Use puredat         ! From libpuredat
  USE user_interaction

  Implicit None

  !****************************************************************************
  !** Declarations ************************************************************

  !-- Chain Variables ---------------------------------------------------------
  Real(Kind=pd_rk)            :: gstart_time, gend_time

  Integer                     :: num_args, un

  Character(len=pd_mcl)       :: outpath, outfile, arg
  Logical                     :: dump_data = .FALSE.

  Type(tBranch)               :: tree

  !== Code ====================================================================

  Call Cpu_time(gstart_time)

  !============================================================================

  num_args = command_argument_count()

  If ((num_args < 3) .or. (num_args > 4)) then
     WRITE(*, FMT_TXT_SEP)
     WRITE(*, FMT_TXT) "Usage:"
     WRITE(*, FMT_TXT) "arg 1: Output Meta path and basename of tree to be logged."
     WRITE(*, FMT_TXT) "arg 2: Project name - currently only 'results'."
     WRITE(*, FMT_TXT) "arg 3: Path of output file."
     WRITE(*, FMT_TXT) "arg 3: Name of output file."
     WRITE(*, FMT_TXT_SEP)
     STOP
  End If

  call get_command_argument(1, pro_path)
  call get_command_argument(2, pro_name)
  call get_command_argument(3, outpath)
  call get_command_argument(4, outfile)

  if (num_args > 3) then
     call get_command_argument(4, arg)
     Read(arg,*)dump_data
  end if

  write(*, FMT_TXT_SEP)
  write(*, FMT_TXT) "Pro path: ",trim(pro_path)
  write(*, FMT_TXT) "Pro name: ",trim(pro_name)
  write(*, FMT_TXT) "Out path: ",trim(outpath)
  write(*, FMT_TXT) "Out file: ",trim(outfile)
  If (num_args > 4)   write(*,*)"dump data   : ",dump_data
  write(*, FMT_TXT_SEP)

  write(*,'(A,$)')"Reading tree ... "
  tree = read_tree()
  write(*,'(A)')"done"

  if (dump_data) then
     write(*,'(A,$)')"Reading data ... "
     call read_streams(tree)
     write(*,'(A)')"done"
     write(*,'(A,$)')"Connecting pointers ... "
     call connect_pointers(tree%streams, tree)
     write(*,'(A)')"done"
  End if

  un = pd_give_new_unit()
  Open(unit=un,file=trim(outpath)//trim(outfile),action="write", status="replace")

  write(*,'(A,$)')"Logging tree ... "
  call log_tree(tree=tree, unit_lf=un, data=dump_data)
  write(*,'(A)')"done"
  !============================================================================

  Call Cpu_time(gend_time)
  Write(*,*)'=='
  Write(*,'(A,F9.3)')" Used CPU time : ",gend_time-gstart_time
 
End Program pd_dump_tree
