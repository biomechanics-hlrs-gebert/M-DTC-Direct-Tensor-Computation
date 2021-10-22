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
Program pd_dump_tree

  Use puredat         ! From libpuredat

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

  If (num_args < 4) then
     Write(*,'(80("="))')
     Write(*,'(A)')"== Usage:"
     Write(*,'(A)')"== arg 1             : Puredat project path of tree to be logged"
     Write(*,'(A)')"== arg 2             : Puredat project name of tree to be logged"
     Write(*,'(A)')"== arg 3             : Path to output file"
     Write(*,'(A)')"== arg 4             : Name of output file"
     Write(*,'(A)')"== arg 5 [optional]  : Dump data default=.FALSE."
     Write(*,'(80("="))')
     Stop
  End If

  call get_command_argument(1, pro_path)
  call get_command_argument(2, pro_name)
  call get_command_argument(3, outpath)
  call get_command_argument(4, outfile)

  if (num_args > 4) then
     call get_command_argument(5, arg)
     Read(arg,*)dump_data
  end if

  write(*,*)"=="
  write(*,*)"Pro path    : ",trim(pro_path)
  write(*,*)"Pro name    : ",trim(pro_name)
  write(*,*)"Out path    : ",trim(outpath)
  write(*,*)"Out file    : ",trim(outfile)
  If (num_args > 4)   write(*,*)"dump data   : ",dump_data
  write(*,*)"=="

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
  Open(unit=un,file=trim(outpath)//trim(outfile),action="write",&
       status="replace")

  write(*,'(A,$)')"Logging tree ... "
  call log_tree(tree=tree, unit_lf=un, data=dump_data)
  write(*,'(A)')"done"
  !============================================================================

  Call Cpu_time(gend_time)
  Write(*,*)'=='
  Write(*,'(A,F9.3)')" Used CPU time : ",gend_time-gstart_time
 
End Program pd_dump_tree
