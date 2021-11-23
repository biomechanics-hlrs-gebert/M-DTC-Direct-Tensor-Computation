!------------------------------------------------------------------------------
! MODULE: messages_errors
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! @Description:
!> Module containing the formatting of all (error) messages
!------------------------------------------------------------------------------
MODULE messages_errors

   USE global_std
   USE strings
   USE puredat_globals
   USE MPI

IMPLICIT NONE

!------------------------------------------------------------------------------
! Parameters of (error) messages
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: mw = 90  ! Message width, including leading and trailing descriptors

CHARACTER(len=5)   , PARAMETER :: creturn = achar(13)

CHARACTER(Len=*), PARAMETER :: FMT_ERR      = "('EE ', A)"
Character(Len=*), Parameter :: FMT_ERR_STOP = "('EE PROGRAM STOPPED.')"
CHARACTER(Len=*), PARAMETER :: FMT_ERR_SEP  = "('EE ', 76('='))"

Character(Len=*), Parameter :: FMT_ERR_AI0  = "('EE ', *(A,I0))"  

!------------------------------------------------------------------------------
! Text formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_TXT      = "('-- ',A)"
CHARACTER(Len=*), PARAMETER :: FMT_TXT_SEP  = "(80('-'))"


!------------------------------------------------------------------------------
! Message formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_MSG      = "('MM ',A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_SEP  = "(80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AxI0  = "('MM ',A,*(1x, I0))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0  = "('MM ',*(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0A = "('MM ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_2AI0 = "('MM ',2(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3I0 = "('MM ',A,3(1x,I0))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0  = "('MM ',A,1X,F0.6)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0A = "('MM ',A,1X,F0.6,1x,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A2F0 = "('MM ',A,2(1X,F0.6))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3F0 = "('MM ',A,3(1X,F0.6))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AL   = "('MM ',A,1x,L1)"

!------------------------------------------------------------------------------
! Warning formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_WRN      = "('WW ',A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_SEP  = "(80('-'))"

CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0  = "('WW ',A,1X,I0)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0A = "('WW ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AF0  = "('WW ',A,1X,F0.6)"

!------------------------------------------------------------------------------
! Debug formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP = "('#DBG#',75('='))"

!------------------------------------------------------------------------------
! Provide colors on std_out (!) 
! Needs to compile with -fbackslash 
! Use of a requires resetting it.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::  FMT_Blck    = "\x1B[30m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Red     = "\x1B[31m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Green   = "\x1B[32m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Orange  = "\x1B[33m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Blue    = "\x1B[34m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Purple  = "\x1B[35m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Cyan    = "\x1B[36m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Gray    = "\x1B[37m"
CHARACTER(LEN=*), PARAMETER ::  FMT_nocolor = "\x1B[0m"

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: show_title
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Show brief information about the program
!------------------------------------------------------------------------------
SUBROUTINE show_title(revision)

CHARACTER(LEN=*), INTENT(IN) :: revision

WRITE(std_out, FMT_TXT_SEP)
WRITE(std_out, FMT_TXT)     'High Performance Computing Center | Stuttgart (HLRS)'
WRITE(std_out, FMT_TXT)     ''
WRITE(std_out, FMT_TXT)     'Directly Discretizing Tensor Computation '//TRIM(ADJUSTL(revision))
WRITE(std_out, FMT_TXT)     ''     
WRITE(std_out, FMT_TXT)     'Author: Dr.-Ing. Ralf Schneider (HLRS, NUM)'
WRITE(std_out, FMT_TXT)     'Author: Johannes Gebert, M.Sc.  (HLRS, NUM)'
WRITE(std_out, FMT_TXT_SEP)
END SUBROUTINE show_title

!------------------------------------------------------------------------------
! FUNCITON: usage
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Print program usage. 
!
!> @param[in] err Optional suppression of program abortion
!------------------------------------------------------------------------------  
SUBROUTINE usage()

WRITE(std_out, FMT_TXT_SEP)
WRITE(std_out, FMT_TXT) 'Directly Discretizing Tensor Computation | Usage:'
WRITE(std_out, FMT_TXT_SEP)
WRITE(std_out, FMT_TXT) './ddtc_vx.y.z_x86_64 »flags« »basename.meta«'
WRITE(std_out, FMT_TXT) ''
WRITE(std_out, FMT_TXT) '-h/ --help      This message.'
WRITE(std_out, FMT_TXT) '-v/ --version   Version of the program'
WRITE(std_out, FMT_TXT) '--restart       Overwrite restart keyword'
WRITE(std_out, FMT_TXT) '--no-restart    Overwrite restart keyword'
WRITE(std_out, FMT_TXT) ''
WRITE(std_out, FMT_TXT) 'The meta-file must be the last command argument.'
WRITE(std_out, FMT_TXT_SEP)

END SUBROUTINE usage


!------------------------------------------------------------------------------
! FUNCITON: determine_stout
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Check whether the program can access a std_out. Usually given by the
!> environment settings. Fall back always results in writing to a file.
!
!> @return fh_std_out File handle of the »real« std_out
!------------------------------------------------------------------------------  
FUNCTION determine_stout() RESULT(fh_std_out)

INTEGER(KIND=ik) :: fh_std_out, stat
CHARACTER(LEN=scl) :: use_std_out


CALL GET_ENVIRONMENT_VARIABLE(NAME='USE_STD_OUT', VALUE=use_std_out, STATUS=stat)

IF ((stat == 0) .AND. (use_std_out == 'YES')) THEN
    fh_std_out = 6 ! Standard std_out
ELSE
    fh_std_out = give_new_unit()
END IF
 
END FUNCTION determine_stout



!------------------------------------------------------------------------------
! SUBROUTINE: give_new_unit
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Function which returns a new free unit
!------------------------------------------------------------------------------
function give_new_unit() result(new_unit)

Integer :: new_unit
Integer :: ii
Logical :: unit_is_open

Do ii = 300, huge(new_unit)-1

    inquire(unit=ii, opened=unit_is_open)

    if( .not.unit_is_open ) then
        new_unit = ii
        Exit
    end if

End Do

if ( unit_is_open ) then
    mssg = 'Something bad and unexpected happened during search for free unit: &
    &Could not find a new unit between 100 and huge(Int(kind=4))'
    CALL print_err_stop(std_out, mssg, 1)
END IF

End function give_new_unit


!------------------------------------------------------------------------------
! SUBROUTINE: print_err_stop
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print and handle error messages. Mpi handling written by 
!> Ralf Schneider. Accepts readily written strings. *No* numbers!
!
!> @Description
!> Aborts with errors.
!> Errors:  err > 0
!> Warings: err < 0
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!> @param[in] err Errorcode / status of the message
!------------------------------------------------------------------------------  
SUBROUTINE print_err_stop(fh, txt, err) ! , pro_path, pro_name

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt
INTEGER(KIND=ik) :: err

INTEGER(KIND=mpi_ik) :: ierr

! if err == 0 is a warning --> all "mpi-errors" are warnings :-)
IF (err == 0) THEN
   CONTINUE
ELSE

   CALL stop_all(ierr)

   WRITE(fh, FMT_ERR) TRIM(txt)
   WRITE(fh, FMT_ERR_STOP)

   IF ( ierr /= 0 ) WRITE(*,'(A)') "MPI_FINALIZE did not succeed"
   STOP 

END IF

END SUBROUTINE print_err_stop




!------------------------------------------------------------------------------
! SUBROUTINE: stop_all
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Stop slaves properly
!
!> @param(out) ierr MPI Error code
!------------------------------------------------------------------------------  
SUBROUTINE stop_all(ierr)

Integer(mpi_ik) :: ierr

CALL MPI_BCAST(pro_path, INT(mcl, mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(pro_name, INT(mcl, mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

!** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
CALL MPI_BCAST(-1_ik, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)   

CALL MPI_FINALIZE(ierr)

END SUBROUTINE stop_all

END MODULE messages_errors