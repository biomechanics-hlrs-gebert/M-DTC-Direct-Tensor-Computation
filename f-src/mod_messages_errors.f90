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
   USE puredat_globals
   USE strings
   USE MPI

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: FMT_WRN_SO = "('\x1B[33m','WW ','\x1B[0m',A, T76,'\x1B[33m',' WW','\x1B[0m')" ! std_out
CHARACTER(LEN=*), PARAMETER :: FMT_ERR_SO = "('\x1B[31m','EE ','\x1B[0m',A, T76,'\x1B[31m',' EE','\x1B[0m')" ! std_out

!------------------------------------------------------------------------------
! Standard formats
!------------------------------------------------------------------------------
CHARACTER(len=5)   , PARAMETER :: creturn = achar(13)

! Seperators
CHARACTER(Len=*), PARAMETER :: SEP_STD= "(68('-'))" ! FMT_HY_SEP
CHARACTER(Len=*), PARAMETER :: FMT_EQ_SEP  = "(68('='))"
CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP = "('#DBG#',95('='))"

! Character constants for nice output ---------------------------------------
Character(LEN=*), Parameter :: fmt_inpsep = "('+',99('-'))"

Character(Len=*), Parameter :: FMT_MSG     = "('MM ',A,T97,' MM')"
Character(Len=*), Parameter :: FMT_MSG_BS  = "('MM ',A,T88,' ... ',$)"
Character(Len=*), Parameter :: FMT_MSG_BE  = "('done MM')"

CHARACTER(LEN=*), PARAMETER :: FMT_WRN     = "('WW ',A,T68,' WW')"  
CHARACTER(LEN=*), PARAMETER :: FMT_ERR     = "('EE ',A,T68,' EE')"


Character(Len=*), Parameter :: FMT_ERR_AI0 = "('EE ',*(A,I0))"  
CHARACTER(Len=*), PARAMETER :: FMT_ERR_A   = "('EE ',A)"

Character(Len=*), Parameter :: FMT_STOP    = "('EE PROGRAM STOPPED ..... ',&
                                             &T97,' EE',/,'<',97('='),'>')"

Character(Len=*), Parameter :: FMT_TIME = "('MM ',A,1X,F0.6,' sec')"


! Warning formats
CHARACTER(Len=*), PARAMETER :: FMT_WRN_A    = "('WW ',A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0  = "('WW ',A,1X,I0)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0A = "('WW ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AF0  = "('WW ',A,1X,F0.6)"

! Message formats
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0  = "('MM ',*(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0A = "('MM ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A8I5 = "('MM ',A,1X,8(',',I5))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_2AI0 = "('MM ',2(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3I0 = "('MM ',A,3(',',I0))"

CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0  = "('MM ',A,F0.6)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0A = "('MM ',A,F0.6,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A2F0 = "('MM ',A,2(',',F0.6))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3F0 = "('MM ',A,3(',',F0.6))"

CHARACTER(Len=*), PARAMETER :: FMT_MSG_AL  = "('MM ',A,L1)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A   = "('MM ',A)"



! PureDat Formatters
Character(Len=*), Parameter :: PDF_E_A    = "('EE ',A)"
Character(Len=*), Parameter :: PDF_E_AI0  = "('EE ',*(A,1X,I0))"
Character(Len=*), Parameter :: PDF_E_STOP = "('EE PROGRAM STOPPED ..... ',/,'<',98('='),'>')"

Character(Len=*), Parameter :: PDF_W_A    = "('WW ',A)"
Character(Len=*), Parameter :: PDF_W_AI0  = "('WW ',*(A,1X,I0))"

Character(Len=*), Parameter :: PDF_M_A    = "('MM ',A)"
Character(Len=*), Parameter :: PDF_M_AI0  = "('MM ',A,1X,I0)"

Character(Len=*), Parameter :: PDF_TIME   = "('MM ',A,1X,F0.6,' sec')"

Character(Len=*), Parameter :: PDF_SEP    = "('<',98('='),'>')"

! Provide colors on std_out (!)
CHARACTER(LEN=*), PARAMETER ::  FMT_Blck  = "\x1B[30m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Red   = "\x1B[31m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Grn   = "\x1B[32m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Orng  = "\x1B[33m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Blue  = "\x1B[34m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Prpl  = "\x1B[35m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Cyan  = "\x1B[36m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Gray  = "\x1B[37m"
CHARACTER(LEN=*), PARAMETER ::  FMT_noc   = "\x1B[0m"

CONTAINS


!------------------------------------------------------------------------------
! SUBROUTINE: print_message
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Accepts a text and a formatter. Prints the message with a proper formatting 
!> and width. Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] fmt Message formatting
!> @param[in] txt Text to print
!------------------------------------------------------------------------------  
SUBROUTINE print_message(fh, fmt, txt) ! , pro_path, pro_name

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: fmt
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=scl) :: sub_mssg

CHARACTER(LEN=mcl)   :: delim, tokens(100), path_tokens(50)
CHARACTER(LEN=mcl+1) :: next_token
INTEGER  (KIND=ik)   :: ntokens, path_ntokens
INTEGER  (KIND=ik)   :: ii, jj, sw, mode

text = TRIM(ADJUSTL(txt))
   
mode = 0                ! Absolute or relative path
sw = 2                  ! Whether it's the beginning or within a path
ntokens = 0             ! Amount of words in message
path_ntokens = 0        ! Amount of words in a path
delim = '/'
ii = 1
jj = 1

IF (txt  /= '') THEN
   ! Parse message
   CALL parse(str=text, delims=' ', args = tokens, nargs=ntokens)

   ! next_token  = tokens(1) 
   next_token = ''
   WRITE(fh, fmt) REPEAT('-',scl)

   DO WHILE (ii .LT. ntokens) 
   
         sub_mssg = REPEAT(' ', scl)
         sub_mssg = TRIM(next_token)

         DO           
            ! path_ntokens = 1
            IF (sw==2) CALL parse(str = tokens(ii), delims='/', args = path_tokens, nargs = path_ntokens)

            IF (path_ntokens .GT. 1) sw=1
            
            IF (sw == 1) THEN
               IF (TRIM(ADJUSTL(path_tokens(1))) =='') mode = 2
               
               IF ((mode == 2) .AND. (jj == 1)) jj = jj + 1
               IF ((mode == 2) .AND. (jj == 2)) THEN
                     delim = ' /'
               ELSE
                     delim = '/'
               END IF

               next_token = TRIM(delim)//ADJUSTL(path_tokens(jj))

               jj = jj + 1                         
               IF (jj .GT. path_ntokens) THEN
                     sw   = 2
                     jj   = 1
                     mode = 1
                     ii   = ii + 1
               END IF
            ELSE
               next_token = ' '//tokens(ii)
               ii = ii + 1
               IF (ii .GT. ntokens+1) EXIT
            END IF
                        
            IF ((LEN_TRIM(ADJUSTL(sub_mssg)) + LEN_TRIM(next_token)) .LT. scl) THEN
               sub_mssg = TRIM(ADJUSTL(sub_mssg))//TRIM(next_token)
            ELSE
               EXIT ! Finishes the current line with scl characters
            END IF
         END DO

         WRITE(fh, fmt) sub_mssg

   END DO
   FLUSH(fh)
END IF
END SUBROUTINE print_message

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

INTEGER(KIND=ik) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt
INTEGER(KIND=ik) :: err

INTEGER(KIND=mpi_ik) :: ierr = 0
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt
CHARACTER(LEN=scl) :: sub_mssg, errnmbr

CHARACTER(LEN=mcl)   :: delim, tokens(100), path_tokens(50)
CHARACTER(LEN=mcl+1) :: next_token
INTEGER  (KIND=ik) :: ntokens, path_ntokens
INTEGER  (KIND=ik) :: ii, jj, sw, mode

text = TRIM(ADJUSTL(txt))
fmt  = FMT_ERR
   
mode = 0                ! Absolute or relative path
sw = 2                  ! Whether it's the beginning or within a path
ntokens = 0             ! Amount of words in message
path_ntokens = 0        ! Amount of words in a path
delim = '/'
ii = 1
jj = 1

! if err == 0 is a warning --> all "mpi-errors" are warnings :-)
IF (err == 0) THEN
   CONTINUE
ELSE
   ! Colorized text in case of printing to std out which might be a terminal.
   ! Someone may write it in a nice way :-)
   IF ((err >=  0) .AND. (fh == std_out)) fmt = FMT_ERR_SO

   CALL print_message(fh, fmt, txt)
  
   IF (err .GT. 0) THEN ! An error occured   
      WRITE(fh, fmt) REPEAT('-',scl)

      WRITE(errnmbr, '(I15)') INT(err, KIND=ik)

      mssg = "with error code "//TRIM(ADJUSTL(errnmbr))//"."

      IF (fh == std_out) THEN
            WRITE(fmt,'(A,I5,A)') "('\x1B[31m','EE Program halted ','\x1B[0m',A, T",&
               scl+LEN('\x1B[31m'//'EE '//'\x1B[0m'),",'\x1B[31m',' EE','\x1B[0m')"
            WRITE(fh, fmt) TRIM(mssg)
      ELSE
            WRITE(fmt,'(A,I5,A)') "('EE Program halted ',A, T",scl + LEN('EE ') + 1,",' EE')"
            WRITE(fh, fmt) TRIM(mssg)
      END IF

      CALL stop_slaves()

      CALL MPI_FINALIZE(ierr)

      IF ( ierr /= 0 ) WRITE(std_out,'(A)') "MPI_FINALIZE did not succeed"
      STOP

   END IF

END IF

END SUBROUTINE print_err_stop


!------------------------------------------------------------------------------
! SUBROUTINE: print_warning
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print warning.
!
!> @Description
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!------------------------------------------------------------------------------  
SUBROUTINE print_warning(fh, txt)

INTEGER(KIND=ik) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=mcl) :: fmt

fmt = FMT_WRN

IF (fh == std_out) fmt = FMT_WRN_SO

CALL print_message(fh, fmt, txt)
END SUBROUTINE print_warning

!------------------------------------------------------------------------------
! SUBROUTINE: stop_slaves
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Stop slaves properly
!------------------------------------------------------------------------------  
SUBROUTINE stop_slaves()

Integer(mpi_ik) :: ierr

CALL MPI_BCAST(pro_path, INT(mcl, mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(pro_name, INT(mcl, mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

!** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
CALL MPI_BCAST(-1_ik, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)   

CALL MPI_FINALIZE(ierr)

IF ( ierr /= 0 ) WRITE(*,'(A)') "MPI_FINALIZE did not succeed"
STOP 
END SUBROUTINE stop_slaves

END MODULE messages_errors