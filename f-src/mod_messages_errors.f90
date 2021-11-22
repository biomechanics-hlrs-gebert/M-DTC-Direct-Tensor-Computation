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

CHARACTER(LEN=*), PARAMETER :: FMT_WRN_SO = "('\x1B[33m','WW ','\x1B[0m',A, T76,'\x1B[33m',' WW','\x1B[0m')" ! std_out
CHARACTER(LEN=*), PARAMETER :: FMT_ERR_SO = "('\x1B[31m','EE ','\x1B[0m',A, T76,'\x1B[31m',' EE','\x1B[0m')" ! std_out

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
! SUBROUTINE: print_text
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
SUBROUTINE print_text(fh, fmt, txt) ! , pro_path, pro_name

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: fmt
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mw)  :: sub_mssg

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

CALL parse(str=text, delims=' ', args = tokens, nargs=ntokens)

IF(ntokens<=1) THEN
   WRITE(fh, fmt) TRIM(txt)
ELSE
   next_token = ''

   DO WHILE (ii .LT. ntokens) 

         sub_mssg = ''
         sub_mssg = TRIM(next_token)

         DO           

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
                        
            !------------------------------------------------------------------------------
            ! Concatenate strings until the next word does not fit into the message width 
            !------------------------------------------------------------------------------  
            IF (((LEN_TRIM(sub_mssg) + LEN_TRIM(next_token)) ).LT. (mw)) THEN
               sub_mssg = TRIM(ADJUSTL(sub_mssg))//TRIM(next_token)
            ELSE
               EXIT ! Finishes the current line with scl characters
            END IF
         END DO

         WRITE(fh, fmt) sub_mssg

   END DO
END IF
FLUSH(fh)

END SUBROUTINE print_text

!------------------------------------------------------------------------------
! SUBROUTINE: print_tag
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print stuff.
!
!> @Description
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!------------------------------------------------------------------------------  
SUBROUTINE print_tag(fh, txt, tag_label, tag_color)

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt
CHARACTER(LEN=*), INTENT(IN) :: tag_label
CHARACTER(LEN=*), INTENT(IN) :: tag_color

CHARACTER(LEN=mcl) :: fmt, leading_fmt, trailing_fmt
CHARACTER(LEN=scl) :: imw ! Internal message width


IF (fh == std_out) THEN
   leading_fmt = "('"//TRIM(tag_color)//TRIM(tag_label)//"  "//TRIM(FMT_nocolor)//"', A, T"
   trailing_fmt = ", '"//TRIM(tag_color)//"  "//TRIM(tag_label)//TRIM(FMT_nocolor)//"')"  
   WRITE(imw, '(I0)') mw + LEN_TRIM(TRIM(tag_color)) + LEN_TRIM(tag_label) + 1
ELSE
   leading_fmt = "('"//TRIM(tag_label)//"  "//"', A, T"
   trailing_fmt = ", '"//"  "//TRIM(tag_label)//"')"  
   WRITE(imw, '(I0)') mw + LEN_TRIM(tag_label) + 1
END IF

fmt = TRIM(leading_fmt)//TRIM(imw)//TRIM(trailing_fmt)

CALL print_text(fh, fmt, txt)

END SUBROUTINE print_tag

!------------------------------------------------------------------------------
! SUBROUTINE: print_std
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print a spearator.
!
!> @Description
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!------------------------------------------------------------------------------  
SUBROUTINE print_std(fh, txt)

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=scl) :: TXT_LBL
CHARACTER(LEN=mcl) :: fmt
CHARACTER(LEN=scl) :: imw ! Internal message width

TXT_LBL   = "--"

WRITE(imw, '(I0)') mw + LEN_TRIM(TXT_LBL) + 1

fmt = "('"//TRIM(TXT_LBL)//"  "//"', A, T"//TRIM(imw)//")"

CALL print_text(fh, fmt, txt)
END SUBROUTINE print_std

!------------------------------------------------------------------------------
! SUBROUTINE: print_sep
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print a spearator.
!
!> @Description
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!------------------------------------------------------------------------------  
SUBROUTINE print_sep(fh)
INTEGER(KIND=ik), INTENT(IN) :: fh 
WRITE (fh, '(A)') REPEAT('-', mw + 2) ! +2 for spaces in formatting
END SUBROUTINE print_sep

!------------------------------------------------------------------------------
! SUBROUTINE: print_message
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print a message.
!
!> @Description
!> Colorized output in case it's printing to std_out.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!------------------------------------------------------------------------------  
SUBROUTINE print_message(fh, txt)

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=scl) :: WRN_COLOR
CHARACTER(LEN=scl) :: WRN_LBL

WRN_LBL   = "MM"
WRN_COLOR = FMT_nocolor

CALL print_tag(fh, txt, WRN_LBL, WRN_COLOR)

END SUBROUTINE print_message


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

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=scl) :: WRN_COLOR
CHARACTER(LEN=scl) :: WRN_LBL

WRN_LBL   = "WW"
WRN_COLOR = FMT_Orange

CALL print_tag(fh, txt, WRN_LBL, WRN_COLOR)

END SUBROUTINE print_warning

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

INTEGER(KIND=mpi_ik) :: ierr = 0
CHARACTER(LEN=scl) :: errnmbr

CHARACTER(LEN=mcl) :: delim
INTEGER  (KIND=ik) :: ntokens, path_ntokens
INTEGER  (KIND=ik) :: ii, jj, sw, mode

CHARACTER(LEN=mcl) :: leading_fmt, trailing_fmt
CHARACTER(LEN=scl) :: imw ! Internal message width
CHARACTER(LEN=scl) :: ERR_COLOR, ERR_COLOR_LONG
CHARACTER(LEN=scl) :: ERR_LBL, ERR_LBL_LONG

ERR_LBL   = "EE"
ERR_COLOR = FMT_Red

ERR_LBL_LONG   = "EE Program halted"
ERR_COLOR_LONG = FMT_Red


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

   CALL print_tag(fh, txt, ERR_LBL, ERR_COLOR)

   WRITE(errnmbr, '(I0)') INT(err, KIND=ik)

   mssg = "with error code "//TRIM(ADJUSTL(errnmbr))//"."

   IF (fh == std_out) THEN
      leading_fmt = "('"//TRIM(ERR_COLOR)//TRIM(ERR_LBL_LONG)//"  "//TRIM(FMT_nocolor)//"', A, T"
      trailing_fmt = ", '"//TRIM(ERR_COLOR)//"  "//TRIM(ERR_LBL)//TRIM(FMT_nocolor)//"')"  
      WRITE(imw, '(I0)') mw + LEN_TRIM(TRIM(ERR_COLOR)) + LEN_TRIM(ERR_LBL) + 1
   ELSE
      leading_fmt = "('"//TRIM(ERR_LBL_LONG)//"  "//"', A, T"
      trailing_fmt = ", '"//"  "//TRIM(ERR_LBL)//"')"  
      WRITE(imw, '(I0)') mw + LEN_TRIM(ERR_LBL) + 1
   END IF

   CALL stop_slaves()

   CALL MPI_FINALIZE(ierr)

   IF ( ierr /= 0 ) WRITE(std_out,'(A)') "MPI_FINALIZE did not succeed"
   STOP

END IF

END SUBROUTINE print_err_stop




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