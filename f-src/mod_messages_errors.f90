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


!------------------------------------------------------------------------------
! Still not obsolete formatters
!------------------------------------------------------------------------------
! Error formats
Character(Len=*), Parameter :: FMT_ERR_AI0 = "('EE ',*(A,I0))"  
CHARACTER(Len=*), PARAMETER :: FMT_ERR_A   = "('EE ',A)"

! Message formats
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0  = "('MM ',*(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_2AI0 = "('MM ',2(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3I0 = "('MM ',A,3(',',I0))"

CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0  = "('MM ',A,F0.6)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3F0 = "('MM ',A,3(',',F0.6))"


!------------------------------------------------------------------------------
! Provide colors on std_out (!) 
! Needs to compile with -fbackslash 
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

CALL print_sep (std_out)
CALL print_std (std_out, 'High Performance Computing Center | Stuttgart (HLRS)')
CALL print_sep (std_out)
CALL print_std (std_out, 'Directly Discretizing Tensor Computation '//TRIM(ADJUSTL(revision)))
CALL print_std (std_out, '')
CALL print_std (std_out, 'Author: Dr.-Ing. Ralf Schneider (HLRS, NUM)')
CALL print_std (std_out, 'Author: Johannes Gebert, M.Sc.  (HLRS, NUM)')
CALL print_sep (std_out)
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

CALL print_sep (std_out)
CALL print_std (std_out, 'Directly Discretizing Tensor Computation | Usage:')
CALL print_sep (std_out)
CALL print_std (std_out, './ddtc_vx.y.z_x86_64 »flags« »basename.meta«')
CALL print_std (std_out, '')
CALL print_std (std_out, '-h/ --help      This message.')
CALL print_std (std_out, '-v/ --version   Version of the program')
CALL print_std (std_out, '--restart       Overwrite restart keyword')
CALL print_std (std_out, '--no-restart    Overwrite restart keyword')
CALL print_std (std_out, '')
CALL print_std (std_out, 'The meta-file must be the last command argument.')
CALL print_sep (std_out)

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


IF (fh == 6) THEN
   leading_fmt = "('"//TRIM(tag_color)//TRIM(tag_label)//"  "//TRIM(FMT_nocolor)//"', A, T"
   trailing_fmt = ", '"//TRIM(tag_color)//"  "//TRIM(tag_label)//TRIM(FMT_nocolor)//"')"  
   WRITE(imw, '(I0)') mw + LEN_TRIM(TRIM(tag_color)) + LEN_TRIM(tag_label) + 1
ELSE
   leading_fmt = "('"//TRIM(tag_label)//"  "//"', A, T"
   trailing_fmt = ", '"//"  "//TRIM(tag_label)//"')"  
   WRITE(imw, '(I0)') mw + LEN_TRIM(tag_label) + 2 - LEN_TRIM(tag_color)
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
SUBROUTINE print_sep(fh, chara)

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: chara

CHARACTER(LEN=1) :: separator

separator='-'
IF (PRESENT(chara)) separator = TRIM(chara)

WRITE (fh, '(A)') REPEAT(separator, mw + 2) ! +2 for spaces in formatting
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
! SUBROUTINE: print_debug
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
SUBROUTINE print_debug(fh, txt)

INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt

CHARACTER(LEN=scl) :: DBG_COLOR
CHARACTER(LEN=scl) :: DBG_LBL

DBG_LBL   = "DBG"
DBG_COLOR = FMT_nocolor

CALL print_tag(fh, txt, DBG_LBL, DBG_COLOR)

END SUBROUTINE print_debug


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