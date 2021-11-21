!------------------------------------------------------------------------------
! MODULE: meta
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! @Description:
!> Module containing all meta file read/write routines.
!
! REVISION HISTORY:
! 21 10 2021 - Initial refactored version
!------------------------------------------------------------------------------
MODULE error_handling

   USE global_std
   USE puredat_globals
   USE strings
   USE MPI

IMPLICIT NONE

CONTAINS

   !------------------------------------------------------------------------------
! SUBROUTINE: handle_err
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print and handle error messages. Mpi handling written by 
!> Ralf Schneider.
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
SUBROUTINE handle_err(fh, txt, err) ! , pro_path, pro_name

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

! Best while containing scl characters! scl = 64 + "EE " + 2
CHARACTER(LEN=*), PARAMETER :: FMT_WRN = "('WW ',A,T68,' WW')"  
CHARACTER(LEN=*), PARAMETER :: FMT_ERR = "('EE ',A,T68,' EE')"
CHARACTER(LEN=*), PARAMETER :: FMT_WRN_SO = "('\x1B[33m','WW ','\x1B[0m',A, T76,'\x1B[33m',' WW','\x1B[0m')" ! std_out
CHARACTER(LEN=*), PARAMETER :: FMT_ERR_SO = "('\x1B[31m','EE ','\x1B[0m',A, T76,'\x1B[31m',' EE','\x1B[0m')" ! std_out

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
   IF ((err == -1) .AND. (fh /= std_out)) fmt = FMT_WRN
   IF ((err == -1) .AND. (fh == std_out)) fmt = FMT_WRN_SO

   IF (txt  /= '') THEN
      ! Parse error message
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

   END IF ! (txt  /= '') THEN

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

END SUBROUTINE handle_err



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

END MODULE error_handling