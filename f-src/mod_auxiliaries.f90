!------------------------------------------------------------------------------
! HLRS - Numerical Methods and Libraries, Nobelstraße 19, 70569 Stumgart
!------------------------------------------------------------------------------
!
! MODULE: auxiliaries
!
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! DESCRIPTION: 
!> Module containing useful routines.
!
! REVISION HISTORY:
! 27 09 2021 - Initial Version
!------------------------------------------------------------------------------

MODULE auxiliaries

USE global_std
USE strings
USE MPI

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: show_title
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Show brief information about the program
!------------------------------------------------------------------------------
SUBROUTINE show_title(revision, hash)

CHARACTER(LEN=*), INTENT(IN) :: revision
CHARACTER(LEN=*), INTENT(IN) :: hash

WRITE(std_out, SEP_STD)
WRITE(std_out,'(A)') '-- High Performance Computing Center | Stuttgart (HLRS)'
WRITE(std_out, SEP_STD)
WRITE(std_out,'( A)') '-- Directly Discretizing Tensor Computation'
WRITE(std_out,'( A)') '--'
WRITE(std_out,'( A)') '-- Author: Dr.-Ing. Ralf Schneider (HLRS, NUM)'
WRITE(std_out,'( A)') '-- Author: Johannes Gebert, M.Sc.  (HLRS, NUM)'
WRITE(std_out,'( A)') '--'
WRITE(std_out,'(2A)') '-- Revision: ', TRIM(ADJUSTL(revision))
WRITE(std_out,'(2A)') '-- Git revision hash: ', TRIM(ADJUSTL(hash))
WRITE(std_out,'( A)') '--'
WRITE(std_out, SEP_STD)
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
SUBROUTINE usage(err)

INTEGER, INTENT(IN) :: err

WRITE(std_out, SEP_STD)
WRITE(std_out, '(A)') 'Directly Discretizing Tensor Computation | Usage:'
WRITE(std_out, SEP_STD)
WRITE(std_out, '(A)') './ddtc_vx.y.z_x86_64 »flags« »basename.meta«'
WRITE(std_out, '(A)') ''
WRITE(std_out, '(A)') '-h/ --help      This message.'
WRITE(std_out, '(A)') '-v/ --version   Version of the program'
WRITE(std_out, '(A)') '--restart       Overwrite restart keyword'
WRITE(std_out, '(A)') '--no-restart    Overwrite restart keyword'
WRITE(std_out, '(A)') ''
WRITE(std_out, '(A)') 'The meta-file must be the last command argument.'
WRITE(std_out, SEP_STD)

CALL handle_err(std_out, '', err)

END SUBROUTINE usage



!------------------------------------------------------------------------------
! SUBROUTINE: date_time
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print the date and time
!
!> @param[in] da Date
!> @param[in] ti Time
!> @param[in] zo Timezone
!> @param[in] long Long or short notation
!> @param[in] str String to feed back
!------------------------------------------------------------------------------  
SUBROUTINE date_time(da, ti, zo, str)

LOGICAL, INTENT(IN) :: da    
LOGICAL, INTENT(IN) :: ti           
LOGICAL, INTENT(IN) :: zo           
CHARACTER(LEN=scl), INTENT(OUT) :: str           

CHARACTER(LEN=8)  :: date
CHARACTER(LEN=10) :: time
CHARACTER(LEN=5)  :: timezone

CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

str = ''
IF(da) str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
IF(ti) str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
IF(zo) str = TRIM(str)//' '//timezone

END SUBROUTINE date_time


!------------------------------------------------------------------------------
! SUBROUTINE: handle_err
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
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
!> 
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!> @param[in] err Errorcode / status of the message
!------------------------------------------------------------------------------  
SUBROUTINE handle_err(fh, txt, err) ! , pro_path, pro_name
 
INTEGER(KIND=ik)             :: fh 
CHARACTER(LEN=*), INTENT(IN) :: txt
INTEGER(KIND=ik)             :: err
! CHARACTER(LEN=*), INTENT(IN) :: pro_path
! CHARACTER(LEN=*), INTENT(IN) :: pro_name

!> Internal variables 
INTEGER(KIND=mpi_ik) :: ierr = 0
CHARACTER(LEN=mcl)   :: text
CHARACTER(LEN=mcl)   :: fmt
CHARACTER(LEN=scl)   :: sub_mssg, errnmbr

! String parsing 
INTEGER  (KIND=ik)   :: ii, jj, sw, mode

CHARACTER(LEN=mcl)   :: delim
CHARACTER(LEN=mcl)   :: tokens(100), path_tokens(50)
CHARACTER(LEN=mcl+1) :: next_token
INTEGER  (KIND=ik)   :: ntokens    ,  path_ntokens

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

        mssg = " with error code "//TRIM(ADJUSTL(errnmbr))//"."

        IF (fh == std_out) THEN
            WRITE(fmt,'(A,I5,A)') "('\x1B[31m','EE Program halted ','\x1B[0m',A, T",&
                scl+LEN('\x1B[31m'//'EE '//'\x1B[0m'),",'\x1B[31m',' EE','\x1B[0m')"
            WRITE(fh, fmt) TRIM(mssg)
        ELSE
            WRITE(fmt,'(A,I5,A)') "('EE Program halted ',A, T",scl + LEN('EE ') + 1,",' EE')"
            WRITE(fh, fmt) TRIM(mssg)
        END IF

        ! CALL stop_slaves(pro_path, pro_name)
    
        CALL MPI_FINALIZE(ierr)
        IF ( ierr /= 0 ) WRITE(fh,'(A)') "MPI_FINALIZE did not succeed"
        STOP 
    END IF

END IF

END SUBROUTINE handle_err


!------------------------------------------------------------------------------
! SUBROUTINE: stop_slaves
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Stop slaves properly
!
!> @param[in] fh Handle of file to print to
!> @param[in] ierr Error code1
!------------------------------------------------------------------------------  
! SUBROUTINE stop_slaves(pro_path, pro_name)

!     CHARACTER(LEN=*), INTENT(IN) :: pro_path
!     CHARACTER(LEN=*), INTENT(IN) :: pro_name
!     Integer(mpi_ik)              :: ierr

!     CALL MPI_BCAST(pro_path, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

!     CALL MPI_BCAST(pro_name, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik,MPI_COMM_WORLD, ierr)
!     !** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
!     CALL MPI_BCAST(-1_ik, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)   

! END SUBROUTINE stop_slaves


!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Automatically writes "sym" if hide_zeros .EQV. .TRUE. and left triangle = 0
!> Hide zeroes is set as default.
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] dim1 Object to print
!> @param[in] dim2 Object to print
!> @param[in] name Name of the object to print
!> @param[in] mat_real Dimensions of the 2nd rank tensor, double precision
!> @param[in] mat_int  Dimensions of the 2nd rank tensor, integer kind = 4
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] hide_zeros Whether to suppress zeros for printing matrices
!------------------------------------------------------------------------------
SUBROUTINE write_matrix (fh, dim1, dim2, name, fmt, unit, mat_real, mat_int, hide_zeros)

INTEGER(KIND=ik)                       , INTENT(IN)           :: fh   
INTEGER(KIND=ik)                       , INTENT(IN)           :: dim1 
INTEGER(KIND=ik)                       , INTENT(IN)           :: dim2 
CHARACTER(LEN=*)                       , INTENT(IN)           :: name 
REAL   (KIND=rk), DIMENSION(:, :)      , INTENT(IN), OPTIONAL :: mat_real    
INTEGER(KIND=ik), DIMENSION(:, :)      , INTENT(IN), OPTIONAL :: mat_int    
CHARACTER(LEN=*)                       , INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*)                       , INTENT(IN), OPTIONAL :: unit 
LOGICAL                                , INTENT(IN), OPTIONAL :: hide_zeros


! Internal variables 
INTEGER(KIND=ik)   :: prec , fw, nm_fmt_lngth, ii, jj
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
LOGICAL            :: hide_zeros_u
LOGICAL            :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
fmt_u = 'standard'
mssg='' 
text = ''

IF((.NOT. PRESENT(mat_real)) .AND. (.NOT. PRESENT(mat_int))) mssg='At least 1 kind of data type required '
IF((      PRESENT(mat_real)) .AND. (      PRESENT(mat_int))) mssg='Please specify max. 1 kind of data type '
IF(mssg /= '') CALL handle_err(fh, TRIM(mssg)//'to print '//TRIM(name)//' matrix', 1)

IF(PRESENT(mat_real)) THEN
    prec = PRECISION(mat_real)
    fw   = prec+8
END IF

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

!------------------------------------------------------------------------------
! Hide zeros by default only if it is a squared matrix
! Check symmetry in case it is a squared matrix.
!------------------------------------------------------------------------------
    hide_zeros_u = .TRUE.

IF (dim1 .EQ. dim2) THEN
    IF(PRESENT(mat_real)) CALL checksym(mat_in =      mat_real          , status=sym_u)
    IF(PRESENT(mat_int )) CALL checksym(mat_in = REAL(mat_int , KIND=rk), status=sym_u)
ELSE
    hide_zeros_u = .FALSE.
    sym_u        = .FALSE.
END IF

! User can explicitly request hiding the zeros 
! even if the matrix is not squared
IF (PRESENT(hide_zeros)) hide_zeros_u = hide_zeros

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')

        IF (PRESENT(mat_real)) WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-2-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')

        IF (PRESENT(mat_real)) WRITE(fmt_a,  "(3(A,I0),A)") "(",dim2,"(F10.3))"
        WRITE(sep  ,  "(A,I0,A)")    "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-2-LEN_TRIM(name)-LEN_TRIM(text) 

   CASE DEFAULT
        mssg='Invalid formatting handle of matrix '//TRIM(name)//'.'
        CALL handle_err(fh=std_out, txt=mssg, err=1)
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(4('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

IF (PRESENT(mat_int )) WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(I10))"

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

IF (hide_zeros_u) THEN                           
    DO ii=1, dim1
    DO jj=1, dim2

        IF ((sym_u) .AND. (ii==dim1-1) .AND. (jj==2)) THEN

            IF(PRESENT(mat_real)) THEN
                IF ((TRIM(fmt_u) .EQ. 'spl') .OR. (TRIM(fmt_u) .EQ. 'simple')) THEN
                    WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
                ELSE
                    WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
                END IF
            ELSE ! Mat int 
                    WRITE(fh, '(A)', ADVANCE='NO') " symmetric"
            END IF

        ELSE
            IF(PRESENT(mat_real)) THEN

                IF ((mat_real (ii,jj) .NE. 0._rk) .AND. &
                ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 

                    WRITE(fh, fmt_a, ADVANCE='NO') mat_real (ii,jj)
                ELSE
                    ! Can'r trim a string with leading and trailing blanks
                    IF ((TRIM(fmt_u) .EQ. 'spl') .OR. (TRIM(fmt_u) .EQ. 'simple')) THEN
                        WRITE(fh, '(A)', ADVANCE='NO') "      .   "
                    ELSE
                        WRITE(fh, '(A)', ADVANCE='NO') "   .                   "
                    END IF
                END IF


            ELSE ! Mat int 
    
                IF ((mat_int (ii,jj) .NE. 0._rk) .AND. &
                ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
                    WRITE(fh, fmt_a, ADVANCE='NO') mat_int (ii,jj)
                ELSE
                    WRITE(fh, '(A)', ADVANCE='NO') "         ."
                END IF
            END IF
        END IF


    END DO
    WRITE(fh,'(A)') ''
    END DO        
ELSE
    IF (PRESENT(mat_real)) WRITE(fh, fmt_a) TRANSPOSE(mat_real)              ! Matrix
    IF (PRESENT(mat_int )) WRITE(fh, fmt_a) TRANSPOSE(mat_int )              ! Matrix
END IF

WRITE(fh, '(A)') ''                               ! Newline & Carriage return

fmt_u   = 'standard'
text = ''

End Subroutine write_matrix


!------------------------------------------------------------------------------
! SUBROUTINE: checksym
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Check the symmetry of an arbitrarily sized square matrix.
!
!> @Description
!> Evaluates the symmetry of the matrix dependent of the quotient of 
!> the L2-Norm of the subtracted minor diagonals devided by the 
!> L2-norm of the input matrix
!
!> @param[in]  mat_in Input Matrix
!> @param[out] mat_out Output Matrix
!> @param[out] status whether the matrix is symmetric
!> @param[in]  llquo Quotient of the L2 Norms
!------------------------------------------------------------------------------  
SUBROUTINE checksym(mat_in, mat_out, status, llquo)
REAL   (KIND=rk), DIMENSION(:,:) , INTENT(IN)              :: mat_in
REAL   (KIND=rk), DIMENSION(:,:) , INTENT(OUT)  , OPTIONAL :: mat_out
LOGICAL                          , INTENT(OUT)  , OPTIONAL :: status
REAL   (KIND=rk)                 , INTENT(OUT)  , OPTIONAL :: llquo

REAL   (KIND=rk), DIMENSION(:,:), ALLOCATABLE :: mat 
REAL   (KIND=rk), DIMENSION(:,:), ALLOCATABLE :: norm_mat
LOGICAL                                       :: status_u

INTEGER(KIND=ik) :: n, m
INTEGER(KIND=ik) :: ii, jj, q
REAL   (KIND=rk) :: norm_norm_mat, norm_in

! m = ABS(m)
n  = SIZE(mat_in, 1)
m  = SIZE(mat_in, 2)

ALLOCATE(mat(n,m))
mat = 0._rk
ALLOCATE(norm_mat(n,m))
norm_mat = 0._rk

norm_in = NORM2(mat_in)

!------------------------------------------------------------------------------
! Calculate the differences of the minor diagonals
!------------------------------------------------------------------------------
q = 1
DO ii=1,n-1 ! columns
    q = q + 1 ! begins with 2
    DO jj=q, m ! rows
        norm_mat (jj,ii) = mat_in (jj,ii) - mat_in (ii,jj) 
    END DO
END DO

!------------------------------------------------------------------------------
! Calculate the L2-Norms and give feedback
!------------------------------------------------------------------------------
norm_norm_mat = NORM2(norm_mat)

! Give feedback about the symmetry of the matrix
IF (norm_norm_mat/norm_in .LT. 10E-06) status_u = .TRUE.

IF(PRESENT(status)) status = status_u

IF(PRESENT(llquo)) llquo = norm_norm_mat/norm_in

q = 1
DO ii=1,n-1 ! columns
    q = q + 1 ! begins with 2
    DO jj=q, m ! rows
        mat (jj,ii) = norm_mat (ii,jj) 
    END DO
END DO

IF(PRESENT(mat_out)) mat_out = mat

END SUBROUTINE checksym

!------------------------------------------------------------------------------
! SUBROUTINE: checksym6x6
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to check the symmetry of a 6x6 matrix. Hardcoded dimensions
!> since it's way quicker than with if/else branches.
!
!> @param[in]  main Input 6x6 Matrix
!> @param[out] maout Input 6x6 Matrix
!> @param[out] status whether the matrix is symmetric
!------------------------------------------------------------------------------ 
SUBROUTINE checksym6x6(fh, main, maout, status)

INTEGER(KIND=ik)             , INTENT(IN)            :: fh   
REAL(KIND=rk), DIMENSION(6,6), INTENT(IN)            :: main
REAL(KIND=rk), DIMENSION(6,6), INTENT(OUT), OPTIONAL :: maout
LOGICAL                      , INTENT(OUT), OPTIONAL :: status

REAL(KIND=rk), DIMENSION(6,6) :: norm_mat
LOGICAL                       :: status_u

! Initialize 
status_u = .FALSE.
norm_mat = 0._rk
CALL write_matrix (fh, 6, 6, 'main', fmt='spl', unit='MPa', mat_real=main)

norm_mat(2:6,1) = [ main(2,1)-main(1,2), main(3,1)-main(1,3), main(4,1)-main(1,4), main(5,1)-main(1,5), main(6,1)-main(1,6) ]
norm_mat(3:6,2) =                      [ main(3,2)-main(2,3), main(4,2)-main(2,4), main(5,2)-main(2,5), main(6,2)-main(2,6) ]
norm_mat(4:6,3) =                                           [ main(4,3)-main(3,4), main(5,3)-main(3,5), main(6,3)-main(3,6) ]
norm_mat(5:6,4) =                                                                [ main(5,4)-main(4,5), main(6,4)-main(4,6) ]
norm_mat(  6,5) =                                                                                       main(6,5)-main(5,6)

! Give feedback about the symmetry of the matrix
IF (NORM2(norm_mat) .LT. 0.01) status_u = .TRUE.
IF(PRESENT(status)) status = status_u

IF(PRESENT(maout))THEN
    maout = main

    maout(2:6,1) = norm_mat(2:6,1)
    maout(3:6,2) = norm_mat(3:6,2)
    maout(4:6,3) = norm_mat(4:6,3)
    maout(5:6,4) = norm_mat(5:6,4)
    maout(  6,5) = norm_mat(  6,5)
END IF


END SUBROUTINE checksym6x6

!------------------------------------------------------------------------------
! SUBROUTINE: zerothres
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets a scalar=0 in case it is less than 10^(-11) by default
!
!> @param[in]    oneD  Scalar input
!> @param[in]    twoD  2-Dim input
!> @param[in]  threeD  3-Dim input
!> @param[in]  threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres(oneD, twoD, threeD, thres)

REAL   (KIND=rk)                  , INTENT(INOUT) , OPTIONAL  ::   oneD 
REAL   (KIND=rk), DIMENSION(:,:)  , INTENT(INOUT) , OPTIONAL  ::   twoD 
REAL   (KIND=rk), DIMENSION(:,:,:), INTENT(INOUT) , OPTIONAL  :: threeD 
REAL   (KIND=rk)                  , INTENT(IN)    , OPTIONAL  :: thres

REAL   (KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii, jj, kk, mm, nn, oo

IF (PRESENT(thres)) THEN
    thres_u = thres
ELSE
    thres_u = 1E-11
END IF

IF(PRESENT(oneD)) THEN
    IF (oneD .GT. 0._rk) THEN
        IF (oneD .LT.  thres_u) oneD = 0._rk
    ELSE
        IF (oneD .GT. -thres_u) oneD = 0._rk
    END IF
END IF

IF(PRESENT(twoD)) THEN
    mm = SIZE(twoD, 1) 
    nn = SIZE(twoD, 2) 

    DO ii=1, mm
    DO jj=1, nn
        IF (twoD(ii,jj) .GT. 0._rk) THEN
            IF ( twoD(ii,jj) .LT.  thres_u) twoD(ii,jj) = 0._rk
        ELSE
            IF ( twoD(ii,jj) .GT. -thres_u) twoD(ii,jj) = 0._rk
        END IF
    END DO
    END DO
END IF

IF(PRESENT(threeD)) THEN
    mm = SIZE(threeD, 1) 
    nn = SIZE(threeD, 2) 
    oo = SIZE(threeD, 3) 

    DO ii=1, mm
    DO jj=1, nn
    DO kk=1, oo
        IF (threeD(ii,jj,kk) .GT. 0._rk) THEN
            IF ( threeD(ii,jj,kk) .LT.  thres_u) threeD(ii,jj,kk) = 0._rk
        ELSE
            IF ( threeD(ii,jj,kk) .GT. -thres_u) threeD(ii,jj,kk) = 0._rk
        END IF
    END DO
    END DO
    END DO
END IF

END SUBROUTINE zerothres

END MODULE auxiliaries
