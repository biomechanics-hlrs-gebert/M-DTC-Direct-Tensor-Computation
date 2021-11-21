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

   Interface write_matrix

      Module Procedure write_matrix_int
      Module Procedure write_matrix_real 
      
   End Interface write_matrix

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
SUBROUTINE usage()

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

END SUBROUTINE usage

!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_real
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
SUBROUTINE write_matrix_real (fh, name, fmt, unit, mat)

INTEGER(KIND=ik), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
REAL   (KIND=rk), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: prec , fw, nm_fmt_lngth, ii, jj, kk, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
mssg='' 
text = ''

prec = PRECISION(mat)
fw = prec+8
sym_u = .FALSE.

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 .EQ. dim2) CALL check_sym(mat_in = mat, status=sym_u)

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')
        WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-2-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')
        WRITE(fmt_a,  "(3(A,I0),A)") "(",dim2,"(F10.3))"
        WRITE(sep  ,  "(A,I0,A)")    "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-2-LEN_TRIM(name)-LEN_TRIM(text) 

   CASE('wxm', 'wxmaxima')
       WRITE(fmt_a, '(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,'],' )"

       WRITE(fh,"(A,A)")TRIM(name),": matrix("

       DO kk = 1, dim1 - 1
          WRITE(fh, fmt_a) mat(kk,:)
       END DO

       WRITE(fmt_a,'(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,']);' )"

       WRITE(fh, fmt_a) mat(dim1, :)
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(4('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, '(A)')
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1-1) .AND. (jj==2)) THEN
        SELECT CASE(fmt_u)
        CASE('spl', 'simple')
            WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
        CASE('std', 'standard')
            WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
        END SELECT
    ELSE
        IF ((ABS(mat(ii,jj)) >=  10E-08) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            SELECT CASE(fmt_u)
            CASE('spl', 'simple')
                WRITE(fh, '(A)', ADVANCE='NO') "      .   "
            CASE('std', 'standard')
                WRITE(fh, '(A)', ADVANCE='NO') "   .                   "
            END SELECT
        END IF
    END IF

END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''                               ! Newline & Carriage return
End Subroutine write_matrix_real


!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_int
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] name Name of the object to print
!> @param[in] mat_real Dimensions of the 2nd rank tensor, double precision
!> @param[in] mat_int  Dimensions of the 2nd rank tensor, integer kind = 4
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] hide_zeros Whether to suppress zeros for printing matrices
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_int (fh, name, fmt, unit, mat)

INTEGER(KIND=ik), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
INTEGER(KIND=ik), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: nm_fmt_lngth, ii, jj, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
sym_u = .FALSE.
mssg='' 
text = ''

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 .EQ. dim2) THEN
    CALL check_sym(fh, mat_in=REAL(mat, KIND=rk), name, sum_sym=sum_sym)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')

        WRITE(sep, "(A,I0,A)") "(",dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formatting to the right
        nm_fmt_lngth = dim2-4-2-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')

        WRITE(sep, "(A,I0,A)") "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-2-LEN_TRIM(name)-LEN_TRIM(text) 
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(4('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(I10))"

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1-1) .AND. (jj==2)) THEN

        WRITE(fh, '(A)', ADVANCE='NO') " symmetric"

    ELSE
        IF ((mat(ii,jj) /= 0) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            WRITE(fh, '(A)', ADVANCE='NO') "         ."
        END IF
    END IF


END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''
End Subroutine write_matrix_int


!------------------------------------------------------------------------------
! SUBROUTINE: check_sym
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
!> @param[in]  fh File handle to write to
!> @param[in]  mat_in Input Matrix
!> @param[in]  name Name of the matrix to evaluate
!> @param[out] mat_out Output Matrix
!> @param[out] sum_sym Sum of differences of all ii,jj entries
!------------------------------------------------------------------------------  
SUBROUTINE check_sym(fh, mat_in, name, mat_out, sum_sym)
INTEGER(KIND=ik),                 INTENT(IN) :: fh
REAL(KIND=rk),    DIMENSION(:,:), INTENT(IN) :: mat_in
CHARACTER(LEN=*),                 INTENT(IN) , OPTIONAL :: name
REAL(KIND=rk)   , DIMENSION(:,:), INTENT(OUT), OPTIONAL :: sum_sym

REAL(KIND=rk), DIMENSION(:,:), ALLOCATABLE :: mat 
REAL(KIND=rk), DIMENSION(:,:), ALLOCATABLE :: sym

INTEGER(KIND=ik) :: n, m
INTEGER(KIND=ik) :: ii, jj, q
REAL   (KIND=rk) :: norm_sym, norm_in
CHARACTER(LEN=scl) :: name_u


name_u = ''
n  = SIZE(mat_in, 1)
m  = SIZE(mat_in, 2)

ALLOCATE(mat(n,m))
mat = 0._rk
ALLOCATE(sym(n,m))
sym = 0._rk

!------------------------------------------------------------------------------
! Calculate the differences of the minor diagonals
!------------------------------------------------------------------------------
q = 1
DO ii=1,n-1 ! columns
    q = q + 1 ! begins with 2
    DO jj=q, m ! rows
        sym (jj,ii) = mat_in (jj,ii) - mat_in (ii,jj) 
    END DO
END DO

q = 1
DO ii=1,n-1 ! columns
    q = q + 1 ! begins with 2
    DO jj=q, m ! rows
        mat (jj,ii) = sym (ii,jj) 
    END DO
END DO

IF(PRESENT(mat_out)) mat_out = mat

IF(PRESENT(name)) name_u = ' '//TRIM(ADJUSTL(name))

WRITE(fh, '(3A, F0.0)')"--------- check_sym",TRIM(name_u),": ", SUM(sym)

END SUBROUTINE check_sym

!------------------------------------------------------------------------------
! SUBROUTINE: check_sym6x6
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
SUBROUTINE check_sym6x6(fh, main, maout, status)

INTEGER(KIND=ik)             , INTENT(IN)            :: fh   
REAL(KIND=rk), DIMENSION(6,6), INTENT(IN)            :: main
REAL(KIND=rk), DIMENSION(6,6), INTENT(OUT), OPTIONAL :: maout
LOGICAL                      , INTENT(OUT), OPTIONAL :: status

REAL(KIND=rk), DIMENSION(6,6) :: norm_mat
LOGICAL                       :: status_u

! Initialize 
status_u = .FALSE.
norm_mat = 0._rk
CALL write_matrix (fh, 'main', 'spl', 'MPa', main)

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


END SUBROUTINE check_sym6x6

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
