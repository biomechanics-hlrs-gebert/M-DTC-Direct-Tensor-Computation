!------------------------------------------------------------------------------
! HLRS - Numerical Methods and Libraries, Nobelstraße 19, 70569 Stumgart
!------------------------------------------------------------------------------
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
! USE MPI

IMPLICIT NONE
   Interface write_matrix
      Module Procedure write_matrix_int
      Module Procedure write_matrix_real 
   End Interface write_matrix

   Interface zero_thres
      Module Procedure zerothres_num
      Module Procedure zerothres_OnD
      Module Procedure zerothres_TwD
      Module Procedure zerothres_ThD
   End Interface zero_thres

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
REAL(KIND=rk) :: sym_out
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

IF (dim1 .EQ. dim2) THEN
    CALL check_sym(fh, mat, name, sym_out=sym_out)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')
        WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-LEN_TRIM(name)-LEN_TRIM(text)

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
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    



!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, '(A)')
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        SELECT CASE(fmt_u)
        CASE('spl', 'simple')
            WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
        CASE('std', 'standard')
            WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
        END SELECT

    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
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
REAL(KIND=rk) :: sym_out
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
    CALL check_sym(fh, REAL(mat, KIND=rk), name, sym_out=sym_out)
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
        nm_fmt_lngth  = dim2*10-4-LEN_TRIM(name)-LEN_TRIM(text) 
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(I10))"

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        WRITE(fh, '(A)', ADVANCE='NO') " symmetric"
    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
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
!> @param[in]  mi Input Matrix
!> @param[in]  name Name of the matrix to evaluate
!> @param[out] mo Output Matrix
!> @param[out] sym_out Sum of differences of all ii,jj entries
!------------------------------------------------------------------------------  
SUBROUTINE check_sym(fh, mi, name, mo, sym_out)

INTEGER(KIND=ik),                 INTENT(IN) :: fh
REAL(KIND=rk),    DIMENSION(:,:), INTENT(IN) :: mi
CHARACTER(LEN=*),                 INTENT(IN) , OPTIONAL :: name
REAL(KIND=rk)   , DIMENSION(:,:), INTENT(OUT), OPTIONAL :: mo
REAL(KIND=rk)                   , INTENT(OUT), OPTIONAL :: sym_out

REAL(KIND=rk), DIMENSION(:,:), ALLOCATABLE :: mat 
REAL(KIND=rk) :: sym 
INTEGER, DIMENSION(2) :: lb, ub
INTEGER(KIND=ik) :: ii, jj, q
CHARACTER(LEN=scl) :: name_u, sym_out_str

name_u = ''

lb = LBOUND(mi)
ub = UBOUND(mi)

ALLOCATE(mat(lb(1), ub(2)))
mat = 0._rk

!------------------------------------------------------------------------------
! Calculate the differences to get the information of symmetry
!------------------------------------------------------------------------------
sym = 0._rk

Do jj = lb(2), ub(2)
    DO ii = lb(1), ub(1)
        sym = sym + ABS(mi(ii,jj) -  mi(jj,ii))
    End DO
End Do

IF(PRESENT(sym_out)) sym_out = sym

!------------------------------------------------------------------------------
! Write matrix out with zeros to show check_sym
!------------------------------------------------------------------------------
mat = mi
IF(sym <= 10E-06) THEN
    q = 1
    DO jj=lb(2), ub(2)-1 ! columns
        q = q + 1 ! begins with 2
        DO ii=lb(1), ub(2) ! rows
            mat (ii,jj) = 0._rk
        END DO
    END DO
END IF

IF(PRESENT(mo)) mo = mat
DEALLOCATE(mat)

IF(PRESENT(name)) name_u = ' '//TRIM(ADJUSTL(name))

!------------------------------------------------------------------------------
! Format output
!------------------------------------------------------------------------------
WRITE(sym_out_str, '(F40.20)') sym

CALL trimzero(sym_out_str)

WRITE(fh, '(4A)')"-- check_sym",TRIM(name_u),": ", TRIM(ADJUSTL(sym_out_str))

END SUBROUTINE check_sym

!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_num
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets a scalar=0 in case it is less than 10^(-11) by default
!
!> @param[in] num Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_num(num, thres)

REAL(KIND=rk), INTENT(INOUT) :: num 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

IF (num >= 0._rk) THEN
    IF (num <=  thres_u) num = 0._rk
ELSE
    IF (num >= -thres_u) num = 0._rk
END IF

END SUBROUTINE zerothres_num


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_OnD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets a vector=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_OnD(oneD, thres)

REAL(KIND=rk), DIMENSION(:), INTENT(INOUT) :: oneD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO ii=1, SIZE(oneD)
    IF (oneD(ii) >= 0._rk) THEN
        IF (oneD(ii) <=  thres_u) oneD(ii) = 0._rk
    ELSE
        IF (oneD(ii) >= -thres_u) oneD(ii) = 0._rk
    END IF
END DO
END SUBROUTINE zerothres_OnD


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_TwD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets an array=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_TwD(TwD, thres)

REAL(KIND=rk), DIMENSION(:,:), INTENT(INOUT) :: TwD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii, jj

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO jj=1, SIZE(TwD, 2)
DO ii=1, SIZE(TwD, 1)
    IF (TwD(ii, jj) >= 0._rk) THEN
        IF (TwD(ii, jj) <=  thres_u) TwD(ii, jj) = 0._rk
    ELSE
        IF (TwD(ii, jj) >= -thres_u) TwD(ii, jj) = 0._rk
    END IF
END DO
END DO
END SUBROUTINE zerothres_TwD


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_ThD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets an array=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_ThD(ThD, thres)

REAL(KIND=rk), DIMENSION(:, :, :), INTENT(INOUT) :: ThD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii, jj, kk

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO kk=1, SIZE(ThD, 3)
DO jj=1, SIZE(ThD, 2)
DO ii=1, SIZE(ThD, 1)
    IF (ThD(ii, jj, kk) >= 0._rk) THEN
        IF (ThD(ii, jj, kk) <=  thres_u) ThD(ii, jj, kk) = 0._rk
    ELSE
        IF (ThD(ii, jj, kk) >= -thres_u) ThD(ii, jj, kk) = 0._rk
    END IF
END DO
END DO
END DO
END SUBROUTINE zerothres_ThD

!------------------------------------------------------------------------------
! SUBROUTINE: calc_real_time
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Returns the elapsed real time in hh:mm:ss.sss format
!
!> @description
!> The subroutines calculates the elapsed real time from two parameter sets
!> returned by the intinsic date_and_time(values)
!
!> @param[in] tstart Measured start time
!> @param[in] tend Measured end time
!> @param[in] elapsed Time elapsed
!> @param[in] echo Print elapsed time to a file handle
!------------------------------------------------------------------------------ 
subroutine calc_real_time(fh, tstart, tend, elapsed, echo)

Integer              , intent(in)  :: fh
Integer, Dimension(8), intent(in)  :: tstart
Integer, Dimension(8), intent(in)  :: tend
Character(len=12)    , intent(out) :: elapsed
Logical, optional    , intent(in)  :: echo

Integer                           :: hh, mm, ss, msec
Integer(Kind=8)                   :: msec_s, msec_e, msec_diff

msec_s = 0
msec_s = msec_s + tstart(8) + tstart(7) * 1000 + tstart(6) * 60 * 1000
msec_s = msec_s + tstart(5) * 60 * 60 * 1000

msec_e = 0
msec_e = msec_e + tend(8) + tend(7) * 1000 + tend(6) * 60 * 1000
msec_e = msec_e + tend(5) * 60 * 60 * 1000

if (msec_e < msec_s) msec_e = msec_e + 24*60*60*1000

msec_diff = msec_e - msec_s

hh = msec_diff/(60 * 60 * 1000)
mm = (msec_diff - hh * (60 * 60 * 1000)) / (60 * 1000)
ss = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000) / 1000
msec = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000 - ss * 1000)

write(elapsed,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)") hh,mm,ss,msec         

If (present(echo)) then
    If (echo) Write(fh, FMT_TIME)&
        'Elapsed time was: '//elapsed//' : ',Real(msec_diff)/1000.
End If

End subroutine calc_real_time

END MODULE auxiliaries
