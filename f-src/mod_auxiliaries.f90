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
! SUBROUTINE: check_file_exist
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to evaluate whether a file exists. 
!
!>@Description
!>It basically is a wrapper for the subroutine "handle_err".
!> The routine explicitly allows not to abort. Sometimes it simply is not 
!> required and sometimes a more specific error message may be issued.
!
!> stat = 1 --> opposite of target_val
!
!> @param[in] fh File handle to decide where to store the data
!> @param[in] target_val True/False if file shall exist
!> @param[in] filename Filename which is to be checked
!> @param[in] abrt Abort or do not.
!> @param[in] stat Gives Feedback
!------------------------------------------------------------------------------
SUBROUTINE check_file_exist(fh, filename, target_val, abrt, pmssg, stat)

INTEGER(KIND=ik)  , INTENT(IN)              :: fh 
LOGICAL           , INTENT(IN)              :: target_val    
CHARACTER(len=*)  , INTENT(IN)              :: filename      
INTEGER(KIND=ik)  , INTENT(IN)   , OPTIONAL :: abrt
LOGICAL           , INTENT(IN)   , OPTIONAL :: pmssg
INTEGER(KIND=ik)  , INTENT(OUT)  , OPTIONAL :: stat

!-- Internal Variable
LOGICAL                                     :: exist=.FALSE. 
INTEGER(KIND=ik)                            :: abrt_u=1
LOGICAL                                     :: pmssg_u=.TRUE.

! Initialize
IF(PRESENT(abrt))  abrt_u  = abrt
IF(PRESENT(pmssg)) pmssg_u = pmssg
IF(PRESENT(stat))  stat    = 0

INQUIRE (FILE = TRIM(filename), EXIST = exist)

! Since exist can take both states, the meaning of "stat" in context automatically follows.
IF (exist .NEQV. target_val)  THEN
    IF (target_val .EQV. .FALSE.) THEN
        mssg='The file '//TRIM(filename)//' already exists.'
        IF(PRESENT(stat))stat = 1
    ELSE
        mssg='The file '//TRIM(filename)//' does not exist.'
        IF(PRESENT(stat))stat = 1
    END IF

    IF(pmssg_u .EQV. .FALSE.) mssg = ''
    ! Only raise an error if the program shall stop.
    CALL handle_err(fh=fh, txt=TRIM(mssg), err=abrt_u)
END IF
    
END SUBROUTINE check_file_exist


!------------------------------------------------------------------------------
! SUBROUTINE: check_and_open
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to check the existance of a file and to eventually close it.
!
!> @param[in] fh File handle to decide where to store the data
!> @param[in] filename Filename which is to be checked
!> @param[in] restart Whether to restart or not to.
!> @param[in] abrt Abort or do not.
!> @param[in] stat Gives Feedback
!------------------------------------------------------------------------------  
SUBROUTINE check_and_open(fh, filename, restart, abrt, stat)

INTEGER(KIND=ik)  , INTENT(IN)              :: fh 
CHARACTER(len=*)  , INTENT(IN)              :: filename  
LOGICAL           , INTENT(IN)   , OPTIONAL :: restart
INTEGER(KIND=ik)  , INTENT(IN)   , OPTIONAL :: abrt
INTEGER(KIND=ik)  , INTENT(OUT)  , OPTIONAL :: stat

INTEGER(KIND=ik)                            :: abrt_u=1

!-- Internal Variable
LOGICAL                                     :: opened=.FALSE.
LOGICAL                                     :: restart_u=.FALSE.
CHARACTER(len=mcl)                          :: filename_u=''

! In case of doubt, abort.
IF(PRESENT(restart)) restart_u=restart
IF(PRESENT(abrt)) abrt_u=abrt
IF(PRESENT(stat)) stat=0

CALL check_file_exist(fh, filename_u, .FALSE., abrt_u, .TRUE., stat)

INQUIRE(UNIT=fh, OPENED=opened)

IF (opened .EQV. .FALSE.) THEN
   !------------------------------------------------------------------------------
   ! Open the file
   !------------------------------------------------------------------------------
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='READWRITE', ACCESS='SEQUENTIAL', STATUS='OLD')
   IF(PRESENT(stat)) stat=0_ik
END IF
 
END SUBROUTINE check_and_open

  
!------------------------------------------------------------------------------
! SUBROUTINE: check_file_exist
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to check the existance of a file and to eventually close it.
!
!> @param[in] fh File handle to decide where to store the data
!> @param[in] filename Filename which is to be checked
!> @param[in] abrt Abort or do not.
!> @param[in] stat Gives Feedback
!------------------------------------------------------------------------------  
SUBROUTINE check_and_close(fh, filename, abrt, stat)

INTEGER(KIND=ik)  , INTENT(IN)              :: fh 
CHARACTER(len=*)  , INTENT(IN)   , OPTIONAL :: filename      
LOGICAL           , INTENT(IN)   , OPTIONAL :: abrt
INTEGER(KIND=ik)  , INTENT(OUT)  , OPTIONAL :: stat

LOGICAL                                     :: abrt_u=.TRUE.
INTEGER(KIND=ik)                            :: stat_u

!-- Internal Variable
LOGICAL                                     :: opened=.FALSE. 
CHARACTER(len=mcl)                          :: filename_u

! Initialize
stat_u = 0
filename_u=''

! In case of doubt, abort.
IF(PRESENT(abrt)) abrt_u = abrt
IF(PRESENT(stat)) stat_u = stat
IF(PRESENT(filename)) filename_u = TRIM(filename)


INQUIRE(UNIT=fh, OPENED=opened)

IF (opened .EQV. .TRUE.) THEN
   CLOSE (fh)
ELSE
    mssg='The file »'//TRIM(filename_u)//'« was closed already.'

    ! Whether to stop the program has to be decided via the call.
    IF (abrt_u .EQV. .TRUE.) stat_u = 1
    CALL handle_err(fh, TRIM(mssg), stat_u)
END IF
 
END SUBROUTINE check_and_close


!------------------------------------------------------------------------------
! SUBROUTINE: date_time
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print the date and time
!
!> @param[in] fh File handle to write to
!> @param[in] da Date
!> @param[in] ti Time
!> @param[in] zo Timezone
!> @param[in] long Long or short notation
!> @param[in] mssgdt Message to print before
!------------------------------------------------------------------------------  
SUBROUTINE date_time(fh, da, ti, zo, mssgdt)

    INTEGER(KIND=ik)  , INTENT(IN), OPTIONAL :: fh 
    LOGICAL           , INTENT(IN), OPTIONAL :: da    
    LOGICAL           , INTENT(IN), OPTIONAL :: ti           
    LOGICAL           , INTENT(IN), OPTIONAL :: zo           
    CHARACTER(LEN=*)  , INTENT(IN), OPTIONAL :: mssgdt           

    CHARACTER(LEN=8)                         :: date
    CHARACTER(LEN=10)                        :: time
    CHARACTER(LEN=5)                         :: timezone

CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

IF(PRESENT(mssgdt)) WRITE(fh, "('MM ',A)", ADVANCE='NO') TRIM(mssgdt)

IF(da) WRITE(fh, FMT_DA, ADVANCE='NO') date(7:8), date(5:6), date(1:4)
IF(ti) WRITE(fh, FMT_TI, ADVANCE='NO') time(1:2), time(3:4), time(5:10)
IF(zo) WRITE(fh, FMT_ZO) timezone

END SUBROUTINE date_time

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
!> @param[in] pro_path Project Path
!> @param[in] pro_name Project Name
!> @param[in] ierr Error code1
!------------------------------------------------------------------------------  
SUBROUTINE stop_slaves(pro_path, pro_name)
 
    CHARACTER(LEN=*), INTENT(IN) :: pro_path
    CHARACTER(LEN=*), INTENT(IN) :: pro_name
    Integer(mpi_ik)              :: ierr

    CALL MPI_BCAST(pro_path, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik, MPI_COMM_WORLD, ierr)

    CALL MPI_BCAST(pro_name, INT(mcl,mpi_ik), MPI_CHAR, 0_mpi_ik,MPI_COMM_WORLD, ierr)
    !** Bcast Serial_root_size = -1 ==> Signal for slave to stop ***
    CALL MPI_BCAST(-1_ik, 1_mpi_ik, MPI_INTEGER8, 0_mpi_ik, MPI_COMM_WORLD, ierr)   

END SUBROUTINE stop_slaves

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
SUBROUTINE handle_err(fh, txt, err)
 
INTEGER  (KIND=ik)                        :: fh 
CHARACTER(LEN=*)  , INTENT(IN)            :: txt
INTEGER  (KIND=ik)                        :: err

!> Internal variables 
INTEGER(KIND=mpi_ik)                      :: ierr = 0
CHARACTER(LEN=mcl)                        :: text
CHARACTER(LEN=mcl)                        :: fmt
CHARACTER(LEN=scl)                        :: sub_mssg, errnmbr

! String parsing 
INTEGER  (KIND=ik)                        :: ii, jj, sw, mode

CHARACTER(LEN=mcl)                        :: delim
CHARACTER(LEN=mcl)                        :: tokens(100), path_tokens(50)
CHARACTER(LEN=mcl+1)                      :: next_token
INTEGER  (KIND=ik)                        :: ntokens    ,  path_ntokens

text = TRIM(ADJUSTL(txt))
fmt  = FMT_ERR
     
mode = 0                ! Absolute or relative path
sw=2
ntokens         = 0
path_ntokens    = 0
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
    
               


        CALL stop_slaves(out%path, out%bsnm)

        CALL MPI_FINALIZE(ierr)
        IF ( ierr /= 0 ) WRITE(fh,'(A)') "MPI_FINALIZE did not succeed"
        STOP 
    END IF

END IF

END SUBROUTINE handle_err

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
WRITE(std_out, '(A)') './ddtc_vx.y.z_x86_64 »flags« »meta filename«'
WRITE(std_out, '(A)') '--restart       Overwrite restart keyword'
WRITE(std_out, '(A)') '--no-restart    Overwrite restart keyword'
WRITE(std_out, '(A)') '-h              This message.'
WRITE(std_out, '(A)') '*.meta          Meta input file.'
WRITE(std_out, '(A)') ''
WRITE(std_out, '(A)') 'The meta-file must be the last command argument.'
WRITE(std_out, SEP_STD)

CALL handle_err(std_out, '', err)

END SUBROUTINE usage


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
INTEGER(KIND=ik)                                              :: prec , fw, nm_fmt_lngth, ii, jj
CHARACTER(LEN=mcl)                                            :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl)                                            :: text
CHARACTER(LEN=mcl)                                            :: fmt_u
LOGICAL                                                       :: hide_zeros_u
LOGICAL                                                       :: sym_u

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

IF (hide_zeros_u .EQV. .TRUE.) THEN                           
    DO ii=1, dim1
    DO jj=1, dim2

        IF ((sym_u .EQV. .TRUE.) .AND. (ii==dim1-1) .AND. (jj==2)) THEN

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
                ((sym_u .EQV. .FALSE.) .OR. ((sym_u .EQV. .TRUE.) .AND. (jj .GE. ii)))) THEN 

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
                ((sym_u .EQV. .FALSE.) .OR. ((sym_u .EQV. .TRUE.) .AND. (jj .GE. ii)))) THEN 
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
REAL   (KIND=rk), DIMENSION(:,:)             , INTENT(IN)              :: mat_in
REAL   (KIND=rk), DIMENSION(:,:)             , INTENT(OUT)  , OPTIONAL :: mat_out
LOGICAL                                      , INTENT(OUT)  , OPTIONAL :: status
REAL   (KIND=rk)                             , INTENT(OUT)  , OPTIONAL :: llquo

REAL   (KIND=rk), DIMENSION(:,:), ALLOCATABLE                          :: mat 
REAL   (KIND=rk), DIMENSION(:,:), ALLOCATABLE                          :: norm_mat
LOGICAL                                                                :: status_u

INTEGER(KIND=ik)                                                       :: n, m
INTEGER(KIND=ik)                                                       :: ii, jj, q
REAL   (KIND=rk)                                                       :: norm_norm_mat, norm_in

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
SUBROUTINE checksym6x6(main, maout, status)

REAL    (KIND=rk), DIMENSION(6,6), INTENT(IN)            :: main
REAL    (KIND=rk), DIMENSION(6,6), INTENT(OUT), OPTIONAL :: maout
LOGICAL                          , INTENT(OUT), OPTIONAL :: status

REAL    (KIND=rk), DIMENSION(6,6)                        :: norm_mat
LOGICAL                                                  :: status_u

! Initialize 
status_u = .FALSE.
norm_mat = 0._rk
CALL write_matrix (fhr, 6, 6, 'main', fmt='spl', unit='MPa', mat_real=main)

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

REAL   (KIND=rk)                                              :: thres_u
INTEGER(KIND=ik)                                              :: ii, jj, kk, mm, nn, oo

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
