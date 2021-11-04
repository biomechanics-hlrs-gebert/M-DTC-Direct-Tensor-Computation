




!------------------------------------------------------------------------------
! MODULE: meta
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! DESCRIPTION:
!> Module containing all meta file read/write routines.
!
! REVISION HISTORY:
! 21 10 2021 - Initial refactored version
!------------------------------------------------------------------------------
MODULE meta

USE global_std
USE auxiliaries
USE strings

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: meta_handle_lock_file
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to encapsule the lock file handling
!
!> @param[in] restart Whether to restart or not to.
!---------------------------------------------------------------------------  
SUBROUTINE meta_handle_lock_file(restart)

CHARACTER         , INTENT(IN)          :: restart

LOGICAL                                 :: exist=.FALSE.
INTEGER  (KIND=ik)                      :: ios
CHARACTER(LEN=mcl)                      :: lockname

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
lockname=TRIM(in%path)//'.'//TRIM(in%bsnm)//lock_suf

INQUIRE (FILE = TRIM(lockname), EXIST = exist)

IF((restart .EQ. 'N') .AND. (exist .EQV. .TRUE.)) THEN
   mssg='The .*.lock file is set and a restart prohibited by default or the user.'
   CALL handle_err(std_out, TRIM(ADJUSTL(mssg)), err=1_ik)
END IF

IF(((restart .EQ. 'Y') .AND. (exist .EQV. .FALSE.)) .OR. ((restart .EQ. 'N') .AND. (exist .EQV. .FALSE.))) THEN
   CALL execute_command_line ('touch '//TRIM(lockname), CMDSTAT=ios)
   CALL handle_err(std_out, 'The .*.lock file could not be set.', err=ios)
END IF

IF((restart .EQ. 'Y') .AND. (exist .EQV. .TRUE.)) CONTINUE

END SUBROUTINE meta_handle_lock_file



!------------------------------------------------------------------------------
! SUBROUTINE: meta_append
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to open a meta file to append data/ keywords
!
!> @param[inout] meta_as_rry Meta data written into a character array
!---------------------------------------------------------------------------  
SUBROUTINE meta_append(meta_as_rry)

CHARACTER(LEN=mcl), DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: meta_as_rry      

! Internal Variables
CHARACTER(LEN=mcl)                                           :: line
INTEGER  (KIND=ik)                                           :: ios, lines, ii
CHARACTER(LEN=mcl)                                           :: tokens(30)
INTEGER  (KIND=ik)                                           :: ntokens

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
CALL check_file_exist(std_out, filename=in%full, target_val=.TRUE., abrt=1, stat=ios)

CALL parse( str=in%full, delims=".", args=tokens, nargs=ntokens)

IF ( '.'//TRIM(tokens(ntokens)) .EQ. meta_suf) THEN
   !------------------------------------------------------------------------------
   ! Parse all basename and path details.
   !------------------------------------------------------------------------------
   in%p_n_bsnm = in%full(1:LEN_TRIM(in%full)-LEN_TRIM(meta_suf)) 
   
   CALL parse( str=TRIM(in%p_n_bsnm), delims="/", args=tokens, nargs=ntokens)

   in%path      = in%p_n_bsnm(1:LEN_TRIM(in%p_n_bsnm) - LEN_TRIM(tokens(ntokens)))     
   in%bsnm      = TRIM(tokens(ntokens))
   in%bsnm      = TRIM(tokens(ntokens))

   CALL parse( str=TRIM(in%bsnm), delims="_", args=tokens, nargs=ntokens)

   in%dataset   = TRIM(tokens(1))
   in%type      = TRIM(tokens(2))
   in%purpose   = TRIM(tokens(3))
   in%app       = TRIM(tokens(4))
   in%features  = TRIM(tokens(5))

   out = in  
ELSE
   ! File is not a meta file
   CALL handle_err(std_out, "The input file is not a *"//meta_suf//" file.", 1)
END IF

!------------------------------------------------------------------------------
! Open the meta input file
!------------------------------------------------------------------------------
OPEN(UNIT=fhmei, FILE=TRIM(in%full), ACTION='READWRITE', ACCESS='SEQUENTIAL', STATUS='OLD')

! Count lines of input file
lines = 0_ik
DO
   READ(fhmei, '(A)', iostat=ios) line
   IF ( ios .NE. 0_ik ) EXIT
   
   lines = lines + 1_ik
END DO

! Reposition to first line of file
REWIND (fhmei)

ALLOCATE(meta_as_rry(lines))
!------------------------------------------------------------------------------
! Read all lines into the file
!------------------------------------------------------------------------------
DO ii=1, lines
   READ(fhmei,'(A)') meta_as_rry(ii)
END DO

!------------------------------------------------------------------------------
! Parse and check basename
!------------------------------------------------------------------------------
CALL parse(str=TRIM(in%bsnm), delims='_', args=tokens, nargs=ntokens)

! Check if the basename consists of exactly the 5 parts.
IF(ntokens /= 5_ik) THEN   
   mssg='The basename »'//TRIM(in%bsnm)//'« of the meta-file was ill-defined. It may be parsed wrong.'
   CALL handle_err(std_out, TRIM(ADJUSTL(mssg)), 0)
END IF

!------------------------------------------------------------------------------
! Alter the meta file name
! The variable »alter« must be given and must be true, 
! because its a dangerous operation which may lead to data loss.
!------------------------------------------------------------------------------
CALL meta_io (fhmon, 'NEW_BSNM_FEATURE', '', meta_as_rry, chars= out%features, kwabrt=0, stat=ios)
IF((ios ==1)) out%features = in%features           ! if ios = 1 --> no keyword found

CALL meta_io (fhmon, 'NEW_BSNM_PURPOSE', '', meta_as_rry, chars= out%purpose  , kwabrt=0, stat=ios)
IF((ios ==1)) out%purpose = in%purpose             ! if ios = 1 --> no keyword found


IF ((out%purpose == in%purpose) .AND. (out%features == in%features)) THEN
   mssg='The basename did not change. When in doubt, please check your meta file.'
   CALL handle_err(std_out, mssg, 0)
END IF


!------------------------------------------------------------------------------
! Build the new outfile path
!------------------------------------------------------------------------------
! Nomenclature: dataset_type_purpose_app_features
! This assignment requres the out = in assignment before
out%bsnm =     TRIM(out%dataset)//&
          '_'//TRIM(out%type)//&
          '_'//TRIM(out%purpose)//&
          '_'//TRIM(out%app)//&        ! out%app shall be defined in the main program!
          '_'//TRIM(out%features)

out%p_n_bsnm = TRIM(out%path)//&
               TRIM(out%bsnm)

out%full = TRIM(out%p_n_bsnm)//meta_suf

!------------------------------------------------------------------------------
! System call to update the file name of the meta file
!------------------------------------------------------------------------------
CALL execute_command_line ('cp '//TRIM(in%full)//' '//TRIM(out%full), CMDSTAT=ios)
CALL handle_err(std_out, 'The update of the meta filename went wrong.', ios)

!------------------------------------------------------------------------------
! Open the meta output file
!------------------------------------------------------------------------------
OPEN(UNIT=fhmeo, FILE=TRIM(out%full), ACTION='READWRITE', ACCESS='SEQUENTIAL', STATUS='OLD')

WRITE(fhmeo, '(A)')
WRITE(fhmeo, FMT_HY_SEP)
WRITE(fhmeo, '(A)')

END SUBROUTINE meta_append




!------------------------------------------------------------------------------
! SUBROUTINE: meta_add_ascii
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to deal with the logging and renaming of the log file in context 
!> of the meta data approach. 
!
!> @Description
!> Meta log only gets called __after__ meta append or meta_close respectively. 
!> This way, a log file is optional.
!
!> If Restart isn's set or .FALSE., then th program aborts in case a meta or
!> a temporary log file is present as they might contain important data or
!> debugging information. The temporaray logfile is a hidden one.
!>
!> If the variable restart is not explicitly set .TRUE., the program will not 
!> restart.
!> If the variable in/out are not set, the program will not start/stop 
!> accordingly.
!
!> The Name of the temporary logfile is hardcoded!
!> »'.temporary.'//log_suf«
!
!> @param[in] fh File handle of the input
!> @param[in] suf Suffix of the file
!> @param[in] st Metafile of the input
!> @param[in] restart Logfiles (temporary and permanent)
!---------------------------------------------------------------------------  
SUBROUTINE meta_add_ascii(fh, suf, st, restart)

INTEGER  (KIND=ik)              , INTENT(IN)              :: fh
CHARACTER(LEN=*)                , INTENT(IN)              :: suf
CHARACTER(LEN=*)                , INTENT(IN)              :: st
CHARACTER                       , INTENT(IN)   , OPTIONAL :: restart

CHARACTER(LEN=mcl)                                        :: temp_f_suf, perm_f_suf
INTEGER  (KIND=ik)                                        :: ios, stat_t, stat_p
CHARACTER                                                 :: restart_u='N'


! The temporaray file is a hidden one.
temp_f_suf = TRIM(out%path)//'.temporary'//TRIM(suf)
perm_f_suf = TRIM(out%p_n_bsnm)//TRIM(suf)

!------------------------------------------------------------------------------
! Create the file
!------------------------------------------------------------------------------
IF (st == 'start') THEN

   !------------------------------------------------------------------------------
   ! Check restart prerequisites
   !------------------------------------------------------------------------------
   ! Check for restart default value
   IF(PRESENT(restart)) restart_u=restart

   ! Check for temporary file
   CALL check_file_exist(std_out, filename=temp_f_suf, target_val=.FALSE., pmssg=.FALSE., abrt=0, stat=stat_t)

   ! Check for a permanent file
   CALL check_file_exist(std_out, filename=out%p_n_bsnm//TRIM(suf), target_val=.FALSE., pmssg=.FALSE., abrt=0, stat=stat_p)

   !------------------------------------------------------------------------------
   ! What happens when a restart is requested.
   !------------------------------------------------------------------------------
   IF (restart_u .EQ. 'Y') THEN
      ! if target_val if check_file_exist = .FALSE. and stat_*l = 0 - the file does not exist
      IF(stat_t == 1) THEN
         CALL execute_command_line ('rm -r '//TRIM(temp_f_suf), CMDSTAT=ios)   
         CALL handle_err(std_out, '»'//TRIM(temp_f_suf)//'« not deletable.',ios)
      END IF

      IF(stat_p == 1) THEN
         CALL execute_command_line ('rm -r '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)
         CALL handle_err(std_out, '»'//TRIM(out%full)//'« not deletable.', ios)
      END IF

   ELSE ! restart_u .EQ. 'N'
      !------------------------------------------------------------------------------
      ! If no restart is requested (default)
      !------------------------------------------------------------------------------
      IF    ((stat_t /= 0) .OR. (stat_p /= 0)) THEN
         IF ((stat_t == 1) .AND. (stat_p == 1)) THEN 
            mssg='The file '//TRIM(perm_f_suf)//' and the file '//TRIM(temp_f_suf)//' already exists.'
         ELSE IF  (stat_t == 1) THEN
            mssg='The file '//TRIM(temp_f_suf)//' already exists.'
         ELSE ! (stat_p == 1) 
            mssg='The file '//TRIM(perm_f_suf)//' already exists.'
         END IF

         CALL handle_err(std_out, mssg, 1)     
      END IF
   END IF

   OPEN(UNIT=fh, FILE=TRIM(temp_f_suf), ACTION='WRITE', ACCESS='SEQUENTIAL', STATUS='NEW')

END IF !  (st == 'start') THEN

!------------------------------------------------------------------------------
! Close the file
!------------------------------------------------------------------------------
IF (TRIM(st) == 'stop') THEN
   CLOSE (fh)

   ! The temporary log file must be renamed to a permanent one
   CALL execute_command_line ('mv '//TRIM(temp_f_suf)//' '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)

   IF(ios /= 0_ik) THEN
      mssg='Can not rename the suffix_file from »'//TRIM(temp_f_suf)//'« to the proper basename.'
      CALL handle_err(std_out, mssg, 0)
   END IF
END IF

END SUBROUTINE meta_add_ascii



!------------------------------------------------------------------------------
! SUBROUTINE: meta_io
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse and to write information in context of a keyword. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
!> The presence of m_in determines whether to read or write.
!> Depending on the datatype given, the routine reads and writes with the
!> appropriate type. The variable data types are not given as INTEN(INOUT)
!> since unspecified - and therefore allocated - variables in context of
!> the call lead to an error. 
!> Please ensure, the datatype and the content of a meta input file match.
!> The routine will not check against it and the program may crash.
!> 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] unit Unit of the value in case a text output is requested.
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars    Optional datatype to read in
!> @param[in] int__0D  Optional datatype to read in
!> @param[in] real_0D  Optional datatype to read in
!> @param[in] int__1D2 Optional datatype to read in
!> @param[in] real_1D2 Optional datatype to read in
!> @param[in] int__1D3 Optional datatype to read in
!> @param[in] real_1D3 Optional datatype to read in
!> @param[in] stat Status
!> @param[in] kwabrt Default -> .TRUE.; do abort
!> @param[in] nd Suppress the time and date of the output
!---------------------------------------------------------------------------
SUBROUTINE meta_io (fh, keyword, unit, m_in, chars, &
                                              int_0D,  &
                                             real_0D,  & 
                                              int_1D2, & 
                                             real_1D2, & 
                                              int_1D3, & 
                                             real_1D3, &
                                             kwabrt, stat, nd)
   
INTEGER  (KIND=ik)              , INTENT(IN)              :: fh 
CHARACTER(LEN=*)                                          :: keyword
CHARACTER(LEN=*)                               , OPTIONAL :: unit
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(INOUT), OPTIONAL :: m_in      
CHARACTER(LEN=*)                               , OPTIONAL :: chars 
INTEGER  (KIND=ik)                             , OPTIONAL ::  int_0D 
REAL     (KIND=rk)                             , OPTIONAL :: real_0D 
INTEGER  (KIND=ik), DIMENSION(2)               , OPTIONAL ::  int_1D2
REAL     (KIND=rk), DIMENSION(2)               , OPTIONAL :: real_1D2
INTEGER  (KIND=ik), DIMENSION(3)               , OPTIONAL ::  int_1D3 
REAL     (KIND=rk), DIMENSION(3)               , OPTIONAL :: real_1D3
INTEGER  (KIND=ik),               INTENT(IN)   , OPTIONAL :: kwabrt
INTEGER  (KIND=ik),               INTENT(INOUT), OPTIONAL :: stat
LOGICAL           ,               INTENT(IN)   , OPTIONAL :: nd 

! Internal variables
CHARACTER(LEN=kcl)                                        :: kywd_lngth
CHARACTER(LEN=ucl)                                        :: unit_lngth
CHARACTER(LEN=mcl)                                        :: tokens(30), line
INTEGER  (KIND=ik)                                        :: ntokens, datatype, ii, do_loop_counter
INTEGER  (KIND=ik)                                        :: cntr, kywd_found, maxchars
INTEGER  (KIND=ik)                                        :: unit_post_fill, kwabrt_u
LOGICAL                                                   :: ndu

!------------------------------------------------------------------------------
! Initialize variables
! Should be defined from new each time the program enters the routine. 
! Otherwise some errors may occur.
!------------------------------------------------------------------------------
cntr           = 0
kywd_found     = 0
kwabrt_u = 1
ndu      = .FALSE.

IF(PRESENT(kwabrt)) kwabrt_u = kwabrt
IF(PRESENT(stat))   stat     = 0
IF(PRESENT(nd))     ndu      = nd

!------------------------------------------------------------------------------
! Check unit keyword for convention and proper formatting
!------------------------------------------------------------------------------
IF(LEN_TRIM(keyword) .GT. LEN(kywd_lngth)) THEN
   WRITE(fh, '(A)') ''
   mssg = "The keyword »"//TRIM(keyword)//"« is longer than the convention allows and therefore truncated!"
   CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)

   kywd_lngth = keyword(1:LEN(kywd_lngth))
ELSE
   kywd_lngth = keyword
END IF

!------------------------------------------------------------------------------
! Write unit if it is given.
!------------------------------------------------------------------------------
IF (PRESENT(unit)) THEN
   ! Check unit length for convention and proper formatting
   IF(LEN_TRIM(unit) .GT. LEN(unit_lngth)) THEN
      mssg = "The unit "//TRIM(unit)//" is longer than the convention allows and therefore truncated!"
      CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)
      unit_lngth = unit(1:LEN(unit_lngth))
   ELSE
      unit_lngth = unit

      ! Length MUST match Charater declaration of the length of the variable. Otherwise the calculation and the formatting are corrupted.
      unit_post_fill = 10 - LEN_TRIM(unit)
   END IF
END IF

IF(PRESENT(int_0D )) THEN
   datatype = 1 
   cntr = cntr +1_ik
END IF
IF(PRESENT(real_0D )) THEN
   datatype = 2 
   cntr = cntr +1_ik
END IF
IF(PRESENT(int_1D2)) THEN
   datatype = 3 
   cntr = cntr +1_ik
END IF
IF(PRESENT(real_1D2)) THEN
   datatype = 4 
   cntr = cntr +1_ik
END IF
IF(PRESENT(int_1D3)) THEN
   datatype = 5 
   cntr = cntr +1_ik
END IF
IF(PRESENT(real_1D3)) THEN
   datatype = 6 
   cntr = cntr +1_ik
END IF
IF(PRESENT(chars)) THEN
   datatype = 7
   cntr = cntr +1_ik
END IF

!------------------------------------------------------------------------------
! Check amount of given datatypes 
!------------------------------------------------------------------------------
mssg=''
IF (cntr .LT. 1_ik) mssg = "The datatype of keyword »"//TRIM(keyword)//"« was not defined!"
IF (cntr .GT. 1_ik) mssg = "Too many datatypes for keyword »"//TRIM(keyword)//"« were defined!"
IF (mssg .NE. '') CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 1)

!------------------------------------------------------------------------------
! Read meta input
!------------------------------------------------------------------------------
IF (PRESENT(m_in) .EQV. .TRUE.) THEN

   IF(.NOT. PRESENT(m_in) ) THEN
      mssg = 'No array of lines to parse keyword »'//TRIM(keyword)//'« given. Check subroutine »read_write_meta«.' 
      CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 1)
   ELSE
      do_loop_counter = SIZE(m_in)
   END IF

   !------------------------------------------------------------------------------
   ! Parse Data out of the input array
   !------------------------------------------------------------------------------
   lineloop: DO ii=SIZE(m_in), 1_ik, -1_ik 
      line = m_in(ii)
   
      IF (line(1:1) .EQ. '*') THEN ! it's a keyword

         CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

         IF (tokens(2) .EQ. TRIM(keyword)) THEN
            kywd_found = 1

            SELECT CASE( datatype )
               CASE(1); READ(tokens(3  ), ifmt)   int_0D 
               CASE(2); READ(tokens(3  ), rfmt)  real_0D 
               CASE(3); READ(tokens(3:4), ifmt)   int_1D2
               CASE(4); READ(tokens(3:4), rfmt)  real_1D2
               CASE(5); READ(tokens(3:5), ifmt)   int_1D3
               CASE(6); READ(tokens(3:5), rfmt)  real_1D3
               CASE(7); chars = TRIM(ADJUSTL(tokens(3))) !(stdspc-LEN_TRIM(tokens(3)) : stdspc)
            END SELECT
            
            ! Exit the loop after parsing the first occurance as its the last mentioning of the keyword. (File is read beginning from the last line)
            EXIT ! NOT CYCLE!!
         END IF
      END IF
   END DO lineloop

   ! If the keyword is not in the file: kwabrt controls whether to stop the program or not
   IF (kywd_found .EQ. 0_ik)  THEN

      mssg = "The keyword »"//TRIM(ADJUSTL(keyword))//"« was not found in the meta file!"
      CALL handle_err(fh, TRIM(ADJUSTL(mssg)), kwabrt_u)
      stat = 1
   END IF
END IF


!------------------------------------------------------------------------------
! Write meta output
!------------------------------------------------------------------------------
IF (PRESENT(m_in) .EQV. .FALSE.) THEN

   maxchars = stdspc

   WRITE(fh, '(A)', ADVANCE='NO') kywd_lngth

   ! Build format specifier and write the output
   SELECT CASE( datatype )
      CASE(1); WRITE(fh, "(30(' '),   I15   )", ADVANCE='NO')  int_0D
      CASE(2); WRITE(fh, "(30(' '),   F15.7 )", ADVANCE='NO') real_0D
      CASE(3); WRITE(fh, "(15(' '), 2(I15)  )", ADVANCE='NO')  int_1D2
      CASE(4); WRITE(fh, "(15(' '), 2(F15.7))", ADVANCE='NO') real_1D2
      CASE(5); WRITE(fh, "(         3(I15)  )", ADVANCE='NO')  int_1D3
      CASE(6); WRITE(fh, "(         3(F15.7))", ADVANCE='NO') real_1D3
      CASE(7)
         ! The width is adjusted to 3x 20 chars (3 dimensions with 15 places each).
         ! above 60 chars, the width is essentially overflowing.
         IF (LEN_TRIM(ADJUSTL(chars)) .GT. stdspc) maxchars = LEN_TRIM(ADJUSTL(chars))+1 
         WRITE(fh, "(2A)", ADVANCE='NO') REPEAT(' ', maxchars-LEN_TRIM(ADJUSTL(chars))), TRIM(ADJUSTL(chars))
   END SELECT


   IF ((PRESENT(unit)) .AND. (datatype /= 7)) THEN
      WRITE(fh, '(A)', ADVANCE='NO') ' '
      WRITE(fh, '(A)', ADVANCE='NO') unit_lngth
   ELSE
      IF (((maxchars-stdspc .GT. 0) .AND. (maxchars-stdspc .LE. ucl)) .OR. (maxchars == stdspc)) THEN
         WRITE(fh, '(A)', ADVANCE='NO') REPEAT(' ', LEN(unit_lngth)+1-(maxchars-stdspc))
      END IF
   END IF
     
   IF ((ndu .EQV. .FALSE.) .AND. (maxchars-stdspc-ucl .LT. 0)) THEN
      CALL date_time(fh, da=.TRUE., ti=.TRUE., zo=.TRUE.)
   ELSE
      WRITE(fh, '(A)') '' ! In this case, a linebreak is required
   END IF
END IF

END SUBROUTINE meta_io


!------------------------------------------------------------------------------
! SUBROUTINE: meta_close
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file and to alter its name.
!> Assign out = in before calling this routine. Also define a new app_name :-)
!
!------------------------------------------------------------------------------
SUBROUTINE meta_close()

!------------------------------------------------------------------------------
! Check and close files - Routine: (fh, filename, abrt, stat)
!------------------------------------------------------------------------------
CALL check_and_close(fhmei, TRIM(in%full) , .FALSE.)
CALL check_and_close(fhmeo, TRIM(out%full), .FALSE.)

END SUBROUTINE meta_close

END MODULE meta