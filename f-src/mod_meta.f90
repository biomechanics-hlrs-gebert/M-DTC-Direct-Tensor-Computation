




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

  !> Interface: meta_read
  !> \author Johannes Gebert
  !> \date 10.11.2021
  Interface meta_read

     Module Procedure meta_read_C 
     Module Procedure meta_read_I0D 
     Module Procedure meta_read_I1D
     Module Procedure meta_read_R0D
     Module Procedure meta_read_R1D

  End Interface meta_read

  !> Interface: meta_write
  !> \author Johannes Gebert
  !> \date 10.11.2021
  Interface meta_write

     Module Procedure meta_write_C 
   !   Module Procedure meta_write_I0D 
   !   Module Procedure meta_write_R0D 
   !   Module Procedure meta_write_I1D2
   !   Module Procedure meta_write_R1D2
   !   Module Procedure meta_write_I1D3 
   !   Module Procedure meta_write_R1D3

  End Interface meta_write

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

CHARACTER, INTENT(IN) :: restart

LOGICAL               :: exist=.FALSE.
INTEGER  (KIND=ik)    :: ios
CHARACTER(LEN=mcl)    :: lockname

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
CHARACTER(LEN=mcl) :: line
INTEGER  (KIND=ik) :: ios, lines, ii
CHARACTER(LEN=mcl) :: tokens(30)
INTEGER  (KIND=ik) :: ntokens

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
CALL meta_read (fhmon, 'NEW_BSNM_FEATURE', meta_as_rry, out%features)
IF((ios ==1)) out%features = in%features           ! if ios = 1 --> no keyword found

CALL meta_read (fhmon, 'NEW_BSNM_PURPOSE', meta_as_rry, out%purpose)
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
OPEN(UNIT=fhmeo, FILE=TRIM(out%full), ACTION='WRITE', ACCESS='APPEND', STATUS='OLD')

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

INTEGER  (KIND=ik), INTENT(IN)            :: fh
CHARACTER(LEN=*)  , INTENT(IN)            :: suf
CHARACTER(LEN=*)  , INTENT(IN)            :: st
CHARACTER         , INTENT(IN) , OPTIONAL :: restart

CHARACTER(LEN=mcl) :: temp_f_suf, perm_f_suf
INTEGER  (KIND=ik) :: ios, stat_t, stat_p
CHARACTER          :: restart_u='N'


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
! SUBROUTINE: check_keyword
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Routine to truncate a keyword which was too long. Could do other stuff 
!> as well.
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword to check
!------------------------------------------------------------------------------  
SUBROUTINE check_keyword(fh, keyword)

INTEGER  (KIND=ik) :: fh 
CHARACTER(LEN=*)   :: keyword
CHARACTER(LEN=kcl) :: kywd_lngth

IF(LEN_TRIM(keyword) .GT. LEN(kywd_lngth)) THEN
   WRITE(fh, '(A)') ''
   mssg = "The keyword »"//TRIM(keyword)//"« is longer than the &
   &convention allows and therefore truncated!"
   CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)

   kywd_lngth = keyword(1:LEN(kywd_lngth))
ELSE
   kywd_lngth = keyword
END IF

END SUBROUTINE check_keyword

!------------------------------------------------------------------------------
! SUBROUTINE: check_unit
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Routine to truncate a keyword which was too long. Could do other stuff 
!> as well.
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword to check
!------------------------------------------------------------------------------  
SUBROUTINE check_unit(fh, unit)

INTEGER  (KIND=ik) :: fh 
CHARACTER(LEN=*)   :: unit
CHARACTER(LEN=ucl) :: unit_lngth

! Check unit length for convention and proper formatting
IF(LEN_TRIM(unit) .GT. LEN(unit_lngth)) THEN
   mssg = "The unit "//TRIM(unit)//" is longer than the convention allows and therefore truncated!"
   CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)
   unit_lngth = unit(1:LEN(unit_lngth))
ELSE
   unit_lngth = unit
END IF

END SUBROUTINE check_unit


!------------------------------------------------------------------------------
! SUBROUTINE: keyword_error
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper of handle_err
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword 
!------------------------------------------------------------------------------  
SUBROUTINE keyword_error(fh, keyword)

INTEGER  (KIND=ik) :: fh 
CHARACTER(LEN=*)   :: keyword

   mssg = "The keyword »"//TRIM(ADJUSTL(keyword))//"« was not found in the meta file!"
   CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 1)

END SUBROUTINE keyword_error


!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_I0D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
! 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars    Optional datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_I0D (fh, keyword, m_in, int_0D)
     
INTEGER  (KIND=ik)              , INTENT(IN)  :: fh 
CHARACTER(LEN=*)                , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
INTEGER  (KIND=ik)              , INTENT(OUT) :: int_0D 

! Internal variables
CHARACTER(LEN=mcl)   :: tokens(30), line
INTEGER  (KIND=ik)   :: ntokens, ii
INTEGER  (KIND=ik)   :: kywd_found

!------------------------------------------------------------------------------
! Initialize variables
! Should be defined from new each time the program enters the routine. 
! Otherwise some errors may occur.
!------------------------------------------------------------------------------
kywd_found     = 0

CALL check_keyword(fh, keyword)

!------------------------------------------------------------------------------
! Parse Data out of the input array
!------------------------------------------------------------------------------
DO ii=SIZE(m_in), 1_ik, -1_ik 
   line = m_in(ii)

   IF (line(1:1) .EQ. '*') THEN ! it's a keyword

      CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

      IF (tokens(2) .EQ. TRIM(keyword)) THEN
         kywd_found = 1

         READ(tokens(3), ifmt) int_0D 
         
         !------------------------------------------------------------------------------
         ! Exit the loop after parsing the first occurance as its the last 
         ! mentioning of the keyword. (File is read beginning from the last line)
         !------------------------------------------------------------------------------
         EXIT
      END IF
   END IF
END DO

IF (kywd_found .EQ. 0_ik)  CALL keyword_error(fh, keyword)
END SUBROUTINE meta_read_I0D


!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_R0D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
! 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] real_0D    Optional datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_R0D (fh, keyword, m_in, real_0D)
     
INTEGER  (KIND=ik)              , INTENT(IN)  :: fh 
CHARACTER(LEN=*)                , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
REAL     (KIND=rk)              , INTENT(OUT) :: real_0D 

! Internal variables
CHARACTER(LEN=mcl)   :: tokens(30), line
INTEGER  (KIND=ik)   :: ntokens, ii
INTEGER  (KIND=ik)   :: kywd_found

kywd_found = 0

CALL check_keyword(fh, keyword)

DO ii=SIZE(m_in), 1_ik, -1_ik 
   line = m_in(ii)

   IF (line(1:1) .EQ. '*') THEN ! it's a keyword

      CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

      IF (tokens(2) .EQ. TRIM(keyword)) THEN
         kywd_found = 1
         READ(tokens(3), rfmt) real_0D 
         EXIT
      END IF
   END IF
END DO

IF (kywd_found .EQ. 0_ik)  CALL keyword_error(fh, keyword)
END SUBROUTINE meta_read_R0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_I1D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
! 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] int_1D    Optional datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_I1D (fh, keyword, m_in, int_1D)

INTEGER  (KIND=ik)              , INTENT(IN)  :: fh 
CHARACTER(LEN=*)                , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
INTEGER  (KIND=ik), DIMENSION(:), INTENT(OUT) :: int_1D 

! Internal variables
CHARACTER(LEN=mcl)   :: tokens(30), line
INTEGER  (KIND=ik)   :: ntokens, ii
INTEGER  (KIND=ik)   :: kywd_found

kywd_found = 0

CALL check_keyword(fh, keyword)

DO ii=SIZE(m_in), 1_ik, -1_ik 
   line = m_in(ii)

   IF (line(1:1) .EQ. '*') THEN ! it's a keyword

      CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

      IF (tokens(2) .EQ. TRIM(keyword)) THEN
         kywd_found = 1
         READ(tokens(3:2+SIZE(int_1D)), ifmt) int_1D
         EXIT
      END IF
   END IF
END DO

IF (kywd_found .EQ. 0_ik)  CALL keyword_error(fh, keyword)
END SUBROUTINE meta_read_I1D


!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_R1D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
! 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] real_1D    Optional datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_R1D (fh, keyword, m_in, real_1D)

INTEGER  (KIND=ik)              , INTENT(IN)  :: fh 
CHARACTER(LEN=*)                , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
REAL     (KIND=rk), DIMENSION(:), INTENT(OUT) :: real_1D 

! Internal variables
CHARACTER(LEN=mcl)   :: tokens(30), line
INTEGER  (KIND=ik)   :: ntokens, ii
INTEGER  (KIND=ik)   :: kywd_found

kywd_found = 0

CALL check_keyword(fh, keyword)

DO ii=SIZE(m_in), 1_ik, -1_ik 
   line = m_in(ii)

   IF (line(1:1) .EQ. '*') THEN ! it's a keyword

      CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

      IF (tokens(2) .EQ. TRIM(keyword)) THEN
         kywd_found = 1
         READ(tokens(3:2+SIZE(real_1D)), rfmt) real_1D
         EXIT
      END IF
   END IF
END DO

IF (kywd_found .EQ. 0_ik)  CALL keyword_error(fh, keyword)
END SUBROUTINE meta_read_R1D



!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_C
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to parse keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
! 
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars    Optional datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_C (fh, keyword, m_in, chars)
   
INTEGER  (KIND=ik)              , INTENT(IN)  :: fh 
CHARACTER(LEN=*)                , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
CHARACTER(LEN=*)                , INTENT(OUT) :: chars 

! Internal variables
CHARACTER(LEN=mcl)   :: tokens(30), line
INTEGER  (KIND=ik)   :: ntokens, ii
INTEGER  (KIND=ik)   :: kywd_found

kywd_found     = 1

CALL check_keyword(fh, keyword)

!------------------------------------------------------------------------------
! Parse Data out of the input array
!------------------------------------------------------------------------------
DO ii = SIZE(m_in), 1_ik, -1_ik 
   line = m_in(ii)

   IF (line(1:1) .EQ. '*') THEN ! it's a keyword

      CALL parse(str=line, delims=' ', args=tokens, nargs=ntokens)

      IF (tokens(2) .EQ. TRIM(keyword)) THEN
         kywd_found = 0

         chars = TRIM(ADJUSTL(tokens(3))) !(stdspc-LEN_TRIM(tokens(3)) : stdspc)
         
         !------------------------------------------------------------------------------
         ! Exit the loop after parsing the first occurance as its the last 
         ! mentioning of the keyword. (File is read beginning from the last line)
         !------------------------------------------------------------------------------
         EXIT 
      END IF
   END IF
END DO

!------------------------------------------------------------------------------
! If the keyword is not in the file:
!------------------------------------------------------------------------------
IF (kywd_found .EQ. 1_ik)  THEN
   mssg = "The keyword »"//TRIM(ADJUSTL(keyword))//"« was not found in the meta file!"
   CALL handle_err(fh, TRIM(ADJUSTL(mssg)), kywd_found)
END IF
END SUBROUTINE meta_read_C


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_C
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords with character output. 
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] chars    Datatype to write
!> @param[in] kwabrt Default -> .TRUE.; do abort
!---------------------------------------------------------------------------

SUBROUTINE meta_write_C (fh, keyword, chars)
   
INTEGER  (KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*)  , INTENT(IN) :: keyword
CHARACTER(LEN=*)  , INTENT(IN) :: chars 

! Internal variables
CHARACTER(LEN=ucl)   :: unit_lngth
INTEGER  (KIND=ik)   :: maxchars

maxchars = stdspc

CALL check_keyword(fh, keyword)

WRITE(fh, '(2A)', ADVANCE='NO') "* ", keyword

! The width is adjusted to 3x 20 chars (3 dimensions with 15 places each).
! above 60 chars, the width is essentially overflowing.
IF (LEN_TRIM(ADJUSTL(chars)) .GT. stdspc) maxchars = LEN_TRIM(ADJUSTL(chars))+1 
WRITE(fh, "(2A)", ADVANCE='NO') REPEAT(' ', maxchars-LEN_TRIM(ADJUSTL(chars))), TRIM(ADJUSTL(chars))


IF (((maxchars-stdspc .GT. 0) .AND. (maxchars-stdspc .LE. ucl)) .OR. (maxchars == stdspc)) THEN
   WRITE(fh, '(A)', ADVANCE='NO') REPEAT(' ', LEN(unit_lngth)+1-(maxchars-stdspc))
END IF

CALL date_time(fh, da=.TRUE., ti=.TRUE., zo=.TRUE.)

END SUBROUTINE meta_write_C

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_I0D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type integer kind 1. 
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Data to read from the meta files
!> @param[in] unit Unit of the value in case a text output is requested.
!> @param[in] int__0D  Optional datatype to read in
!> @param[in] kwabrt Default -> .TRUE.; do abort
!---------------------------------------------------------------------------

! SUBROUTINE meta_write_I0D (fh, keyword, unit, int_0D, kwabrt)
   
! INTEGER  (KIND=ik), INTENT(IN)              :: fh 
! CHARACTER(LEN=*)  , INTENT(IN)              :: keyword
! CHARACTER(LEN=*)  , INTENT(IN)   , OPTIONAL :: unit
! INTEGER  (KIND=ik), INTENT(IN)              :: int_0D 
! INTEGER  (KIND=ik), INTENT(IN)   , OPTIONAL :: kwabrt

! ! Internal variables
! CHARACTER(LEN=kcl)   :: kywd_lngth
! CHARACTER(LEN=ucl)   :: unit_lngth
! INTEGER  (KIND=ik)   :: maxchars
! INTEGER  (KIND=ik)   :: unit_post_fill, kwabrt_u

! IF(PRESENT(kwabrt)) kwabrt_u = kwabrt

! maxchars = stdspc

! !------------------------------------------------------------------------------
! ! Check unit keyword for convention and proper formatting
! !------------------------------------------------------------------------------
! IF(LEN_TRIM(keyword) .GT. LEN(kywd_lngth)) THEN
!    mssg = "»"//TRIM(keyword)//"« is to long and therefore truncated before reading!"
!    CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)

!    kywd_lngth = keyword(1:LEN(kywd_lngth))
! ELSE
!    kywd_lngth = keyword
! END IF

! !------------------------------------------------------------------------------
! ! Write unit if it is given.
! !------------------------------------------------------------------------------
! IF (PRESENT(unit)) THEN
!    ! Check unit length for convention and proper formatting
!    IF(LEN_TRIM(unit) .GT. LEN(unit_lngth)) THEN
!       mssg = "The unit "//TRIM(unit)//" is to long and therefore truncated!"
!       CALL handle_err(fh, TRIM(ADJUSTL(mssg)), 0)
!       unit_lngth = unit(1:LEN(unit_lngth))
!    ELSE
!       unit_lngth = unit

!       ! Length MUST match Charater declaration of the length of the variable. Otherwise the calculation and the formatting are corrupted.
!       unit_post_fill = 10 - LEN_TRIM(unit)
!    END IF
! END IF

! WRITE(fh, '(2A)', ADVANCE='NO') "* ", kywd_lngth

! WRITE(fh, "(30(' '),I15)", ADVANCE='NO')  int_0D

! IF ((PRESENT(unit)) THEN
!    WRITE(fh, '(A)', ADVANCE='NO') ' '
!    WRITE(fh, '(A)', ADVANCE='NO') unit_lngth
! ELSE
!    IF (((maxchars-stdspc .GT. 0) .AND. (maxchars-stdspc .LE. ucl)) .OR. (maxchars == stdspc)) THEN
!       WRITE(fh, '(A)', ADVANCE='NO') REPEAT(' ', LEN(unit_lngth)+1-(maxchars-stdspc))
!    END IF
! END IF
   
! CALL date_time(fh, da=.TRUE., ti=.TRUE., zo=.TRUE.)

! END SUBROUTINE meta_write_I0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_close
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file and to alter its name.
!> Assign out = in before calling this routine. Also define a new app_name :-)
!
!> @param[in] m_in Array of lines of ascii meta file
!------------------------------------------------------------------------------
SUBROUTINE meta_close()

CALL meta_write (fhmeo, 'COMPUTATION_FINISHED' , '')

WRITE(fhmeo, '(A)')
WRITE(fhmeo, SEP_STD)

!------------------------------------------------------------------------------
! Check and close files - Routine: (fh, filename, abrt, stat)
!------------------------------------------------------------------------------
CALL check_and_close(fhmei, TRIM(in%full) , .FALSE.)
CALL check_and_close(fhmeo, TRIM(out%full), .FALSE.)

END SUBROUTINE meta_close

END MODULE meta