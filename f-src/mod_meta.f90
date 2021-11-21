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
MODULE meta

   USE strings
   USE error_handling

IMPLICIT NONE

   !------------------------------------------------------------------------------
   ! Provide versioning information for transparent data tracking
   !------------------------------------------------------------------------------  
   INCLUDE 'include_f90/revision_meta.f90'

   INTEGER, PARAMETER :: meta_ik = 8
   INTEGER, PARAMETER :: meta_rk = 8

   ! Character lengths
   INTEGER, PARAMETER :: kcl    = 25   ! Keyword character  length
   INTEGER, PARAMETER :: ucl    = 10   ! Unit    character  length
   INTEGER, PARAMETER :: stdspc = 45   ! Keyword standard space

   CHARACTER(LEN=kcl), PARAMETER :: meta_program_keyword
   CHARACTER(LEN=kcl), PARAMETER :: meta_prgrm_mstr_app

   ! Standard files
   INTEGER(KIND=ik), PARAMETER :: fh_meta_in  = 20, fhmei  = 20
   INTEGER(KIND=ik), PARAMETER :: fh_meta_put = 21, fhmeo  = 21
   INTEGER(KIND=ik), PARAMETER :: fh_mon      = 25, fhmon  = 25
   INTEGER(KIND=ik), PARAMETER :: fh_out      = 30, fho    = 30
   INTEGER(KIND=ik), PARAMETER :: fh_log      = 35, fhl    = 35
   INTEGER(KIND=ik), PARAMETER :: fh_res      = 40, fhr    = 40
   INTEGER(KIND=ik), PARAMETER :: fh_csv      = 45, fhc    = 45
   INTEGER(KIND=ik), PARAMETER :: fh_head     = 50, fhh    = 50
   CHARACTER(LEN=*), PARAMETER :: log_suf  = '.log'
   CHARACTER(LEN=*), PARAMETER :: lock_suf = '.lock'
   CHARACTER(LEN=*), PARAMETER :: head_suf = '.head'
   CHARACTER(LEN=*), PARAMETER :: meta_suf = '.meta'
   CHARACTER(LEN=*), PARAMETER :: mon_suf  = '.mon'
   CHARACTER(LEN=*), PARAMETER :: res_suf  = '.result'
   CHARACTER(LEN=*), PARAMETER :: csv_suf  = '.csv'

   ! Meta data basename handling
   TYPE basename
      ! For the use in filenames, a max. length of a part of a basename of kcl characters must suffice.
      ! Nomenclature: dataset_type_purpose_app_features
      CHARACTER(LEN=mcl) :: full     = '' ! Including suffix and path
      CHARACTER(LEN=mcl) :: path     = '' ! Only the path to the file
      CHARACTER(LEN=mcl) :: p_n_bsnm = '' ! Just the path and the basename
      CHARACTER(LEN=mcl) :: bsnm     = '' ! Just the basename
      CHARACTER(LEN=kcl) :: dataset  = '' ! For example FH01-1 (Femoral Head 1, Scan1)
      CHARACTER(LEN=2)   :: type     = '' ! 'cl' - clinical or 'mu' - microfocus
      CHARACTER(LEN=3)   :: purpose  = '' ! 'Dev' or 'Pro' (Development or Production)
      CHARACTER(LEN=kcl) :: app      = '' ! Application. For example "Binarization"
      CHARACTER(LEN=kcl) :: features = '' ! Features. For example the parametrization
   END TYPE basename

   ! Always provide in/out for meta driven environments
   TYPE(basename) :: in, out


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
      Module Procedure meta_write_I0D 
      Module Procedure meta_write_R0D 
      Module Procedure meta_write_I1D
      Module Procedure meta_write_R1D

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

LOGICAL :: exist=.FALSE.
INTEGER  (KIND=ik) :: ios
CHARACTER(LEN=mcl) :: lockname

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
lockname=TRIM(in%path)//'.'//TRIM(in%bsnm)//lock_suf

INQUIRE (FILE = TRIM(lockname), EXIST = exist)

IF((restart .EQ. 'N') .AND. (exist)) THEN
   mssg='The .*.lock file is set and a restart prohibited by default or the user.'
   CALL handle_err(std_out, TRIM(ADJUSTL(mssg)), err=1_ik)
END IF

IF(((restart .EQ. 'Y') .AND. (.NOT. exist)) .OR. ((restart .EQ. 'N') .AND. (.NOT. exist))) THEN
   CALL execute_command_line ('touch '//TRIM(lockname), CMDSTAT=ios)
   CALL handle_err(std_out, 'The .*.lock file could not be set.', err=ios)
END IF

IF((restart .EQ. 'Y') .AND. (exist)) CONTINUE

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

LOGICAL :: exist

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
INQUIRE (FILE = TRIM(in%full), EXIST = exist)
IF (.NOT. exist) CALL handle_err(std_out, "The file "//TRIM(in%full)//" does not exist.", 1)

CALL parse( str=in%full, delims=".", args=tokens, nargs=ntokens)

IF ( '.'//TRIM(tokens(ntokens)) .EQ. meta_suf) THEN
   !------------------------------------------------------------------------------
   ! Parse all basename and path details.
   !------------------------------------------------------------------------------
   in%p_n_bsnm = in%full(1:LEN_TRIM(in%full)-LEN_TRIM(meta_suf)) 
   
   CALL parse( str=TRIM(in%p_n_bsnm), delims="/", args=tokens, nargs=ntokens)

   in%path = in%p_n_bsnm(1:LEN_TRIM(in%p_n_bsnm) - LEN_TRIM(tokens(ntokens)))     
   in%bsnm = TRIM(tokens(ntokens))

   CALL parse( str=TRIM(in%bsnm), delims="_", args=tokens, nargs=ntokens)

   in%dataset = TRIM(tokens(1))
   in%type    = TRIM(tokens(2))
   in%purpose = TRIM(tokens(3))
   in%app      = TRIM(tokens(4))
   in%features = TRIM(tokens(5))

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
          '_'//TRIM(meta_prgrm_mstr_app)//&
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
! SUBROUTINE: meta_start_ascii
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
!> If Restart isn's set or .FALSE., then the program aborts in case a meta or
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
!> @param[in] restart Logfiles (temporary and permanent)
!---------------------------------------------------------------------------  
SUBROUTINE meta_start_ascii(fh, suf, restart)

INTEGER  (KIND=ik), INTENT(IN) :: fh
CHARACTER(LEN=*), INTENT(IN) :: suf
CHARACTER, INTENT(IN), OPTIONAL :: restart

CHARACTER(LEN=mcl) :: temp_f_suf, perm_f_suf
INTEGER  (KIND=ik) :: ios
CHARACTER :: restart_u

LOGICAL :: exist_temp, exist_perm

restart_u='N'

! The temporaray file is a hidden one.
temp_f_suf = TRIM(out%path)//'.temporary'//TRIM(suf)
perm_f_suf = TRIM(out%p_n_bsnm)//TRIM(suf)
  
!------------------------------------------------------------------------------
! Check for a temporary file
! Check for a permanent file
!------------------------------------------------------------------------------
INQUIRE (FILE = temp_f_suf, EXIST = exist_temp)
INQUIRE (FILE = out%p_n_bsnm//TRIM(suf), EXIST = exist_perm)

!------------------------------------------------------------------------------
! Check restart
!------------------------------------------------------------------------------
IF(PRESENT(restart)) restart_u=restart

IF (restart_u .EQ. 'Y') THEN
   ! if target_val if inquire(exist) = .FALSE. and stat_*l = 0 - the file does not exist
   IF(exist_temp) THEN
      CALL execute_command_line ('rm -r '//TRIM(temp_f_suf), CMDSTAT=ios)   
      CALL handle_err(std_out, '»'//TRIM(temp_f_suf)//'« not deletable.',ios)
   END IF

   IF(exist_perm) THEN
      CALL execute_command_line ('rm -r '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)
      CALL handle_err(std_out, '»'//TRIM(out%full)//'« not deletable.', ios)
   END IF

ELSE ! restart_u .EQ. 'N'
   IF ((exist_temp) .OR. (exist_perm)) THEN

      IF ((exist_temp) .AND. (exist_perm)) THEN 
         mssg='The file '//TRIM(perm_f_suf)//' and the file '//TRIM(temp_f_suf)//' already exists.'
      ELSE IF  (exist_temp) THEN
         mssg='The file '//TRIM(temp_f_suf)//' already exists.'
      ELSE ! (exist_perm) 
         mssg='The file '//TRIM(perm_f_suf)//' already exists.'
      END IF

      CALL handle_err(std_out, TRIM(mssg)//' Previous job maybe was aborted.', 1)     
   END IF
END IF

OPEN(UNIT=fh, FILE=TRIM(temp_f_suf), ACTION='WRITE', ACCESS='SEQUENTIAL', STATUS='NEW')

END SUBROUTINE meta_start_ascii

!------------------------------------------------------------------------------
! SUBROUTINE: meta_stop_ascii
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
SUBROUTINE meta_stop_ascii(fh, suf)

INTEGER  (KIND=ik), INTENT(IN) :: fh
CHARACTER(LEN=*), INTENT(IN) :: suf

CHARACTER(LEN=mcl) :: temp_f_suf, perm_f_suf
INTEGER  (KIND=ik) :: ios

temp_f_suf = TRIM(out%path)//'.temporary'//TRIM(suf)
perm_f_suf = TRIM(out%p_n_bsnm)//TRIM(suf)

CLOSE (fh)

!------------------------------------------------------------------------------
! The temporary log file must be renamed to a permanent one
!------------------------------------------------------------------------------
CALL execute_command_line ('mv '//TRIM(temp_f_suf)//' '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)

IF(ios /= 0_ik) THEN
   mssg='Can not rename the suffix_file from »'//TRIM(temp_f_suf)//'« to the proper basename.'
   CALL handle_err(std_out, mssg, 0)
END IF

END SUBROUTINE meta_stop_ascii

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

kywd_lngth = ''

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
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_C (fh, keyword, m_in, chars)
   
INTEGER  (KIND=ik), INTENT(IN)  :: fh 
CHARACTER(LEN=*)  , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
CHARACTER(LEN=*), INTENT(OUT) :: chars 

! Internal variables
CHARACTER(LEN=mcl) :: tokens(30)
INTEGER(KIND=ik) :: ntokens, ii
INTEGER(KIND=ik) :: kywd_found
LOGICAL :: override

kywd_found = 0
override = .FALSE.

CALL check_keyword(fh, keyword)

!------------------------------------------------------------------------------
! Parse Data out of the input array
!------------------------------------------------------------------------------
DO ii = , SIZE(m_in) 
   CALL parse(str=m_in(ii), delims=' ', args=tokens, nargs=ntokens)

! ACCEPT IF NOT WRITTEN IN PROGRAMS OWN SCOPE

   SELECT CASE(tokens(1))
      CASE('*', 'r', 'w')
         IF (tokens(2) == TRIM(keyword)) THEN
            kywd_found = 1
            chars = TRIM(ADJUSTL(tokens(3))) !(stdspc-LEN_TRIM(tokens(3)) : stdspc)

            !------------------------------------------------------------------------------
            ! Keep keyword if it appears in the programs scope
            !------------------------------------------------------------------------------
            IF ((override) .AND. (tokens(1) /= 'w')) EXIT

         END IF
      CASE('p')
         !------------------------------------------------------------------------------
         ! User tells, that the scope of another program begins.
         !------------------------------------------------------------------------------
         IF ((.NOT. override) .AND. (tokens(2) == TRIM(meta_program_keyword))) THEN
            override = .TRUE.
         ELSE
            !------------------------------------------------------------------------------
            ! Leave, if the scope of the next program begins.
            !------------------------------------------------------------------------------
            EXIT
         END IF

   END SELECT
END DO

IF (kywd_found == 0)  CALL keyword_error(fh, keyword)

END SUBROUTINE meta_read_C

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
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
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
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_R0D (fh, keyword, m_in, real_0D)
     
INTEGER(KIND=ik), INTENT(IN)  :: fh 
CHARACTER(LEN=*), INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
REAL(KIND=rk), INTENT(OUT) :: real_0D 

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

         !------------------------------------------------------------------------------
         ! Exit the loop after parsing the first occurance as its the last 
         ! mentioning of the keyword. (File is read beginning from the last line)
         !------------------------------------------------------------------------------
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
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_I1D (fh, keyword, m_in, int_1D)

INTEGER(KIND=ik) , INTENT(IN)  :: fh 
CHARACTER(LEN=*) , INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
INTEGER(KIND=ik), DIMENSION(:), INTENT(OUT) :: int_1D 

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

         !------------------------------------------------------------------------------
         ! Exit the loop after parsing the first occurance as its the last 
         ! mentioning of the keyword. (File is read beginning from the last line)
         !------------------------------------------------------------------------------
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
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_read_R1D (fh, keyword, m_in, real_1D)

INTEGER(KIND=ik), INTENT(IN)  :: fh 
CHARACTER(LEN=*), INTENT(IN)  :: keyword
CHARACTER(LEN=mcl), DIMENSION(:), INTENT(IN)  :: m_in      
REAL(KIND=rk), DIMENSION(:), INTENT(OUT) :: real_1D 

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

         !------------------------------------------------------------------------------
         ! Exit the loop after parsing the first occurance as its the last 
         ! mentioning of the keyword. (File is read beginning from the last line)
         !------------------------------------------------------------------------------
         EXIT
      END IF
   END IF
END DO

IF (kywd_found .EQ. 0_ik)  CALL keyword_error(fh, keyword)
END SUBROUTINE meta_read_R1D



!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_keyword
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write finalized strings of keywords
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill String with data
!> @param[in] unit Unit of the value
!---------------------------------------------------------------------------
SUBROUTINE meta_write_keyword (fh, keyword, stdspcfill, unit)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: stdspcfill 
CHARACTER(LEN=*), INTENT(IN) :: unit

CHARACTER(LEN=scl) :: fmt, str
CHARACTER(LEN=8)  :: date
CHARACTER(LEN=10) :: time
CHARACTER(LEN=5)  :: timezone

CALL check_keyword(fh, keyword)
CALL check_unit(fh, unit)

WRITE(fmt, '(A,I0,A)') "(2A, T", kcl, ")"
WRITE(fh, fmt, ADVANCE='NO') "w ", keyword

WRITE(fmt, '(A,I0,A)') "(A, T", stdspc+1, ")"
WRITE(fh, fmt, ADVANCE='NO') TRIM(ADJUSTL(stdspcfill))

WRITE(fmt, '(A,I0,A)') "(A, T", ucl+1, ")"
WRITE(fh, fmt, ADVANCE='NO') unit
   
CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

str = ''
str = date(7:8)//'.'//date(5:6)//'.'//date(1:4)
str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
str = TRIM(str)//' '//timezone


WRITE(fh, '(A)') TRIM(str)

END SUBROUTINE meta_write_keyword



!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_sha256sum
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write finalized strings of keywords
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill String with data
!> @param[in] unit Unit of the value
!---------------------------------------------------------------------------
SUBROUTINE meta_write_sha256sum (binary_name)
   
CHARACTER(LEN=*), INTENT(IN) :: binary_name

CHARACTER(LEN=kcl-1) :: keyword = ''
CHARACTER(LEN=scl) :: fmt, stdspcfill
INTEGER(KIND=ik), DIMENSION(5) :: stat = 0
INTEGER(KIND=ik) :: ios

LOGICAL :: exist

!---------------------------------------------------------------------------
! Write "Keyword"
!---------------------------------------------------------------------------
keyword = "* SHA256SUM_OF_BINARY"

WRITE(fmt, '(A,I0,A)') "(2A, T", kcl, ")"
WRITE(fhmeo, fmt, ADVANCE='NO') keyword

!---------------------------------------------------------------------------
! Check the buffer file
!---------------------------------------------------------------------------
INQUIRE(FILE = 'temp_buffer', EXIST = exist)

IF (exist) THEN
   CALL EXECUTE_COMMAND_LINE ('rm -r temp_buffer', CMDSTAT=ios)   

   IF(ios /= 0_ik) THEN
      mssg='Can not delete the temp_buffer'
      CALL handle_err(std_out, mssg, 0)
      stat(1) = 1      
   END IF
END IF

!---------------------------------------------------------------------------
! Check for auxiliary programs
!---------------------------------------------------------------------------
CALL EXECUTE_COMMAND_LINE("which cut > /dev/null 2> /dev/null", CMDSTAT=stat(2))
CALL EXECUTE_COMMAND_LINE("which sha256sum > /dev/null 2> /dev/null", CMDSTAT=stat(3))

!---------------------------------------------------------------------------
! Deal with the buffer file
!---------------------------------------------------------------------------
IF(SUM(stat)==0) THEN
   OPEN(UNIT=9, FILE='temp_buffer', ACTION='READWRITE', STATUS='NEW')

   CALL EXECUTE_COMMAND_LINE("sha256sum "//TRIM(ADJUSTL(binary_name))//" | cut -d ' ' -f 1 >> 'temp_buffer'", CMDSTAT=stat(4))

   READ(9, '(A)', iostat=stat(5)) stdspcfill

   CLOSE(9)
END IF

IF (SUM(stat) == 0) THEN
   WRITE(fhmeo, fmt) TRIM(ADJUSTL(stdspcfill))
ELSE
   WRITE(fhmeo, fmt) "Could not get sha256sum. One of the previious system calls failed."
END IF

INQUIRE(FILE='temp_buffer', EXIST=exist)

IF (exist) CALL EXECUTE_COMMAND_LINE ('rm -r temp_buffer', CMDSTAT=ios)   

END SUBROUTINE meta_write_sha256sum


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_C
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords with character output. 
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill Characters to write
!---------------------------------------------------------------------------
SUBROUTINE meta_write_C (fh, keyword, stdspcfill)
   
INTEGER  (KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*)  , INTENT(IN) :: keyword
CHARACTER(LEN=*)  , INTENT(IN) :: stdspcfill 

CALL meta_write_keyword (fh, keyword, stdspcfill, '')

END SUBROUTINE meta_write_C

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_I0D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type integer dim 0.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] int_0D Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_write_I0D (fh, keyword, unit, int_0D)
   
INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
INTEGER(KIND=ik), INTENT(IN) :: int_0D 

CHARACTER(LEN=scl) :: stdspcfill

WRITE(stdspcfill, '(I0)') int_0D

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_I0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_R0D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type Real dim 0.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] real_0D Datatype to read in
!---------------------------------------------------------------------------
SUBROUTINE meta_write_R0D (fh, keyword, unit, real_0D)
   
INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
REAL(KIND=ik), INTENT(IN) :: real_0D 

CHARACTER(LEN=scl) :: stdspcfill

WRITE(stdspcfill, '(F0.0)') real_0D

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_R0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_I1D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type integer dim 1.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] int_0D Datatype
!---------------------------------------------------------------------------
SUBROUTINE meta_write_I1D (fh, keyword, unit, int_1D)
   
INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
INTEGER(KIND=ik), INTENT(IN), DIMENSION(:) :: int_1D 

CHARACTER(LEN=scl) :: stdspcfill, str
INTEGER  (KIND=ik) :: ii

stdspcfill = ''
str = ''

DO ii=1, SIZE(int_1D)
   str = ''
   WRITE(str, '(I0)') int_1D(ii)
   stdspcfill = TRIM(stdspcfill)//' '//TRIM(str)
END DO

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_I1D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_R1D
!---------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type Real dim 1.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] real_1D Datatype
!---------------------------------------------------------------------------
SUBROUTINE meta_write_R1D (fh, keyword, unit, real_1D)
   
INTEGER(KIND=ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
REAL(KIND=ik), INTENT(IN), DIMENSION(:) :: real_1D 

CHARACTER(LEN=scl) :: stdspcfill, str
INTEGER  (KIND=ik) :: ii

stdspcfill = ''
str = ''

DO ii=1, SIZE(real_1D)
   str = ''
   WRITE(str, '(F0.0)') real_1D(ii)
   stdspcfill = TRIM(stdspcfill)//' '//TRIM(str)
END DO

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_R1D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_signing
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file.
!
!> @description
!> Requires a "revision.meta" or similar inclusion of verisoning info, 
!> provided by a makefile. Furhermore, it requires a global_stds file.
!------------------------------------------------------------------------------
SUBROUTINE meta_signing(binary_name)

CHARACTER(LEN=*), INTENT(IN)  :: binary_name

WRITE(fhmeo, '(A)')
CALL meta_write (fhmeo, 'PROGRAM_VERSION' , revision)
CALL meta_write (fhmeo, 'PROGRAM_GIT_HASH' , hash)

CALL meta_write_sha256sum (binary_name)

CALL meta_write (fhmeo, 'COMPUTATION_FINISHED' , 'Succesfully')

END SUBROUTINE meta_signing


!------------------------------------------------------------------------------
! SUBROUTINE: meta_close
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file.
!
!> @description
!> provided by a makefile. Furhermore, it requires a global_stds file.
!------------------------------------------------------------------------------
SUBROUTINE meta_close()

LOGICAL :: opened

WRITE(fhmeo, '(A)')
WRITE(fhmeo, SEP_STD)

!------------------------------------------------------------------------------
! Check and close files - Routine: (fh, filename, abrt, stat)
!------------------------------------------------------------------------------
INQUIRE(UNIT=fhmei, OPENED=opened)
IF(opened) CLOSE (fhmei)

INQUIRE(UNIT=fhmeo, OPENED=opened)
IF(opened) CLOSE (fhmeo)

 
END SUBROUTINE meta_close

END MODULE meta