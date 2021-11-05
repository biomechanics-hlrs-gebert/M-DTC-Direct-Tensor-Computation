!------------------------------------------------------------------------------
! MODULE: standards
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and numbers.
!------------------------------------------------------------------------------
MODULE global_std

IMPLICIT NONE

! General constants
INTEGER            , PARAMETER :: sik           = 2    ! INTEGER Kind
INTEGER            , PARAMETER :: ik            = 8    ! INTEGER Kind
INTEGER            , PARAMETER :: rk            = 8    ! Real    Kind
INTEGER            , PARAMETER :: mcl           = 512  ! Maximal character  length
INTEGER            , PARAMETER :: hcl           = 256  ! Half    character  length
INTEGER            , PARAMETER :: scl           = 64   ! Short   character  length
INTEGER            , PARAMETER :: kcl           = 25   ! Keyword character  length
INTEGER            , PARAMETER :: ucl           = 10   ! Unit    character  length
INTEGER            , PARAMETER :: stdspc        = 45   ! Keyword standard space

!-- File handles, debug_lvl and suffix
INTEGER(KIND=ik)   , PARAMETER :: timer_level   = 3 ! 1 ! 2
INTEGER(KIND=ik)   , PARAMETER :: dbg_lvl       = 1
CHARACTER(LEN=mcl)             :: mssg          = ''

INTEGER(KIND=ik)   , PARAMETER :: std_in        = 5
INTEGER(KIND=ik)   , PARAMETER :: std_out       = 6
INTEGER(KIND=ik)   , PARAMETER :: std_err       = 0

! Standard files
INTEGER(KIND=ik)   , PARAMETER :: fh_meta_in    = 20, fhmei  = 20
INTEGER(KIND=ik)   , PARAMETER :: fh_meta_put   = 21, fhmeo  = 21
INTEGER(KIND=ik)   , PARAMETER :: fh_mon        = 25, fhmon  = 25
INTEGER(KIND=ik)   , PARAMETER :: fh_out        = 30, fho    = 30
INTEGER(KIND=ik)   , PARAMETER :: fh_log        = 35, fhl    = 35
INTEGER(KIND=ik)   , PARAMETER :: fh_res        = 40, fhr    = 40
INTEGER(KIND=ik)   , PARAMETER :: fh_csv        = 45, fhc    = 45
INTEGER(KIND=ik)   , PARAMETER :: fh_head       = 50, fhh    = 50
CHARACTER(LEN=*)   , PARAMETER :: log_suf       = '.log'
CHARACTER(LEN=*)   , PARAMETER :: lock_suf      = '.lock'
CHARACTER(LEN=*)   , PARAMETER :: head_suf      = '.head'
CHARACTER(LEN=*)   , PARAMETER :: meta_suf      = '.meta'
CHARACTER(LEN=*)   , PARAMETER :: mon_suf       = '.mon'
CHARACTER(LEN=*)   , PARAMETER :: res_suf       = '.result'
CHARACTER(LEN=*)   , PARAMETER :: csv_suf       = '.csv'
 
!------------------------------------------------------------------------------
! Standard formats
!------------------------------------------------------------------------------
CHARACTER(len=5)   , PARAMETER :: creturn       = achar(13)
CHARACTER(len=5)   , PARAMETER :: ifmt          = '(I10)'       ! general integer    format
CHARACTER(len=8)   , PARAMETER :: rfmt          = '(F30.10)'    ! general real       format
CHARACTER(len=10)  , PARAMETER :: sfmt          = '(E23.15E2)'  ! general scientific format

! Character constants for nice output ---------------------------------------
Character(Len=*), Parameter :: fmt_sep    = "('<',97('='),'>')"
Character(LEN=*), Parameter :: fmt_inpsep = "('+',99('-'))"

Character(Len=*), Parameter :: FMT_MSG     = "('MM ',A,T97,' MM')"
Character(Len=*), Parameter :: FMT_MSG_BS  = "('MM ',A,T88,' ... ',$)"
Character(Len=*), Parameter :: FMT_MSG_BE  = "('done MM')"

Character(Len=*), Parameter :: FMT_WRN     = "('WW ',A,T97,' WW')"  
Character(Len=*), Parameter :: FMT_ERR     = "('EE ',A,T97,' EE')"
CHARACTER(LEN=*), PARAMETER :: FMT_WRN_SO  = "('\x1B[35m','WW ','\x1B[0m',A,T97)" ! std_out
CHARACTER(LEN=*), PARAMETER :: FMT_ERR_SO  = "('\x1B[31m','EE ','\x1B[0m',A,T97)" ! std_out
CHARACTER(LEN=*), PARAMETER :: FMT_ERR_SOC = "('\x1B[31m','EE ',A,T97)" ! std_out - complete color

Character(Len=*), Parameter :: FMT_ERR_AI0 = "('EE ',*(A,I0))"  
CHARACTER(Len=*), PARAMETER :: FMT_ERR_A   = "('EE ',A)"

Character(Len=*), Parameter :: FMT_STOP    = "('EE PROGRAM STOPPED ..... ',&
                                             &T97,' EE',/,'<',97('='),'>')"

Character(Len=*), Parameter :: FMT_TIME = "('MM ',A,1X,F0.6,' sec')"

! Time/ date formats
CHARACTER(Len=*), PARAMETER :: FMT_DA       = "(A,'.',A,'.',A    )"
CHARACTER(Len=*), PARAMETER :: FMT_TI       = "(x,A,':',A,':',A,x)"
CHARACTER(Len=*), PARAMETER :: FMT_ZO       = "(A)"

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

! Seperators
CHARACTER(Len=*), PARAMETER :: FMT_HY_SEP  = "(112('-'))"
CHARACTER(Len=*), PARAMETER :: FMT_EQ_SEP  = "(112('='))"
CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP = "('#DBG#',95('='))"

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
CHARACTER(LEN=*), PARAMETER ::  FMT_Blck  = "\x1B[30m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Red   = "\x1B[31m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Grn   = "\x1B[32m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Orng  = "\x1B[33m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Blue  = "\x1B[34m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Prpl  = "\x1B[35m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Cyan  = "\x1B[36m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Gray  = "\x1B[37m"
CHARACTER(LEN=*), PARAMETER ::  FMT_noc   = "\x1B[0m"

!-- Mpi-specific kinds
INTEGER         , PARAMETER :: mpi_ik     = 4             ! MPI INTEGER Kind; Compile with corresponding mpi!!

! Meta data basename handling
TYPE basename
   ! For the use in filenames, a max. length of a part of a basename of kcl characters must suffice.
   ! Nomenclature: dataset_type_purpose_app_features
   CHARACTER(LEN=mcl)         :: full          = ''            ! Including suffix and path
   CHARACTER(LEN=mcl)         :: path          = ''            ! Only the path to the file
   CHARACTER(LEN=mcl)         :: p_n_bsnm      = ''            ! Just the path and the basename
   CHARACTER(LEN=mcl)         :: bsnm          = ''            ! Just the basename
   CHARACTER(LEN=kcl)         :: dataset       = ''            ! For example FH01-1 (Femoral Head 1, Scan1)
   CHARACTER(LEN=2)           :: type          = ''            ! 'cl' - clinical or 'mu' - microfocus
   CHARACTER(LEN=3)           :: purpose       = ''            ! 'Dev' or 'Pro' (Development or Production)
   CHARACTER(LEN=kcl)         :: app           = ''            ! Application. For example "Binarization"
   CHARACTER(LEN=kcl)         :: features      = ''            ! Features. For example the parametrization
END TYPE basename

! Always provide in/out for meta driven environments
TYPE(basename)                :: in, out

END MODULE global_std


!------------------------------------------------------------------------------
! MODULE: mechanical
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and vmbers.
!
! REVISION HISTORY:
! 26 09 2021 - Initial Version
!------------------------------------------------------------------------------
MODULE mechanical

USE global_std

IMPLICIT NONE

! Add other parameters if necessary.
TYPE materialcard
    REAL(KIND=rk)               :: E
    REAL(KIND=rk)               :: nu
    ! For use with effective nummerical stiffness calculations 
    REAL(KIND=rk), DIMENSION(3) :: pdsize ! Physical domain/ Macro element size 
END TYPE materialcard

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION: lamee_lambda
!------------------------------------------------------------------------------ 
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the lamé constant lambda
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return lambda
!------------------------------------------------------------------------------  
FUNCTION lamee_lambda(E, v) RESULT (lambda)
 
   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     lambda

   lambda = E*v/((1._rk+v)*(1._rk-2._rk*v))

END FUNCTION lamee_lambda

!------------------------------------------------------------------------------
! FUNCTION: lamee_mu_shear
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the lamé constant µ / the shear modulus G
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return G
!------------------------------------------------------------------------------  
FUNCTION lamee_mu_shear(E, v) RESULT (G)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     G ! (shear modulus) 

   G = E / (2._rk*(1._rk+v))

   ! Just for info
   ! mu = G 

END FUNCTION lamee_mu_shear


!------------------------------------------------------------------------------
! FUNCTION: bulk_modulus
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the bulk modulus
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return k
!------------------------------------------------------------------------------  
FUNCTION bulk_modulus(E, v) RESULT (k)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     k ! (shear modulus) 

   k = E / (3._rk*(1._rk-(2._rk*v)))

END FUNCTION bulk_modulus


!------------------------------------------------------------------------------
! FUNCTION: iso_compliance_voigt
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank compliance tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso_inv
!------------------------------------------------------------------------------  
FUNCTION iso_compliance_voigt(E, v) RESULT (t_iso_inv)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     t_iso_inv

   REAL (KIND=rk)                   ::     fctr

   fctr = 1._rk/E

!  sum of symmetric components - therefore: eps_12+eps_21 => 2*eps_12 etc.
   t_iso_inv(:,1)=(/ 1._rk,    -v,    -v,  .0_rk,           .0_rk,                  .0_rk    /)
   t_iso_inv(:,2)=(/    -v, 1._rk,    -v,  .0_rk,           .0_rk,                  .0_rk    /)
   t_iso_inv(:,3)=(/    -v,    -v, 1._rk,  .0_rk,           .0_rk,                  .0_rk    /)
   t_iso_inv(:,4)=(/ .0_rk, .0_rk, .0_rk,  2._rk*(1._rk+v), .0_rk,                  .0_rk    /)
   t_iso_inv(:,5)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,           2._rk*(1._rk+v),        .0_rk    /)
   t_iso_inv(:,6)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,           .0_rk,           2._rk*(1._rk+v) /)

   t_iso_inv = t_iso_inv*fctr

END FUNCTION iso_compliance_voigt


!------------------------------------------------------------------------------
! FUNCTION: iso_compliance_kelvin
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank compliance tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso_inv
!------------------------------------------------------------------------------  
FUNCTION iso_compliance_kelvin(E, v) RESULT (t_iso_inv)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     t_iso_inv

   REAL (KIND=rk)                   ::     fctr

   fctr = 1._rk/E

   t_iso_inv(:,1)=(/ 1._rk,    -v,    -v,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,2)=(/    -v, 1._rk,    -v,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,3)=(/    -v,    -v, 1._rk,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,4)=(/ .0_rk, .0_rk, .0_rk,  1._rk+v, .0_rk,   .0_rk   /)
   t_iso_inv(:,5)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,   1._rk+v, .0_rk   /)
   t_iso_inv(:,6)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,   .0_rk,   1._rk+v /)

   t_iso_inv = t_iso_inv*fctr

END FUNCTION iso_compliance_kelvin



!------------------------------------------------------------------------------
! FUNCTION: iso_stiffness_voigt
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank stiffness tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso
!------------------------------------------------------------------------------  
FUNCTION iso_stiffness_voigt(E, v) RESULT (t_iso)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     t_iso

   REAL (KIND=rk)                   ::     fctr, fctr_shear
   
   fctr       = E / ((1._rk+v)*(1._rk-(2._rk*v)))
   fctr_shear = 1._rk-(2._rk*v)
   t_iso(:,1)=(/ 1._rk-v,  v,     v,       .0_rk, .0_rk, .0_rk /)
   t_iso(:,2)=(/    v , 1._rk-v,  v,       .0_rk, .0_rk, .0_rk /)
   t_iso(:,3)=(/    v ,     v, 1._rk-v,    .0_rk, .0_rk, .0_rk /)
   t_iso(:,4)=(/ .0_rk, .0_rk, .0_rk, fctr_shear/2, .0_rk, .0_rk /)
   t_iso(:,5)=(/ .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear/2, .0_rk /)
   t_iso(:,6)=(/ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear/2 /)

   t_iso = t_iso*fctr ! Elementwise operation
END FUNCTION iso_stiffness_voigt


!------------------------------------------------------------------------------
! FUNCTION: iso_stiffness_kelvin
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank stiffness tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso
!------------------------------------------------------------------------------  
FUNCTION iso_stiffness_kelvin(E, v) RESULT (t_iso)

   REAL (KIND=rk)                   ::     E, v
   REAL (KIND=rk), DIMENSION(6,6)   ::     t_iso

   REAL (KIND=rk)                   ::     fctr, fctr_shear
   
   fctr       = E / ((1._rk+v)*(1._rk-(2._rk*v))) 
   fctr_shear = 1._rk-(2._rk*v)

   ! Factor of 2 not required as for normal components: sigma=E* eps
   ! Factor of 2 not required as for shear  components: sigma=E*2eps
   t_iso(:,1)=(/ 1._rk-v,  v,     v,       .0_rk, .0_rk, .0_rk /)
   t_iso(:,2)=(/    v , 1._rk-v,  v,       .0_rk, .0_rk, .0_rk /)
   t_iso(:,3)=(/    v ,     v, 1._rk-v,    .0_rk, .0_rk, .0_rk /)
   t_iso(:,4)=(/ .0_rk, .0_rk, .0_rk, fctr_shear, .0_rk, .0_rk /)
   t_iso(:,5)=(/ .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear, .0_rk /)
   t_iso(:,6)=(/ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear /)

   t_iso = t_iso*fctr ! Elementwise operation

END FUNCTION iso_stiffness_kelvin

END MODULE mechanical