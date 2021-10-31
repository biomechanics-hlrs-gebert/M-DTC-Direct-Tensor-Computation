!------------------------------------------------------------------------------
! MODULE: standards
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and numbers.
!------------------------------------------------------------------------------
MODULE global_std

USE ISO_FORTRAN_ENV

IMPLICIT NONE

! General constants
INTEGER            , PARAMETER :: sik           = 2    ! INTEGER Kind
INTEGER            , PARAMETER :: ik            = 8    ! INTEGER Kind
INTEGER            , PARAMETER :: rk            = 8    ! Real Kind
INTEGER            , PARAMETER :: mcl           = 512  ! Maximal character  length
INTEGER            , PARAMETER :: hcl           = 256  ! Half    character  length
INTEGER            , PARAMETER :: scl           = 64   ! Short   character  length
INTEGER            , PARAMETER :: kcl           = 30   ! Keyword character  length
INTEGER            , PARAMETER :: kdescl        = 45   ! Keyword descriptor length

! Puredat constants
INTEGER            , PARAMETER :: pd_ik         = 8    ! Puredat Integer kind parameter
INTEGER            , PARAMETER :: pd_rk         = 8    ! Puredat Real    kind parameter
INTEGER            , PARAMETER :: pd_mpi_ik     = 4    ! Puredat Integer MPI kind parameter

!-- File handles, debug_lvl and suffix
INTEGER(KIND=ik)   , PARAMETER :: dbg_lvl       = 1
CHARACTER(LEN=mcl)             :: mssg          = ''

INTEGER(KIND=ik)   , PARAMETER :: std_in        = 5
INTEGER(KIND=ik)   , PARAMETER :: std_out       = 6
INTEGER(KIND=ik)   , PARAMETER :: std_err       = 0

! »Standard« data types
INTEGER(KIND=ik)   , PARAMETER :: fh_meta       = 20, fhme  = 20
INTEGER(KIND=ik)   , PARAMETER :: fh_mon        = 25, fhmo  = 25
INTEGER(KIND=ik)   , PARAMETER :: fh_out        = 30, fho   = 30
INTEGER(KIND=ik)   , PARAMETER :: fh_log        = 35, fhl   = 35
INTEGER(KIND=ik)   , PARAMETER :: fh_res        = 40, fhr   = 40
INTEGER(KIND=ik)   , PARAMETER :: fh_csv        = 45, fhc   = 45
CHARACTER(LEN=4)   , PARAMETER :: log_suf       = '.log'
CHARACTER(LEN=5)   , PARAMETER :: lock_suf      = '.lock'
CHARACTER(LEN=5)   , PARAMETER :: meta_suf      = '.meta'
CHARACTER(LEN=4)   , PARAMETER :: mon_suf       = '.mon'
CHARACTER(LEN=7)   , PARAMETER :: res_suf       = '.result'
CHARACTER(LEN=4)   , PARAMETER :: csv_suf       = '.csv'


CHARACTER(len=5)               :: creturn       = achar(13)
CHARACTER(len=5)               :: ifmt          = '(I10)'       ! general integer    format
CHARACTER(len=8)               :: rfmt          = '(F30.10)'    ! general real       format
CHARACTER(len=10)              :: sfmt          = '(E23.15E2)'  ! general scientific format

!-- Mpi-specific kinds
INTEGER            , PARAMETER :: mpi_ik        = 4             ! MPI INTEGER Kind; Compile with corresponding mpi!!

! Meta data basename handling
TYPE basename
   ! For the use in filenames, a max. length of a part of a basename of kcl characters must suffice.
   ! Nomenclature: dataset_type_purpose_app_features
   CHARACTER(LEN=mcl)         :: full          = ''            ! Including suffix and path
   CHARACTER(LEN=mcl)         :: path          = ''            ! Only the path to the file
   CHARACTER(LEN=mcl)         :: p_n_bsnm      = ''            ! Just the path and the basename
   CHARACTER(LEN=mcl)         :: bsnm          = ''            ! Just the basename
   !
   CHARACTER(LEN=kcl)         :: dataset       = ''            ! For example FH01-1 (Femoral Head 1, Scan1)
   CHARACTER(LEN=2)           :: type          = ''            ! 'cl' - clinical or 'mu' - microfocus
   CHARACTER(LEN=3)           :: purpose       = ''            ! 'Dev' or 'Pro' (Development or Production)
   CHARACTER(LEN=kcl)         :: app           = ''            ! Application. For example "Binarization"
   CHARACTER(LEN=kcl)         :: features      = ''            ! Features. For example the parametrization
END TYPE basename

! Always provide in/out for meta driven environments
TYPE(basename)                                     :: in, out

END MODULE global_std


!------------------------------------------------------------------------------
! MODULE: mechanical_standards
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and vmbers.
!
! REVISION HISTORY:
! 26 09 2021 - Initial Version
!------------------------------------------------------------------------------
MODULE mechanical_standards

USE global_std

IMPLICIT NONE
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

END MODULE mechanical_standards


!==============================================================================
!> Global constants and parameters for the puredat data handling library
!> \author Ralf Schneider
!> \date 22.01.2010
!>
Module global_pd

Implicit none

!> Number of currently used stream variables in puredat_streams
!>
!> The total number of currently used stream variables which is the  
!> number of arrays defined in puredat_streams independently from 
!> their data type
Integer, Parameter :: no_streams = 7

!> Maximum character length used in puredat library
Integer, Parameter :: pd_mcl = 512
!> Maximum Character Length in pd_ik elements
Integer, Parameter :: pd_ce  = 512/8

! Character constants for nice output ---------------------------------------
Character(Len=*), Parameter :: PDF_E_A    = "('EE ',A)"
Character(Len=*), Parameter :: PDF_E_AI0  = "('EE ',*(A,1X,I0))"
Character(Len=*), Parameter :: PDF_E_STOP = &
      "('EE PROGRAM STOPPED ..... ',/,'<',78('='),'>')"

Character(Len=*), Parameter :: PDF_W_A    = "('WW ',A)"
Character(Len=*), Parameter :: PDF_W_AI0  = "('WW ',*(A,1X,I0))"

Character(Len=*), Parameter :: PDF_M_A    = "('MM ',A)"
Character(Len=*), Parameter :: PDF_M_AI0  = "('MM ',A,1X,I0)"

Character(Len=*), Parameter :: PDF_TIME   = "('MM ',A,1X,F0.6,' sec')"

Character(Len=*), Parameter :: PDF_SEP    = "('<',78('='),'>')"

!> puredat project path
!>
!> Path to puredat project files which means header-, stream- and 
!> log-files
Character(len=pd_mcl) :: pro_path
!> puredat project name
!>
!> Base name of the puredat project files which are ...
Character(len=pd_mcl) :: pro_name

!   !> puredat monitor unit
Integer               :: pd_umon  != OUTPUT_UNIT
  
End Module global_pd
