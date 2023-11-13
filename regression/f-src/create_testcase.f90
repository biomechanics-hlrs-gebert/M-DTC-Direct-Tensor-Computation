!---------------------------------------------------------------------------------------------------
!> \mainpage HLRS - Test cases for SIR (Spatial Image Registration)
!>
!> <hr>
!> \section desc Description
!>
!> This program creates testcases to check whether specific parts work properly.
!>
!> <hr>
!> \section developers Developers
!>
!> Johannes Gebert
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on: 17.09.2023
!---------------------------------------------------------------------------------------------------
PROGRAM create_DTC_testcases

USE global_std
USE user_interaction
USE formatted_plain
USE ser_binary

IMPLICIT NONE

INTEGER(ik), PARAMETER :: std_array_ik = 4, shik = 4
INTEGER(ik), DIMENSION(3) :: dims
INTEGER(ik) :: baseline_HU, upper_hu
INTEGER(INT32), DIMENSION(:,:,:), ALLOCATABLE  :: array

REAL(rk), DIMENSION(3) :: spcng

CHARACTER(mcl) :: fl_out, mode="", purpose=""

LOGICAL :: exist = .FALSE.


!------------------------------------------------------------------------------
! Filename/Metadata:
!------------------------------------------------------------------------------
CALL GET_COMMAND_ARGUMENT(1, mode)
CALL GET_COMMAND_ARGUMENT(2, purpose)

SELECT CASE (TRIM(ADJUSTL(mode)))
CASE ("isotropic")
   CONTINUE
CASE ("orthotropic")
   CONTINUE
CASE DEFAULT
   WRITE(std_out,FMT_ERR) "No valid mode given: cmd arg 1 = <mode>"
   WRITE(std_out,FMT_ERR) "    'isotroic'"
   WRITE(std_out,FMT_ERR) "    'orthotropic'"
   WRITE(std_out,FMT_ERR) "    'monotropic'"
   GOTO 1000
END SELECT

IF (TRIM(ADJUSTL(fl_out))=="") THEN
   WRITE(std_out,FMT_ERR) "No valid basename purpose given:"
   WRITE(std_out,FMT_ERR) "Cmd arg 2 = ./TC01-1_test_dev_dtc_<mode>-<purpose>.raw"
   GOTO 1000
ELSE
   IF (TRIM(ADJUSTL(purpose))=="") THEN
      fl_out = "TC01-1_test_dev_dtc_"//TRIM(ADJUSTL(mode))//".raw"
   ELSE
      fl_out = "TC01-1_test_dev_dtc_"//TRIM(ADJUSTL(mode))//"-"//TRIM(ADJUSTL(purpose))//".raw"
   END IF 
END IF 


CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "])

WRITE(std_out,FMT_TXT) "Modify source code to create and to alter testcases"
WRITE(std_out,FMT_TXT) "Testcase selected:   "//TRIM(ADJUSTL(mode))
WRITE(std_out,FMT_TXT) "Output file defined: "//TRIM(ADJUSTL(fl_out))
WRITE(std_out,FMT_SEP) 

!------------------------------------------------------------------------------
! Standard test structure and allocation
!------------------------------------------------------------------------------
baseline_HU =     0_shik
upper_HU    = 14250_shik

dims  = [ 100_ik , 100_ik , 100_ik ]
spcng = [ 0.01_rk, 0.01_rk, 0.01_rk ]

ALLOCATE(array(dims(1),dims(2),dims(3)))
array = baseline_HU

!------------------------------------------------------------------------------
! Actual test cases
!------------------------------------------------------------------------------
SELECT CASE (TRIM(ADJUSTL(mode)))
CASE ("isotropic")
   array = upper_HU

   ! Write parameters to std_out to inform the user 
   ! about the parameters used for the test case.
   WRITE(std_out, FMT_TXT_AxI0) "Spcng:", dims
   WRITE(std_out, FMT_TXT_AxF0) "Dims: ", spcng
CASE ("orthotropic")
   array(56:65,   :  ,   :  ) = upper_HU
   array(  :  , 56:65,   :  ) = upper_HU

   WRITE(std_out, FMT_TXT_AxI0) "Spcng:", dims
   WRITE(std_out, FMT_TXT_AxF0) "Dims: ", spcng

CASE DEFAULT
   GOTO 1000
END SELECT

!------------------------------------------------------------------------------
! Write to file 
!------------------------------------------------------------------------------
INQUIRE(FILE = TRIM(fl_out), EXIST=exist)

IF(exist) THEN
   WRITE(std_out,FMT_SEP) 
   WRITE(std_out,FMT_ERR) "Output file already exists. Please mv or rm."
   GOTO 1000
ELSE

   CALL ser_write_binary(94_ik, fl_out, array)

END IF 

1000 CONTINUE

IF(ALLOCATED(array)) DEALLOCATE(array)

WRITE(std_out,FMT_SEP) 

END PROGRAM create_DTC_testcases
