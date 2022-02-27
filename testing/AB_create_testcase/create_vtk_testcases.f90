!-------------------------------------------------------------------------------------------------
! HLRS - NUM - vtk STURCUTRED_POINTS data manipulation tool
!
! This program creates testcases to check whether specific parts work properly.
! It's more like a script which needs to be modified according to your specific needs.
!
! Johannes Gebert
! Created:  14.05.2021
! last_mod: 09.02.2022
!-------------------------------------------------------------------------------------------------
PROGRAM create_scalar_vtk_cases

USE ISO_FORTRAN_ENV

IMPLICIT NONE

REAL(KIND=REAL64), DIMENSION (3) :: spcng

INTEGER(KIND=INT32), PARAMETER     :: fl_in_un=22, fl_out_un=23
INTEGER(KIND=INT32), PARAMETER     :: ik=4, rk=8
INTEGER(KIND=INT32), DIMENSION (3) :: dims, flh, cross_dims, low_bnds, upp_bnds
INTEGER(KIND=INT32), DIMENSION (:,:,:), ALLOCATABLE :: array
INTEGER(KIND=INT32) :: ii, xx, yy, zz, xxx, yyy, zzz, intervall
INTEGER(KIND=INT32) :: steps, step_width_sf, step_width_sv
INTEGER(KIND=INT32) :: remainder_sv, remainder_sf, baseline_sv, max_sv, mov_val

CHARACTER(LEN=512) :: fl_in, fl_out, struc

! -----------------------------------------------------------------------------
! Manual parametrization
! -----------------------------------------------------------------------------
fl_out = "/zhome/academic/HLRS/hlrs/hpcgeber/hpcgeber-dtc_hawk/datasets/FH01-2_mu_Dev_ctbi_isotropic_20.raw"
! fl_out = "/home/geb/00_bone_eval_chain/M-DTC-Direct Tensor Computation/&
!    &testing/AB_create_testcase/TC00-0_mu_test_2dcross_ortho.raw"
! TC00-0_mu_test_2dcross_ortho.vtk
! -----------------------------------------------------------------------------
! Structural details
! -----------------------------------------------------------------------------
dims  = [ 2904, 2910, 2104 ]                        ! vtk header
spcng = [ 0.0049800, 0.0049800, 0.0049800 ]         ! vtk header

! -----------------------------------------------------------------------------
! Manual parametrization
! -----------------------------------------------------------------------------
baseline_sv = 0_ik
max_sv = 1_ik
steps  = 1_ik

! Target structure
! struc = "UNIFORM"     ! Uniform scalar field
! struc = "STAIRCASE_X" ! Value increasing from x=0 to x=dim(1)
                        ! Make sure, steps fit into dimensions without fractions :-)
! struc = "CHEQUERED"   ! Voxel by Voxel - lo - hi - lo - hi - lo - ...
struc = "PERCENTAGE_PLANES"
   intervall=5 ! yields in a density of 0.2
! struc = "ZEBRA"
! struc = "3DCROSS"         
! struc = "2DCROSS" 
! struc = "PLANE"   
! struc = "BAR"     
! struc = "InnerCross" 
! struc = "RubiksCube"  
! struc = "RubiksCube2"
! struc = "Test"


! Programs task
ALLOCATE(array(dims(1),dims(2),dims(3)))

! UNIFORM
IF(struc == "UNIFORM") THEN
   array(:,:,:) = baseline_sv
END IF

! STAIRCASE_X
IF(struc == "STAIRCASE_X") THEN
   remainder_sv  = MODULO(max_sv-baseline_sv,steps)
   step_width_sv = ((max_sv-baseline_sv)-remainder_sv)/steps

   remainder_sf  = MODULO(dims(1),steps)
   step_width_sf = (dims(1)-remainder_sf)/steps
   DO ii=1, steps
      array(1_ik+ii*step_width_sf-step_width_sf:ii*step_width_sf,:,:) = (ii-1_ik)*step_width_sv+baseline_sv
   END DO

   IF(remainder_sf/=0_ik) THEN
      array(dims(1)-remainder_sf:dims(1),:,:) = max_sv
   END IF
END IF

WRITE(*,'(A,I6)') "baseline_sv: ",baseline_sv
WRITE(*,'(A,I6)') "     max_sv: ",max_sv

! PERCENTAGE_PLANES
IF(struc == "PERCENTAGE_PLANES") THEN
   array (:,:,:) = baseline_sv

   DO XX=1, dims(3), intervall
      array (xx,:,:) = max_sv
   END DO
END IF

! CHEQUERED
IF(struc == "CHEQUERED") THEN
   array (:,:,:) = max_sv
   DO xxx=1,dims(1)-1_ik
      IF (MODULO(xxx,2_ik) == 1_ik) THEN
         xx=xxx
      ELSE
         xx=xxx+1_ik
      END IF
      DO yyy=1,dims(2)-1_ik
         IF (MODULO(yyy,2_ik) == 1_ik) THEN
            yy=yyy
         ELSE
            yy=yyy+1_ik
         END IF
         DO zzz=1,dims(3)-1_ik,2
            array(xx,yy,zzz) = baseline_sv
         END DO
      END DO
   END DO
END IF

! ZEBRA
IF(struc == "ZEBRA") THEN
   array (:,:,:) = max_sv
   DO xx=1,dims(1),2
      array(xx,:,:)=baseline_sv
   END DO
END IF

! 2DCROSS
IF(struc == '2DCROSS') THEN
   array (:,:,:) = baseline_sv

   ! cross_dims = dims/4_ik ! Integer division - simply removes decimal places
   
   ! low_bnds = FLOOR(cross_dims * 1.1_rk)
   ! upp_bnds = FLOOR(cross_dims * 2.1_rk)

   low_bnds = 8  
   upp_bnds = 11

   array (low_bnds(1):upp_bnds(1), :, :) = INT(max_sv, KIND=4)
   array (:, low_bnds(2):upp_bnds(2), :) = INT(max_sv, KIND=4)
   ! array (9:10, 8:12, :) = INT(max_sv, KIND=4)
   ! array (:, 9:10, :) = INT(max_sv, KIND=4)
END IF

! 3DCROSS
IF(struc == '3DCROSS') THEN
   array (:,:,:) = baseline_sv

   cross_dims = dims/4_ik ! Integer division - simply removes decimal places
   
   low_bnds = FLOOR(cross_dims * 1.5_rk)
   upp_bnds = FLOOR(cross_dims * 2.5_rk)

   array (low_bnds(1):upp_bnds(1), :, :) = INT(max_sv, KIND=2_ik)
   array (:, low_bnds(2):upp_bnds(2), :) = INT(max_sv, KIND=2_ik)
   array (:, :, low_bnds(3):upp_bnds(3)) = INT(max_sv, KIND=2_ik)
END IF

! PLANE
IF(struc == 'PLANE') THEN
   array (:,:,:) = baseline_sv

   cross_dims = dims/4_ik ! Integer division - simply removes decimal places
   
   low_bnds = FLOOR(cross_dims * 1.5_rk)
   upp_bnds = FLOOR(cross_dims * 2.5_rk)

   array (low_bnds(1):upp_bnds(1), :, :) = INT(max_sv, KIND=2_ik)
END IF

! BAR
IF(struc == 'BAR') THEN
   array (:,:,:) = baseline_sv

   cross_dims = dims/4_ik ! Integer division - simply removes decimal places
   
   low_bnds = FLOOR(cross_dims * 1.5_rk)
   upp_bnds = FLOOR(cross_dims * 2.5_rk)
   array (low_bnds(1):upp_bnds(1), low_bnds(2):upp_bnds(2), :) = INT(max_sv, KIND=2_ik)
END IF

! InnerCross
IF(struc == 'InnerCross') THEN
   array (:,:,:) = baseline_sv

   cross_dims = dims/4_ik ! Integer division - simply removes decimal places
   
   low_bnds = FLOOR(cross_dims * 1.5_rk)
   upp_bnds = FLOOR(cross_dims * 2.5_rk)
   array (low_bnds(1):upp_bnds(1), low_bnds(2):upp_bnds(2), :) = INT(max_sv, KIND=2_ik)
   array (:, low_bnds(2):upp_bnds(2), low_bnds(1):upp_bnds(1)) = INT(max_sv, KIND=2_ik)
   array (low_bnds(1):upp_bnds(3), :, low_bnds(1):upp_bnds(1)) = INT(max_sv, KIND=2_ik)

END IF   

! RubiksCube
IF(struc == 'RubiksCube') THEN
   array (:,:,:) = max_sv

   DO xx=2,dims(1),50
   DO yy=2,dims(2),50
   DO zz=3,dims(3),50
   
   
   
      array(xx, : , :)=baseline_sv
      array(: , yy, :)=baseline_sv
      array(: , : ,zz)=baseline_sv

   END DO   
   END DO
   END DO


END IF

! RubiksCube2
IF(struc == 'RubiksCube2') THEN
   array (:,:,:) = max_sv

   DO xx=2,dims(1),5
   DO yy=2,dims(2),25
   DO zz=3,dims(3),50
      
   
      array(xx, : , :)=baseline_sv
      array(: , yy, :)=baseline_sv
      array(: , : ,zz)=baseline_sv

   END DO   
   END DO
   END DO


END IF

! Test
IF(struc == "Test") THEN

   array (:,:,:) = max_sv

   
      DO xxx=1,dims(1)-1_ik
      IF (MODULO(xxx,2_ik) == 1_ik) THEN
         xx=xxx
      ELSE
         xx=xxx+1_ik
      END IF
   
   !DO xx =1*2_ik,dims(1),10
   DO yy=2,dims(2),25
   DO zz=3,dims(3),50
   
   
   
      array(xx, : , :)=baseline_sv
      array(: , yy, :)=baseline_sv
      array(: , : ,zz)=baseline_sv

   END DO   
   END DO
   END DO


END IF

! Write new raw File
OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), &
      ACCESS="stream", FORM="unformatted", STATUS="new")
WRITE(UNIT=fl_out_un) array(:,:,:)
CLOSE(UNIT=fl_out_un)


! Write new vtk File
! OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="new")
! WRITE(fl_out_un,'(A)')            "# vtk DataFile Version 5.1"
! WRITE(fl_out_un,'(A)')            "vtk output"
! WRITE(fl_out_un,'(A)')            "BINARY"
! WRITE(fl_out_un,'(A)')            "DATASET STRUCTURED_POINTS"
! WRITE(fl_out_un,'(A,3(I5,A))')    "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
! WRITE(fl_out_un,'(A,3(F11.6,A))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3)
! WRITE(fl_out_un,'(A)')            "ORIGIN 0 0 0"
! WRITE(fl_out_un,'(A, I11)')       "POINT_DATA", dims(1)*dims(2)*dims(3)
! WRITE(fl_out_un,'(A)')            "SCALARS DICOMImage int"            ! tricky part of this stuff!
! WRITE(fl_out_un,'(A)')            "LOOKUP_TABLE default"
! CLOSE(UNIT=fl_out_un)

! OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), CONVERT='big_endian', &
!       ACCESS="stream", FORM="unformatted", STATUS="old", POSITION="append")
! WRITE(UNIT=fl_out_un) array(:,:,:)
! CLOSE(UNIT=fl_out_un)

! OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="old", POSITION="append")
! WRITE(fl_out_un ,'(A)')          "METADATA"
! WRITE(fl_out_un ,'(A)')          "INFORMATION 0"
! WRITE(fl_out_un ,'(A)')
! CLOSE(UNIT=fl_out_un)

END PROGRAM create_scalar_vtk_cases