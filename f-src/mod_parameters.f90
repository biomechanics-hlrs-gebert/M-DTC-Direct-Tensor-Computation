!------------------------------------------------------------------------------
! MODULE: petsc_opt
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Set PETSc options.
!------------------------------------------------------------------------------
module petsc_opt

  USE global_std
  use petsc
  
  implicit none

contains
  
  Subroutine Set_PETSc_Options()

    Integer(kind=mik) :: petsc_ierr
    
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-mat_type"    , "mpiaij",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_type"     , "jacobi",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_type"    , "cg"    ,petsc_ierr)
    ! Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-start_in_debugger"    , ""  ,petsc_ierr)
    
    IF(out_amount == "DEBUG") THEN
      Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_monitor" , ""      ,petsc_ierr)
    END IF

    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_rtol"    , "1.e-06",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_max_it"  , "500000",petsc_ierr)

    ! CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-malloc_test"  , "true",petsc_ierr)
    ! CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-log_view"  , "log.view",petsc_ierr)
    ! CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-objects_dump" , "all",petsc_ierr)
    ! CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-malloc_debug" , "true",petsc_ierr)
    ! CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-malloc_dump"  , "true",petsc_ierr)

  End Subroutine Set_PETSc_Options
  
end module petsc_opt
