!******************************************************************************
!**                                                                          **
!** petsc_options module                                                     **
!**                                                                          **
!** by Ralf Schneider                                                        **
!**                                                                          **
module petsc_opt

  USE global_std
  use petsc
  
  implicit none

contains
  
  Subroutine Set_PETSc_Options()

    Integer(kind=mpi_ik) :: petsc_ierr
    
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-mat_type"    , "mpiaij",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_type"     , "jacobi",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_type"    , "cg"    ,petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_monitor" , ""      ,petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_rtol"    , "1.e-03",petsc_ierr)
    Call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_max_it"  , "500000",petsc_ierr)

  End Subroutine Set_PETSc_Options
  
end module petsc_opt
