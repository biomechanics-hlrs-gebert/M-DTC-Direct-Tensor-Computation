!******************************************************************************
!**                                                                          **
!** kind module                                                              **
!**                                                                          **
!** by Ralf Schneider                                                        **
!**                                                                          **
!** ToDo: Make this correct with the ISO_FORTRAN_ENV module !                **
module kinds

  implicit none

  Integer, Parameter :: ik     = 8   ! Integer Kind
  Integer, Parameter :: rk     = 8   ! Real Kind
  Integer, Parameter :: mcl    = 512 ! maximum character length

  Integer, Parameter :: mpi_ik = 4   ! MPI Integer kind

end module kinds

