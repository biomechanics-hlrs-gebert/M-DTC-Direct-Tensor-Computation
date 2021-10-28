Module Operating_System
  
  implicit none

  INTERFACE

     Subroutine stat_dir(dir, status) &
          BIND(C, Name="stat_dir")

       use, intrinsic :: iso_c_binding
       use            :: standards
       
       character(kind=c_char), dimension(4*mcl) :: dir
       Integer  (kind=c_int )                   :: status
       
     End Subroutine stat_dir

  END INTERFACE

  INTERFACE

     Subroutine make_dir(dir, status) &
          BIND(C, Name="make_dir")

       use, intrinsic :: iso_c_binding
       use            :: standards
       
       character(kind=c_char), dimension(4*mcl) :: dir
       Integer  (kind=c_int )                   :: status
       
     End Subroutine make_dir

  END INTERFACE

End Module Operating_System
