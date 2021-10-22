Module metis

  use, intrinsic :: iso_c_binding

  implicit none

  INTERFACE

     Subroutine F_METIS_PartMeshDual( ne, nn, eptr, eind, &
          vwgt, vsize, ncommon, nparts, tpwgts, options, objval, &
          epart, npart) BIND(C, Name="F_METIS_PartMeshDual") 

       use, intrinsic :: iso_c_binding

       !> ne : Number of elements
       !> nn : number of all nodes describing the topology
       !>      (E.g: for two connected 8 node bricks this number is 16)
       Integer  (kind=C_INT64_T)                      :: ne, nn
       Integer  (kind=C_INT64_T)  , Dimension(ne+1)   :: eptr
       Integer  (kind=C_INT64_T)  , Dimension(nn)     :: eind
       Integer  (kind=C_INT64_T)  , Dimension(ne)     :: vwgt, vsize
       Integer  (kind=C_INT64_T)                      :: ncommon, nparts
       Real     (kind=c_double)  , Dimension(nparts)  :: tpwgts
       Integer  (kind=C_INT64_T) , Dimension(40)      :: options
       Integer  (kind=C_INT64_T)                      :: objval
       Integer  (kind=C_INT64_T)  , Dimension(ne)     :: epart
       Integer  (kind=C_INT64_T)  , Dimension(nn)     :: npart

     End Subroutine F_METIS_PartMeshDual
 
  END INTERFACE
  
  INTERFACE

     Subroutine F_METIS_PartMeshNodal( ne, nn, eptr, eind, &
          vwgt, vsize, nparts, tpwgts, options, objval, &
          epart, npart) BIND(C, Name="F_METIS_PartMeshNodal") 

       use, intrinsic :: iso_c_binding

       !> ne : Number of elements
       !> nn : number of all odes describing the topology
       !>      (E.g: for two connected 8 node bricks this number is 16)
       Integer  (kind=c_long)                       :: ne, nn
       Integer  (kind=c_long)   , Dimension(ne+1)   :: eptr
       Integer  (kind=c_long)   , Dimension(nn)     :: eind
       Integer  (kind=c_long)   , Dimension(ne)     :: vwgt, vsize
       Integer  (kind=c_long)                       :: nparts
       Real     (kind=c_double), Dimension(nparts)  :: tpwgts
       Integer  (kind=c_long)  , Dimension(40)      :: options
       Integer  (kind=c_long)                       :: objval
       Integer  (kind=c_long)   , Dimension(ne)     :: epart
       Integer  (kind=c_long)   , Dimension(nn)     :: npart

     End Subroutine F_METIS_PartMeshNodal

  END INTERFACE

End Module metis
