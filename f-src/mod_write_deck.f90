Module write_deck

use linFE
use mesh_partitioning
use strings

implicit none

Real(rk), Parameter :: delta_b = 1.E-9_rk
Character(*) , Parameter :: fmt_filename = "(A,I0,A,I0,A)"

INTERFACE determine_prescribed_displ
    MODULE PROCEDURE determine_prescribed_displ_hexe8
    MODULE PROCEDURE determine_prescribed_displ_hexe20
END INTERFACE determine_prescribed_displ

CONTAINS

!------------------------------------------------------------------------------  
! FUNCTION: determine_prescribed_displ_hexe8
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Determine prescribed displacements
!------------------------------------------------------------------------------  
Subroutine determine_prescribed_displ_hexe8(node, xa, xe, elnodes, bnode_vals)
    
    REAL(rk), INTENT(IN), DIMENSION(3) :: node, xa, xe
    INTEGER(ik), INTENT(IN) :: elnodes
    REAL(rk), DIMENSION(:,:), INTENT(OUT) :: bnode_vals

    REAL(rk), DIMENSION(elnodes, 3 ,elnodes*3) :: uu 
    INTEGER(ik) :: ii

    ! Loadcase init
    uu = init_displ_hexe8(elnodes)

    Do ii = 1, elnodes*3
       bnode_vals(ii,1) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,1,ii))
       bnode_vals(ii,2) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,2,ii))
       bnode_vals(ii,3) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,3,ii))
    End Do

  End Subroutine determine_prescribed_displ_hexe8

!------------------------------------------------------------------------------
! FUNCTION: determine_prescribed_displ_hexe20 
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Determine prescribed displacements
!------------------------------------------------------------------------------  
Subroutine determine_prescribed_displ_hexe20(node, xa, xe, elnodes, bnode_vals)

    REAL(rk), INTENT(IN), DIMENSION(3) :: node, xa, xe
    INTEGER(ik), INTENT(IN) :: elnodes
    REAL(rk), DIMENSION(:,:), INTENT(OUT) :: bnode_vals

    REAL(rk), DIMENSION(elnodes, 3 ,elnodes*3) :: uu 
    INTEGER(ik) :: ii

    ! Loadcase init
    uu = init_displ_hexe20(elnodes)

    Do ii = 1, elnodes*3
        bnode_vals(ii,1) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,1,ii))
        bnode_vals(ii,2) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,2,ii))
        bnode_vals(ii,3) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,3,ii))
    End Do

End Subroutine determine_prescribed_displ_hexe20

!------------------------------------------------------------------------------
! FUNCTION: init_displ_hexe8 
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Initialisation of loadcases
!------------------------------------------------------------------------------  
Function init_displ_hexe8(elnodes) Result(uu)

    Integer      , Intent(in)                     :: elnodes

    Real(Kind=rk), dimension(elnodes,3,elnodes*3) :: uu 
    Real(Kind=rk)                                 :: eps=1.E-6_rk

    uu = 0._rk

    uu(:,1,1)   = [ 0._rk, eps, eps, 0._rk, 0._rk, eps, eps, 0._rk ]
    uu(:,2,2)   = [ 0._rk, 0._rk, eps, eps, 0._rk, 0._rk, eps, eps ]
    uu(:,3,3)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps, eps, eps, eps ]
    uu(:,1,4)   = [ eps, eps, 0._rk, 0._rk, eps, eps, 0._rk, 0._rk ]
    uu(:,1,5)   = [ eps, eps, eps, eps, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2,6)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps, eps, eps, eps ]

    uu(:,1, 7)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps   ]
    uu(:,2, 8)   = [ eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2, 9)   = [ 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,10)   = [ eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,11)   = [ 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,12)   = [ 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,13)   = [ 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk ]

    uu(:,1,14)   = [ 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,1,15)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk ]
    uu(:,1,16)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk ]
    uu(:,1,17)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps , 0._rk ]

    uu(:,1,18)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps   ]
    uu(:,3,18)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,-2*eps, 0._rk, 0._rk ]

    uu(:,2,19)   = [ eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,19)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -3*eps,0._rk ]

    uu(:,2,20)   = [ 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,20)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -4*eps]

    uu(:,2,21)   = [ 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2,22)   = [ 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2,23)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps , 0._rk, 0._rk, 0._rk ]

    uu(:,2,24)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps , -eps , 0._rk ]

End Function init_displ_hexe8

!------------------------------------------------------------------------------
! FUNCTION: init_displ_hexe20 
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Initialisation of loadcases
!------------------------------------------------------------------------------  
Function init_displ_hexe20(elnodes) Result(uu)

    INTEGER, INTENT(IN) :: elnodes

    REAL(rk), dimension(elnodes,3,elnodes*3) :: uu 
    REAL(rk), PARAMETER :: eps = 1.E-6_rk, z = 0._rk

    uu = z   

    uu(:,1,1)    = [ z  , eps, eps, z  , z  , eps, eps, z    ]
    uu(:,2,2)    = [ z  , z  , eps, eps, z  , z  , eps, eps ]
    uu(:,3,3)    = [ z  , z  , z  , z  , eps, eps, eps, eps ]
    uu(:,1,4)    = [ eps, eps, z  , z  , eps, eps, z  , z    ]
    uu(:,1,5)    = [ eps, eps, eps, eps, z  , z  , z  , z    ]
    uu(:,2,6)    = [ z  , z  , z  , z  , eps, eps, eps, eps ]

    uu(:,1, 7)   = [ z  , z  , z  , z  , z  , z  , z  , eps   ]
    uu(:,2, 8)   = [ eps , z  , z  , z  , z  , z  , z  , z    ]
    uu(:,2, 9)   = [ z  , eps , z  , z  , z  , z  , z  , z    ]
    uu(:,3,10)   = [ eps , z  , z  , z  , z  , z  , z  , z    ]
    uu(:,3,11)   = [ z  , eps , z  , z  , z  , z  , z  , z    ]
    uu(:,3,12)   = [ z  , z  , eps , z  , z  , z  , z  , z    ]
    uu(:,3,13)   = [ z  , z  , z  , eps , z  , z  , z  , z    ]

    uu(:,1,14)   = [ z  , z  , z  , eps , z  , z  , z  , z    ]
    uu(:,1,15)   = [ z  , z  , z  , z  , eps , z  , z  , z    ]
    uu(:,1,16)   = [ z  , z  , z  , z  , z  , eps , z  , z    ]
    uu(:,1,17)   = [ z  , z  , z  , z  , z  , z  , eps , z    ]

    uu(:,1,18)   = [ z  , z  , z  , z  , z  , z  , z  , eps   ]
    uu(:,3,18)   = [ z  , z  , z  , z  , z  ,-2*eps, z  , z    ]

    uu(:,2,19)   = [ eps , z  , z  , z  , z  , z  , z  , z    ]
    uu(:,3,19)   = [ z  , z  , z  , z  , z  , z  , -3*eps,z    ]

    uu(:,2,20)   = [ z  , eps , z  , z  , z  , z  , z  , z    ]
    uu(:,3,20)   = [ z  , z  , z  , z  , z  , z  , z  , -4*eps]

    uu(:,2,21)   = [ z  , z  , eps , z  , z  , z  , z  , z    ]
    uu(:,2,22)   = [ z  , z  , z  , eps , z  , z  , z  , z    ]
    uu(:,2,23)   = [ z  , z  , z  , z  , eps , z  , z  , z    ]

    uu(:,2,24)   = [ z  , z  , z  , z  , z  , eps , -eps , z    ]

End Function init_displ_hexe20


End Module write_deck
