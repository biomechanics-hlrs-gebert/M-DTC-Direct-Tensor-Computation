Module write_deck

  use linFE
  use mesh_partitioning
  use strings
  
  implicit none

  Real(Kind=rk), Parameter     :: delta_b = 1.E-9_rk

  Character(Len=*) , Parameter :: fmt_filename = "(A,I0,A,I0,A)"

contains

  !****************************************************************************
  !**
  !** Determine prescribed displacements
  !**
  Subroutine determine_prescribed_displ(node, xa, xe, elnodes, bnode_vals)

    Real(Kind=rk)   , Intent(In), Dimension(3)  :: node, xa, xe
    Integer         , Intent(in)                :: elnodes

    Real(Kind=rk), dimension(elnodes,3,elnodes*3) :: uu 

    integer                                       :: ii

    Real(Kind=rk)   , Dimension(:,:), Intent(Out) :: bnode_vals

    !**************************************************************************

    !** Loadcase init *********************************************************  
    uu = init_displ(elnodes)

    Do ii = 1, elnodes*3

       bnode_vals(ii,1) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,1,ii))
       bnode_vals(ii,2) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,2,ii))
       bnode_vals(ii,3) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,3,ii))

    End Do

  End Subroutine determine_prescribed_displ

  !****************************************************************************
  !** Initialisation of loadcases                                 
  !**
  Function init_displ(elnodes) Result(uu)

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
    uu(:,2, 8)   = [ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2, 9)   = [ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,10)   = [ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,11)   = [ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,12)   = [ 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,13)   = [ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk ]

    uu(:,1,14)   = [ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,1,15)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk ]
    uu(:,1,16)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk ]
    uu(:,1,17)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk ]

    uu(:,1,18)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps   ]
    uu(:,3,18)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,-2*eps, 0._rk, 0._rk ]

    uu(:,2,19)   = [ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,19)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -3*eps,0._rk ]

    uu(:,2,20)   = [ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,3,20)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -4*eps]

    uu(:,2,21)   = [ 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2,22)   = [ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk ]
    uu(:,2,23)   = [ 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk ]

    uu(:,2,24)   = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , -eps , 0._rk ]

  End Function init_displ

End Module write_deck
