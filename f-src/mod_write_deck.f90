Module write_deck

use linFE
use mesh_partitioning
use strings

implicit none

Real(rk), Parameter :: delta_b = 1.E-9_rk
Character(*) , Parameter :: fmt_filename = "(A,I0,A,I0,A)"

CONTAINS

!------------------------------------------------------------------------------  
! FUNCTION: determine_prescribed_displ
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Determine prescribed displacements
!------------------------------------------------------------------------------  
Subroutine determine_prescribed_displ(node, xa, xe, elnodes, bnode_vals)
    
    REAL(rk), INTENT(IN), DIMENSION(3) :: node, xa, xe
    INTEGER(ik), INTENT(IN) :: elnodes
    REAL(rk), DIMENSION(:,:), INTENT(OUT) :: bnode_vals

    REAL(rk), DIMENSION(elnodes, 3 ,elnodes*3) :: uu 
    INTEGER(ik) :: ii

    ! Loadcase init

    IF (elnodes == 8) THEN
        uu = init_displ_hexe8(elnodes)

        Do ii = 1, elnodes*3
            bnode_vals(ii,1) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,1,ii))
            bnode_vals(ii,2) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,2,ii))
            bnode_vals(ii,3) = sum(phi_NN_hexe8(t_geom_xi(node,xa,xe))*uu(:,3,ii))
         End Do

    ELSE IF (elnodes == 20) THEN
        uu = init_displ_hexe20(elnodes)

        Do ii = 1, elnodes*3
            bnode_vals(ii,1) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,1,ii))
            bnode_vals(ii,2) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,2,ii))
            bnode_vals(ii,3) = sum(phi_NN_hexe20(t_geom_xi(node,xa,xe))*uu(:,3,ii))
         End Do

    END IF 

End Subroutine determine_prescribed_displ

!------------------------------------------------------------------------------
! FUNCTION: init_displ_hexe8 
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Initialisation of loadcases
!------------------------------------------------------------------------------  
Function init_displ_hexe8(elnodes) Result(uu)

    Integer, Intent(in) :: elnodes

    Real(rk), dimension(elnodes,3,elnodes*3) :: uu 
    Real(rk) :: eps=1.E-6_rk
 
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
    REAL(rk), PARAMETER :: eps = 1.E-6_rk, z = 0._rk, f = 0.99794_rk, oh=1.5_rk, h=0.5_rk
    INTEGER(ik), PARAMETER :: ish = 1_ik ! Index shift for better compatibility with, e.g., Python

    uu  = z

    uu(:,ish+0,ish+0)    = [       z,     eps,     eps,       z,       z,     eps,     eps,       z,       .5_rk*eps,&
         eps,  .5_rk*eps,       z,  .5_rk*eps,     eps,  .5_rk*eps,       z,       z,     eps,     eps,       z ] ! Load Case 1 | X dir | Normal Stress
    uu(:,ish+1,ish+1)    = [       z,       z,     eps,     eps,       z,       z,     eps,     eps,            z,&
      .5_rk*eps,     eps,  .5_rk*eps,       z,  .5_rk*eps,     eps,  .5_rk*eps,       z,       z,     eps,     eps ] ! Load Case 2 | Y dir | Normal Stress
    uu(:,ish+2,ish+2)    = [       z,       z,       z,       z,     eps,     eps,     eps,     eps,            z,&
           z,       z,       z,     eps,     eps,     eps,     eps,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps ] ! Load Case 3 | Z dir | Normal Stress
    uu(:,ish+0,ish+3)    = [     eps,     eps,       z,       z,     eps,     eps,       z,       z,          eps,&
      .5_rk*eps,       z,  .5_rk*eps,     eps,  .5_rk*eps,       z,  .5_rk*eps,     eps,     eps,       z,       z ] ! Load Case 4 | Y in X dir | Shear Stress
    uu(:,ish+0,ish+4)    = [     eps,     eps,     eps,     eps,       z,       z,       z,       z,          eps,&
         eps,     eps,     eps,       z,       z,       z,       z,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps ] ! Load Case 5 | Z in X dir | Shear Stress
    uu(:,ish+1,ish+5)    = [       z,       z,       z,       z,     eps,     eps,     eps,     eps,            z,&
           z,       z,       z,     eps,     eps,     eps,     eps,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps,  .5_rk*eps ] ! Load Case 6 | Z in Y dir | Shear Stress
    uu(:,ish+0,ish+6)    = [       z,       z,       z,       z,       z,       z,       z,     eps,            z,&
           z,       z,       z,       z,       z,  .5_rk*eps,       z,       z,       z,       z,  .5_rk*eps ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+7)    = [  -2._rk*eps,       z,       z,       z,       z,       z,       z,       z,       .5_rk*eps,&
           z,       z,       z,       z,       z,       z,       z,  .5_rk*eps,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+8)    = [       z,    -eps,       z,       z,       z,       z,       z,       z,       .5_rk*eps,&
           z,       z,       z,       z,       z,       z,       z,       z,  .5_rk*eps,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+9)    = [    -eps,       z,       z,       z,       z,       z,       z,       z,       .5_rk*eps,&
           z,       z,  .1_rk*eps,       z,       z,       z,       z,  .5_rk*eps,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+10)   = [       z,  -2._rk*eps,       z,       z,       z,       z,       z,       z,            z,&
      .5_rk*eps,       z,       z,       z,       z,       z,       z,       z,  .5_rk*eps,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+11)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,  .1_rk*eps,       z,       z,       z,       z,       z,       z,       z,  .9*eps,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+12)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,   5._rk*eps ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+0,ish+13)   = [       z,       z,       z,  .2_rk*eps,       z,       z,       z,       z,            z,&
           z,       z,  .5_rk*eps,       z,       z,       z,       z,       z,       z,       z,   5._rk*eps ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+0,ish+14)   = [       z,       z,       z,       z,     eps,       z,       z,       z,            z,&
           z,       z,       z,  .5_rk*eps,       z,       z,       z,       z,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+0,ish+15)   = [       z,       z,       z,       z,       z,   5._rk*eps,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,  .2_rk*eps,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+0,ish+16)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,  .5_rk*eps,       z,       z,       z,       z,  .3_rk*eps,       z ] ! Other std. load cases, inherited by HEXE8  
    uu(:,ish+0,ish+17)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,  .1_rk*eps,       z,       z,       z,       z,  .1_rk*eps ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+17)   = [       z,       z,       z,       z,       z, -.7*eps,       z,       z,            z,&
           z,       z,       z, -.1_rk*eps,    -eps,       z,       z,       z,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+18)   = [       z,       z,       z,       z,       z,       z,       z,       z,       .5_rk*eps,&
           z,       z,  .5_rk*eps,       z,       z,       z,       z,  .5_rk*eps,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+18)   = [       z,       z,       z,       z,       z,       z,  -3._rk*eps,       z,            z,&
           z,       z,       z,       z,-1.5_rk*eps,       z,       z,       z,       z,-1.5_rk*eps,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+19)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,  .1_rk*eps,       z,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+2,ish+19)   = [       z,       z,       z,       z,       z,       z,       z,  -3._rk*eps,            z,&
           z,       z,       z,       z,       z, -.1_rk*eps,       z,       z,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+20)   = [       z,       z,   5._rk*eps,       z,       z,       z,       z,       z,            z,&
      .5_rk*eps,       z,       z,       z,       z,       z,       z,       z,       z,  .5_rk*eps,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+21)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,   5._rk*eps,       z,       z,       z,       z,       z,       z,       z,       z,  .5_rk*eps ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+22)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,  2._rk*eps,       z,       z,       z,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+23)   = [       z,       z,       z,       z,       z,     z,         z,       z,            z,&
           z,       z,       z,   5._rk*eps,       z,       z,       z,       z,       z, -.5_rk*eps,       z ] ! Other std. load cases, inherited by HEXE8
    uu(:,ish+1,ish+24)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  1 | tLC 25
    uu(:,ish+2,ish+24)   = [       z,       z,       z,       z,       z,       z,       z,       z,        2._rk*eps,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  1 | tLC 25
    uu(:,ish+0,ish+25)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  2 | tLC 26
    uu(:,ish+2,ish+25)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
      5._rk* eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  2 | tLC 26
    uu(:,ish+1,ish+26)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,-.1_rk*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  3 | tLC 27
    uu(:,ish+2,ish+26)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,    eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  3 | tLC 27
    uu(:,ish+0,ish+27)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,        z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  4 | tLC 28
    uu(:,ish+2,ish+27)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,   4._rk*eps,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  4 | tLC 28
    uu(:,ish+1,ish+28)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,  .1_rk*eps,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  5 | tLC 29
    uu(:,ish+2,ish+28)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,  -4._rk*eps,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  5 | tLC 29
    uu(:,ish+0,ish+29)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z, -.1_rk*eps,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  6 | tLC 30
    uu(:,ish+2,ish+29)   = [       z,       z,       z,       z,       z,   2._rk*eps,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  6 | tLC 30
    uu(:,ish+1,ish+30)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z, .1_rk*-eps,       z,       z,       z,       z,       z ] ! HEXE20 iLC  7 | tLC 31
    uu(:,ish+2,ish+30)   = [       z,       z,   5._rk*eps,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  7 | tLC 31
    uu(:,ish+0,ish+31)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  8 | tLC 32
    uu(:,ish+2,ish+31)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z, -.5_rk*eps,       z,       z,       z,       z ] ! HEXE20 iLC  8 | tLC 32
    uu(:,ish+0,ish+32)   = [       z,       z,       z,       z,       z,  .1_rk*eps,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,  -1*eps,       z,       z,       z ] ! HEXE20 iLC  9 | tLC 33
    uu(:,ish+1,ish+32)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC  9 | tLC 33
    uu(:,ish+0,ish+33)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 10 | tLC 34
    uu(:,ish+1,ish+33)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,  -4._rk*eps,       z,       z ] ! HEXE20 iLC 10 | tLC 34
    uu(:,ish+0,ish+34)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 11 | tLC 35
    uu(:,ish+1,ish+34)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z ] ! HEXE20 iLC 11 | tLC 35
    uu(:,ish+0,ish+35)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 12 | tLC 36
    uu(:,ish+1,ish+35)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,   2._rk*eps ] ! HEXE20 iLC 12 | tLC 36
    uu(:,ish+0,ish+36)   = [   3._rk*eps,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 13 | tLC 37 | Single Vertex
    uu(:,ish+1,ish+37)   = [       z,  -5._rk*eps,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 14 | tLC 38 | Single Vertex
    uu(:,ish+0,ish+38)   = [       z,       z,  -3._rk*eps,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 15 | tLC 39 | Single Vertex
    uu(:,ish+1,ish+39)   = [       z,       z,       z,  -2._rk*eps,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 16 | tLC 40 | Single Vertex
    uu(:,ish+0,ish+40)   = [       z,       z,       z,       z,   5._rk*eps,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 17 | tLC 41 | Single Vertex
    uu(:,ish+1,ish+41)   = [       z,       z,       z,       z,       z,   3._rk*eps,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 18 | tLC 42 | Single Vertex
    uu(:,ish+0,ish+42)   = [       z,       z,       z,       z,       z,       z,  -4._rk*eps,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 19 | tLC 43 | Single Vertex
    uu(:,ish+1,ish+43)   = [       z,       z,       z,       z,       z,       z,       z,  4._rk* eps,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 20 | tLC 44 | Single Vertex
    uu(:,ish+1,ish+44)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,   2._rk*eps,       z,       z,       z ] ! HEXE20 iLC 21 | tLC 45 | Compress Center
    uu(:,ish+2,ish+45)   = [       z,       z,       z,  -2._rk*eps,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,-1.5_rk*eps,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 22 | tLC 46 | Compress Center
    uu(:,ish+0,ish+46)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,    -eps,    -eps,  2._rk*eps ] ! HEXE20 iLC 23 | tLC 47 | Compress Center
    uu(:,ish+2,ish+47)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ] ! HEXE20 iLC 24 | tLC 48 | Compress Center
    uu(:,ish+0,ish+48)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
     -.9*eps,       z,       z,       z, -.7*eps,       z,     eps,       z,       z,       z,       z ] ! HEXE20 iLC 25 | tLC 49 | Compress Center
    uu(:,ish+1,ish+49)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,     eps,       z,       z,       z,       z,       z ] ! HEXE20 iLC 26 | tLC 50 | Compress Center
    uu(:,ish+0,ish+50)   = [       z,       z,       z,       z,       z,       z,       z,       z,          eps,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 27 | tLC 51 | Move Center
    uu(:,ish+1,ish+51)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
       3._rk*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 28 | tLC 52 | Move Center
    uu(:,ish+0,ish+52)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 29 | tLC 53 | Move Center
    uu(:,ish+1,ish+53)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,   5._rk*eps,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 30 | tLC 54 | Move Center
    uu(:,ish+0,ish+54)   = [       z,       z,       z,-1.5_rk*eps,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 31 | tLC 55 | Move Center
    uu(:,ish+1,ish+55)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,  -5._rk*eps,       z,       z,       z,       z,       z,       z ] ! HEXE20 iLC 32 | tLC 56 | Move Center
    uu(:,ish+0,ish+56)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ] ! HEXE20 iLC 33 | tLC 57 | Move Center
    uu(:,ish+0,ish+57)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z ] ! HEXE20 iLC 34 | tLC 58 | Move Center
    uu(:,ish+2,ish+58)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ] ! HEXE20 iLC 35 | tLC 59 | Move Center
    uu(:,ish+2,ish+59)   = [       z,       z,       z,       z,       z,       z,       z,       z,            z,&
           z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z ] ! HEXE20 iLC 36 | tLC 60 | Move Center
                 
End Function init_displ_hexe20


End Module write_deck
