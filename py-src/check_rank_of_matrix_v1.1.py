import numpy as np
np.set_printoptions(linewidth=np.inf, threshold=np.inf)

#------------------------------------------------------------------------------
# FUNCTION: assemble_UU 
#------------------------------------------------------------------------------  
#> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
#
#> @brief
#> Assemble the square displacement matrix
#------------------------------------------------------------------------------  
def assemble_UU(uu, nodes, dof):

  size = nodes * dof

  UU = np.zeros((size,size))
  for ii in range(size):
    U_row = np.zeros((size))

    ll = 0
    for kk in range(dof):
      for jj in range(nodes):
          U_row[ll] = uu[jj,kk,ii]
          ll += 1

      UU[:,ii] = U_row

  return(UU)

#------------------------------------------------------------------------------
# FUNCTION: init_displ_hexe20 
#------------------------------------------------------------------------------  
#> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
#
#> @brief
#> Initialisation of loadcases
#------------------------------------------------------------------------------  
def init_displ_hexe8() :

    uu = np.zeros((8,3,24))
    z = 0.
    eps = 0.000001

    uu[:,0,0]   = np.array([ z,  eps,  eps, z, z,  eps,  eps, z ])
    uu[:,1,1]   = np.array([ z, z,  eps,  eps, z, z,  eps,  eps ])
    uu[:,2,2]   = np.array([ z, z, z, z,  eps,  eps,  eps,  eps ])
    uu[:,0,3]   = np.array([  eps,  eps, z, z,  eps,  eps, z, z ])
    uu[:,0,4]   = np.array([  eps,  eps,  eps,  eps, z, z, z, z ])
    uu[:,1,5]   = np.array([ z, z, z, z,  eps,  eps,  eps,  eps ])

    uu[:,0,6]   = np.array([ z, z, z, z, z, z, z,  eps   ])
    uu[:,1,7]   = np.array([  eps , z, z, z, z, z, z, z ])
    uu[:,1,8]   = np.array([ z,  eps , z, z, z, z, z, z ])
    uu[:,2,9]   = np.array([  eps , z, z, z, z, z, z, z ])
    uu[:,2,10]   = np.array([ z,  eps , z, z, z, z, z, z ])
    uu[:,2,11]   = np.array([ z, z,  eps , z, z, z, z, z ])
    uu[:,2,12]   = np.array([ z, z, z,  eps , z, z, z, z ])

    uu[:,0,13]   = np.array([ z, z, z,  eps , z, z, z, z ])
    uu[:,0,14]   = np.array([ z, z, z, z,  eps , z, z, z ])
    uu[:,0,15]   = np.array([ z, z, z, z, z,  eps , z, z ])
    uu[:,0,16]   = np.array([ z, z, z, z, z, z,  eps , z ])

    uu[:,0,17]   = np.array([ z, z, z, z, z, z, z,  eps   ])
    uu[:,2,17]   = np.array([ z, z, z, z, z,-2*eps, z, z ])

    uu[:,1,18]   = np.array([  eps , z, z, z, z, z, z, z ])
    uu[:,2,18]   = np.array([ z, z, z, z, z, z, -3*eps,z ])

    uu[:,1,19]   = np.array([ z,  eps , z, z, z, z, z, z ])
    uu[:,2,19]   = np.array([ z, z, z, z, z, z, z, -4*eps])

    uu[:,1,20]   = np.array([ z, z,  eps , z, z, z, z, z ])
    uu[:,1,21]   = np.array([ z, z, z,  eps , z, z, z, z ])
    uu[:,1,22]   = np.array([ z, z, z, z,  eps , z, z, z ])

    uu[:,1,23]   = np.array([ z, z, z, z, z,  eps ,  -eps , z ])

    return(uu)


def init_displ_hexe20() :

    # REAL(rk), dimension(elnodes,3,elnodes*3) :: uu 

    uu = np.zeros((20,3,60))
    z = 0.
    eps = 0.000001
    f  = 0.99794 
    ish = 0 # indice shift


  # uu[ Knoten (1-20), dof (1-3), Load Case (1-60)]
  # uu[:,1,1]    = np.array([       1,       2,       3,       4,       5,       6,       7,       8,            9,      10,      11,      12,      13,      14,      15,      16,      17,      18,      19,      20 ])

    uu[:,ish+0,ish+0]    = np.array([       z,   2*eps,   2*eps,       z,       z,   2*eps,   2*eps,       z,       .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,       z,     eps,     eps,       z ]) # Load Case 1 | X dir | Normal Stress
    uu[:,ish+1,ish+1]    = np.array([       z,       z,   2*eps,   2*eps,       z,       z,   2*eps,   2*eps,            z,  .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,       z,     eps,     eps ]) # Load Case 2 | Y dir | Normal Stress
    uu[:,ish+2,ish+2]    = np.array([       z,       z,       z,       z,  .5*eps,  .5*eps,  .5*eps,  .5*eps,            z,       z,       z,       z,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 3 | Z dir | Normal Stress

    uu[:,ish+0,ish+3]    = np.array([  .5*eps,  .5*eps,       z,       z,  .5*eps,  .5*eps,       z,       z,          eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,     eps,       z,       z ]) # Load Case 4 | Y in X dir | Shear Stress
    uu[:,ish+0,ish+4]    = np.array([  .5*eps,  .5*eps,  .5*eps,  .5*eps,       z,       z,       z,       z,        2*eps,   2*eps,   2*eps,   2*eps,       z,       z,       z,       z,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 5 | Z in X dir | Shear Stress
    uu[:,ish+1,ish+5]    = np.array([       z,       z,       z,       z,   2*eps,   2*eps,   2*eps,   2*eps,            z,       z,       z,       z,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 6 | Z in Y dir | Shear Stress
    
    uu[:,ish+0,ish+6]    = np.array([       z,       z,       z,       z,       z,       z,       z,     eps,            z,       z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+1,ish+7]    = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+1,ish+8]    = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+9]    = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+10]   = np.array([       z,   2*eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+11]   = np.array([       z,       z,     eps,       z,       z,       z,       z,       z,            z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+12]   = np.array([       z,       z,       z,     eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+0,ish+13]   = np.array([       z,       z,       z,   2*eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,   5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+0,ish+14]   = np.array([       z,       z,       z,       z,   2*eps,       z,       z,       z,            z,       z,       z,       z,  .5*eps,       z,       z,  .5*eps,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+0,ish+15]   = np.array([       z,       z,       z,       z,       z,   2*eps,       z,       z,            z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+0,ish+16]   = np.array([       z,       z,       z,       z,       z,       z,   2*eps,       z,            z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+0,ish+17]   = np.array([       z,       z,       z,       z,       z,       z,       z,     eps,            z,       z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+17]   = np.array([       z,       z,       z,       z,       z,  -2*eps,       z,       z,            z,       z,       z,       z,    -eps,    -eps,       z,       z,       z,    -eps,       z,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+1,ish+18]   = np.array([   2*eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+18]   = np.array([       z,       z,       z,       z,       z,       z,  -3*eps,       z,            z,       z,       z,       z,       z,-1.5*eps,-1.5*eps,       z,       z,       z,-1.5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+1,ish+19]   = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,  .5*eps,       z,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+2,ish+19]   = np.array([       z,       z,       z,       z,       z,       z,       z,  -4*eps,            z,       z,       z,       z,       z,       z,  -2*eps,  -2*eps,       z,       z,       z,  -2*eps ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+1,ish+20]   = np.array([       z,       z,   2*eps,       z,       z,       z,       z,       z,            z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+1,ish+21]   = np.array([       z,       z,       z,     eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,ish+1,ish+22]   = np.array([       z,       z,       z,       z,     eps,       z,       z,       z,            z,       z,       z,       z,  .5*eps,       z,       z,  .5*eps,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,ish+1,ish+23]   = np.array([       z,       z,       z,       z,       z,     eps,    -eps,       z,            z,       z,       z,       z,  .5*eps, -.5*eps,       z,       z,       z,  .5*eps, -.5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    # Node 9                                                                                                                                                                                                             # independent Load Case | total Load Case
    uu[:,ish+1,ish+24]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  1 | tLC 25
    uu[:,ish+2,ish+24]   = np.array([       z,  -4*eps,       z,       z,       z,       z,       z,       z,        2*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  1 | tLC 25
    # Node 10
    uu[:,ish+0,ish+25]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  2 | tLC 26
    uu[:,ish+2,ish+25]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,   2*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  2 | tLC 26
    # Node 11
    uu[:,ish+1,ish+26]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  3 | tLC 27
    uu[:,ish+2,ish+26]   = np.array([       z,       z,     eps,       z,       z,       z,       z,       z,            z,       z,   2*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  3 | tLC 27
    # Node 12
    uu[:,ish+0,ish+27]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,  .5*eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  4 | tLC 28
    uu[:,ish+2,ish+27]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,   2*eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  4 | tLC 28
    # Node 13
    uu[:,ish+1,ish+28]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,  .5*eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  5 | tLC 29
    uu[:,ish+2,ish+28]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,  -2*eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  5 | tLC 29
    # Node 14
    uu[:,ish+0,ish+29]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  6 | tLC 30
    uu[:,ish+2,ish+29]   = np.array([       z,       z,       z,       z,       z,   2*eps,       z,       z,            z,       z,       z,       z,       z,  -2*eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  6 | tLC 30
    # Node 15
    uu[:,ish+1,ish+30]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC  7 | tLC 31
    uu[:,ish+2,ish+30]   = np.array([       z,       z,   2*eps,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,  -2*eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC  7 | tLC 31
    # Node 16
    uu[:,ish+0,ish+31]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,  .5*eps,       z,       z,       z,       z ]) # HEXE20 iLC  8 | tLC 32
    uu[:,ish+2,ish+31]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC  8 | tLC 32
    # Node 17
    uu[:,ish+0,ish+32]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,   2*eps,       z,       z,       z ]) # HEXE20 iLC  9 | tLC 33
    uu[:,ish+1,ish+32]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ]) # HEXE20 iLC  9 | tLC 33
    # Node 18
    uu[:,ish+0,ish+33]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z ]) # HEXE20 iLC 10 | tLC 34
    uu[:,ish+1,ish+33]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,   2*eps,       z,       z ]) # HEXE20 iLC 10 | tLC 34
    # Node 19
    uu[:,ish+0,ish+34]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z ]) # HEXE20 iLC 11 | tLC 35
    uu[:,ish+1,ish+34]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z ]) # HEXE20 iLC 11 | tLC 35
    # Node 20
    uu[:,ish+0,ish+35]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,     eps ]) # HEXE20 iLC 12 | tLC 36
    uu[:,ish+1,ish+35]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps ]) # HEXE20 iLC 12 | tLC 36

    uu[:,ish+0,ish+36]   = np.array([   3*eps,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 13 | tLC 37 | Single Vertex
    uu[:,ish+1,ish+37]   = np.array([       z,   3*eps,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 14 | tLC 38 | Single Vertex
    uu[:,ish+0,ish+38]   = np.array([       z,       z,  -3*eps,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 15 | tLC 39 | Single Vertex
    uu[:,ish+1,ish+39]   = np.array([       z,       z,       z,  -3*eps,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 16 | tLC 40 | Single Vertex
    uu[:,ish+0,ish+40]   = np.array([       z,       z,       z,       z,   3*eps,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 17 | tLC 41 | Single Vertex
    uu[:,ish+1,ish+41]   = np.array([       z,       z,       z,       z,       z,   3*eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 18 | tLC 42 | Single Vertex
    uu[:,ish+0,ish+42]   = np.array([       z,       z,       z,       z,       z,       z,  2*-eps,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 19 | tLC 43 | Single Vertex
    uu[:,ish+1,ish+43]   = np.array([       z,       z,       z,       z,       z,       z,       z,  2*-eps,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 20 | tLC 44 | Single Vertex

    uu[:,ish+1,ish+44]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,   2*eps,  .5*eps,    -eps,    -eps ]) # HEXE20 iLC 21 | tLC 45 | Compress Center
    uu[:,ish+2,ish+45]   = np.array([       z,       z,       z,   3*eps,       z,   2*eps,       z,       z,            z,   2*eps,       z,     eps,       z,  -2*eps,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC 22 | tLC 46 | Compress Center
    uu[:,ish+0,ish+46]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,  .5*eps,    -eps,    -eps,   2*eps ]) # HEXE20 iLC 23 | tLC 47 | Compress Center
    uu[:,ish+2,ish+47]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,  .5*eps,       z,    -eps,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 24 | tLC 48 | Compress Center
    uu[:,ish+0,ish+48]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,    -eps,       z,   2*eps,       z,  -2*eps,       z,     eps,       z,       z,       z,       z ]) # HEXE20 iLC 25 | tLC 49 | Compress Center
    uu[:,ish+1,ish+49]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,    -eps,       z,     eps,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 26 | tLC 50 | Compress Center

    uu[:,ish+0,ish+50]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,        2*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 27 | tLC 51 | Move Center
    uu[:,ish+1,ish+51]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,   2*eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 28 | tLC 52 | Move Center
    uu[:,ish+0,ish+52]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 29 | tLC 53 | Move Center
    uu[:,ish+1,ish+53]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 30 | tLC 54 | Move Center
    uu[:,ish+0,ish+54]   = np.array([       z,       z,       z,   2*eps,       z,       z,       z,       z,            z,       z,       z,       z,     eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 31 | tLC 55 | Move Center
    uu[:,ish+1,ish+55]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,  3* eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 32 | tLC 56 | Move Center
    uu[:,ish+0,ish+56]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 33 | tLC 57 | Move Center
    uu[:,ish+0,ish+57]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC 34 | tLC 58 | Move Center
    uu[:,ish+2,ish+58]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ]) # HEXE20 iLC 35 | tLC 59 | Move Center
    uu[:,ish+2,ish+59]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z ]) # HEXE20 iLC 36 | tLC 60 | Move Center

    return(uu)
    
#------------------------------------------------------------------------------  
# Get HEXE8 information
#------------------------------------------------------------------------------  
uu = init_displ_hexe8() 
UU = assemble_UU(uu, 8, 3)

Rank = np.linalg.matrix_rank(UU, tol=None, hermitian=False)

#------------------------------------------------------------------------------  
# Print HEXE8 displacements
#------------------------------------------------------------------------------  
UU *= 1000000
print("Scaled UU:")
# print(np.around(UU.astype(int)))

#------------------------------------------------------------------------------  
# Diagonalize HEXE8 UU
#------------------------------------------------------------------------------  
# DD = np.diag(np.diag(UU))
# print("Diagonalized:")
# print(np.around(DD.astype(int)))

#------------------------------------------------------------------------------  
# Get determinant
#------------------------------------------------------------------------------  
det = np.linalg.det(UU)

print("Rank        of HEXE8 displacements:", Rank)
print("Determinant of HEXE8 displacements:", np.around(det,3))


#------------------------------------------------------------------------------  
# Get Condition
#------------------------------------------------------------------------------  
cond = np.linalg.cond(UU)
print("Condition   of HEXE8 displacements:", np.around(cond,3))

print('\n')





#------------------------------------------------------------------------------  
# Get HEXE20 information
#------------------------------------------------------------------------------  
uu = init_displ_hexe20() 
UU = assemble_UU(uu, 20, 3)

Rank = np.linalg.matrix_rank(UU, tol=None, hermitian=False)

#------------------------------------------------------------------------------  
# Print HEXE20 displacements
#------------------------------------------------------------------------------  
UU *= 1000000
print("Scaled UU:")
# print(np.around(UU,12))

#------------------------------------------------------------------------------  
# Diagonalize HEXE20 UU
#------------------------------------------------------------------------------  
# DD = np.diag(np.diag(UU))
# print("Diagonalized:")
# print(DD)

#------------------------------------------------------------------------------  
# Get determinant
#------------------------------------------------------------------------------  
det = np.linalg.det(UU)

print("Rank        of HEXE20 displacements:", Rank)
print("Determinant of HEXE20 displacements:", np.around(det,3))

#------------------------------------------------------------------------------  
# Get Condition
#------------------------------------------------------------------------------  
cond = np.linalg.cond(UU)
print("Condition   of HEXE20 displacements:", np.around(cond,3))
