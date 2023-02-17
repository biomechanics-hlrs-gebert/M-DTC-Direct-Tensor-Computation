import numpy as np

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

    uu[:,0,0]   = np.array([ z, eps, eps, z, z, eps, eps, z ])
    uu[:,1,1]   = np.array([ z, z, eps, eps, z, z, eps, eps ])
    uu[:,2,2]   = np.array([ z, z, z, z, eps, eps, eps, eps ])
    uu[:,0,3]   = np.array([ eps, eps, z, z, eps, eps, z, z ])
    uu[:,0,4]   = np.array([ eps, eps, eps, eps, z, z, z, z ])
    uu[:,1,5]   = np.array([ z, z, z, z, eps, eps, eps, eps ])

    uu[:,0,6]   = np.array([ z, z, z, z, z, z, z, eps   ])
    uu[:,1,7]   = np.array([ eps , z, z, z, z, z, z, z ])
    uu[:,1,8]   = np.array([ z, eps , z, z, z, z, z, z ])
    uu[:,2,9]   = np.array([ eps , z, z, z, z, z, z, z ])
    uu[:,2,10]   = np.array([ z, eps , z, z, z, z, z, z ])
    uu[:,2,11]   = np.array([ z, z, eps , z, z, z, z, z ])
    uu[:,2,12]   = np.array([ z, z, z, eps , z, z, z, z ])

    uu[:,0,13]   = np.array([ z, z, z, eps , z, z, z, z ])
    uu[:,0,14]   = np.array([ z, z, z, z, eps , z, z, z ])
    uu[:,0,15]   = np.array([ z, z, z, z, z, eps , z, z ])
    uu[:,0,16]   = np.array([ z, z, z, z, z, z, eps , z ])

    uu[:,0,17]   = np.array([ z, z, z, z, z, z, z, eps   ])
    uu[:,2,17]   = np.array([ z, z, z, z, z,-2*eps, z, z ])

    uu[:,1,18]   = np.array([ eps , z, z, z, z, z, z, z ])
    uu[:,2,18]   = np.array([ z, z, z, z, z, z, -3*eps,z ])

    uu[:,1,19]   = np.array([ z, eps , z, z, z, z, z, z ])
    uu[:,2,19]   = np.array([ z, z, z, z, z, z, z, -4*eps])

    uu[:,1,20]   = np.array([ z, z, eps , z, z, z, z, z ])
    uu[:,1,21]   = np.array([ z, z, z, eps , z, z, z, z ])
    uu[:,1,22]   = np.array([ z, z, z, z, eps , z, z, z ])

    uu[:,1,23]   = np.array([ z, z, z, z, z, eps , -eps , z ])

    return(uu)


def init_displ_hexe20() :

    # REAL(rk), dimension(elnodes,3,elnodes*3) :: uu 

    uu = np.zeros((20,3,60))
    z = 0.
    eps = 0.000001

  # uu[ Knoten (1-20), dof (1-3), Load Case (1-60)]
  # uu[:,1,1]    = np.array([      1,      2,      3,      4,      5,      6,      7,      8,           9,     10,     11,     12,     13,     14,     15,     16,     17,     18,     19,     20 ])

    uu[:,0,0]    = np.array([      z,    eps,    eps,      z,      z,    eps,    eps,      z,      .5*eps,    eps, .5*eps,      z, .5*eps,    eps, .5*eps,      z,      z,    eps,    eps,      z ]) # Load Case 1 | X dir | Normal Stress
    uu[:,1,1]    = np.array([      z,      z,    eps,    eps,      z,      z,    eps,    eps,           z, .5*eps,    eps, .5*eps,      z, .5*eps,    eps, .5*eps,      z,      z,    eps,    eps ]) # Load Case 2 | Y dir | Normal Stress
    uu[:,2,2]    = np.array([      z,      z,      z,      z,    eps,    eps,    eps,    eps,           z,      z,      z,      z,    eps,    eps,    eps,    eps, .5*eps, .5*eps, .5*eps, .5*eps ]) # Load Case 3 | Z dir | Normal Stress

    uu[:,0,3]    = np.array([    eps,    eps,      z,      z,    eps,    eps,      z,      z,         eps, .5*eps,      z, .5*eps,    eps, .5*eps,      z, .5*eps,    eps,    eps,      z,      z ]) # Load Case 4 | Y in X dir | Shear Stress
    uu[:,0,4]    = np.array([    eps,    eps,    eps,    eps,      z,      z,      z,      z,         eps,    eps,    eps,    eps,      z,      z,      z,      z, .5*eps, .5*eps, .5*eps, .5*eps ]) # Load Case 5 | Z in X dir | Shear Stress
    uu[:,1,5]    = np.array([      z,      z,      z,      z,    eps,    eps,    eps,    eps,           z,      z,      z,      z,    eps,    eps,    eps,    eps, .5*eps, .5*eps, .5*eps, .5*eps ]) # Load Case 6 | Z in Y dir | Shear Stress
    
    uu[:,0,6]    = np.array([      z,      z,      z,      z,      z,      z,      z,    eps,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,1,7]    = np.array([    eps,      z,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,1,8]    = np.array([      z,    eps,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,9]    = np.array([    eps,      z,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,10]   = np.array([      z,    eps,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,11]   = np.array([      z,      z,    eps,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,12]   = np.array([      z,      z,      z,    eps,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,0,13]   = np.array([      z,      z,      z,    eps,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,0,14]   = np.array([      z,      z,      z,      z,    eps,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,0,15]   = np.array([      z,      z,      z,      z,      z,    eps,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,0,16]   = np.array([      z,      z,      z,      z,      z,      z,    eps,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,0,17]   = np.array([      z,      z,      z,      z,      z,      z,      z,    eps,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,17]   = np.array([      z,      z,      z,      z,      z, -2*eps,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,1,18]   = np.array([    eps,      z,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,18]   = np.array([      z,      z,      z,      z,      z,      z, -3*eps,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,1,19]   = np.array([      z,    eps,      z,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,2,19]   = np.array([      z,      z,      z,      z,      z,      z,      z, -4*eps,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,1,20]   = np.array([      z,      z,    eps,      z,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,1,21]   = np.array([      z,      z,      z,    eps,      z,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])
    uu[:,1,22]   = np.array([      z,      z,      z,      z,    eps,      z,      z,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])

    uu[:,1,23]   = np.array([      z,      z,      z,      z,      z,    eps,   -eps,      z,           z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z,      z ])


    # uu[:,1,1]    = np.array([ z  , eps, eps, z  , z  , eps, eps, z    ])
    # uu[:,2,2]    = np.array([ z  , z  , eps, eps, z  , z  , eps, eps ])
    # uu[:,3,3]    = np.array([ z  , z  , z  , z  , eps, eps, eps, eps ])
    # uu[:,1,4]    = np.array([ eps, eps, z  , z  , eps, eps, z  , z    ])
    # uu[:,1,5]    = np.array([ eps, eps, eps, eps, z  , z  , z  , z    ])
    # uu[:,2,6]    = np.array([ z  , z  , z  , z  , eps, eps, eps, eps ])

    # uu[:,1, 7]   = np.array([ z  , z  , z  , z  , z  , z  , z  , eps   ])
    # uu[:,2, 8]   = np.array([ eps , z  , z  , z  , z  , z  , z  , z    ])
    # uu[:,2, 9]   = np.array([ z  , eps , z  , z  , z  , z  , z  , z    ])
    # uu[:,3,10]   = np.array([ eps , z  , z  , z  , z  , z  , z  , z    ])
    # uu[:,3,11]   = np.array([ z  , eps , z  , z  , z  , z  , z  , z    ])
    # uu[:,3,12]   = np.array([ z  , z  , eps , z  , z  , z  , z  , z    ])
    # uu[:,3,13]   = np.array([ z  , z  , z  , eps , z  , z  , z  , z    ])

    # uu[:,1,14]   = np.array([ z  , z  , z  , eps , z  , z  , z  , z    ])
    # uu[:,1,15]   = np.array([ z  , z  , z  , z  , eps , z  , z  , z    ])
    # uu[:,1,16]   = np.array([ z  , z  , z  , z  , z  , eps , z  , z    ])
    # uu[:,1,17]   = np.array([ z  , z  , z  , z  , z  , z  , eps , z    ])

    # uu[:,1,18]   = np.array([ z  , z  , z  , z  , z  , z  , z  , eps   ])
    # uu[:,3,18]   = np.array([ z  , z  , z  , z  , z  ,-2*eps, z  , z    ])

    # uu[:,2,19]   = np.array([ eps , z  , z  , z  , z  , z  , z  , z    ])
    # uu[:,3,19]   = np.array([ z  , z  , z  , z  , z  , z  , -3*eps,z    ])

    # uu[:,2,20]   = np.array([ z  , eps , z  , z  , z  , z  , z  , z    ])
    # uu[:,3,20]   = np.array([ z  , z  , z  , z  , z  , z  , z  , -4*eps])

    # uu[:,2,21]   = np.array([ z  , z  , eps , z  , z  , z  , z  , z    ])
    # uu[:,2,22]   = np.array([ z  , z  , z  , eps , z  , z  , z  , z    ])
    # uu[:,2,23]   = np.array([ z  , z  , z  , z  , eps , z  , z  , z    ])

    # uu[:,2,24]   = np.array([ z  , z  , z  , z  , z  , eps , -eps , z    ])

    return(uu)
    
uu = init_displ_hexe8() 
Rank = np.linalg.matrix_rank(uu, tol=None, hermitian=False)
print("Rank of HEXE8 displacements:", Rank)
SRank = np.sum(Rank)
print("Sum of Ranks of HEXE8 displacements:", SRank)


uu = init_displ_hexe20() 
Rank = np.linalg.matrix_rank(uu, tol=None, hermitian=False)
print("Rank of HEXE20 displacements:", Rank)
SRank = np.sum(Rank)
print("Sum of Ranks of HEXE20 displacements:", SRank)
