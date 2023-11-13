   uu[:,0,0]    = np.array([       z,     eps,     eps,       z,       z,     eps,     eps,       z,       .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,       z,     eps,     eps,       z ]) # Load Case 1 | X dir | Normal Stress
    uu[:,1,1]    = np.array([       z,       z,     eps,     eps,       z,       z,     eps,     eps,            z,  .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,       z,     eps,     eps ]) # Load Case 2 | Y dir | Normal Stress
    uu[:,2,2]    = np.array([       z,       z,       z,       z,     eps,     eps,     eps,     eps,            z,       z,       z,       z,     eps,     eps,     eps,     eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 3 | Z dir | Normal Stress

    uu[:,0,3]    = np.array([     eps,     eps,       z,       z,     eps,     eps,       z,       z,          eps,  .5*eps,       z,  .5*eps,     eps,  .5*eps,       z,  .5*eps,     eps,     eps,       z,       z ]) # Load Case 4 | Y in X dir | Shear Stress
    uu[:,0,4]    = np.array([     eps,     eps,     eps,     eps,       z,       z,       z,       z,          eps,     eps,     eps,     eps,       z,       z,       z,       z,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 5 | Z in X dir | Shear Stress
    uu[:,1,5]    = np.array([       z,       z,       z,       z,     eps,     eps,     eps,     eps,            z,       z,       z,       z,     eps,     eps,     eps,     eps,  .5*eps,  .5*eps,  .5*eps,  .5*eps ]) # Load Case 6 | Z in Y dir | Shear Stress
    
    uu[:,0,6]    = np.array([       z,       z,       z,       z,       z,       z,       z,     eps,            z,       z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,1,7]    = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,1,8]    = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,9]    = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,10]   = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,11]   = np.array([       z,       z,     eps,       z,       z,       z,       z,       z,            z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,12]   = np.array([       z,       z,       z,     eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8

    uu[:,0,13]   = np.array([       z,       z,       z,     eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,   5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,0,14]   = np.array([       z,       z,       z,       z,     eps,       z,       z,       z,            z,       z,       z,       z,  .5*eps,       z,       z,  .5*eps,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,0,15]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,0,16]   = np.array([       z,       z,       z,       z,       z,       z,     eps,       z,            z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,0,17]   = np.array([       z,       z,       z,       z,       z,       z,       z,     eps,            z,       z,       z,       z,       z,       z,  .5*eps,  .5*eps,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,17]   = np.array([       z,       z,       z,       z,       z,  -2*eps,       z,       z,            z,       z,       z,       z,    -eps,    -eps,       z,       z,       z,    -eps,       z,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,1,18]   = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,       .5*eps,       z,       z,  .5*eps,       z,       z,       z,       z,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,18]   = np.array([       z,       z,       z,       z,       z,       z,  -3*eps,       z,            z,       z,       z,       z,       z,-1.5*eps,-1.5*eps,       z,       z,       z,-1.5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,1,19]   = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,       .5*eps,  .5*eps,       z,       z,       z,       z,       z,  .5*eps,       z,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,2,19]   = np.array([       z,       z,       z,       z,       z,       z,       z,  -4*eps,            z,       z,       z,       z,       z,       z,  -2*eps,  -2*eps,       z,       z,       z,  -2*eps ]) # Other std. load cases, inherited by HEXE8

    uu[:,1,20]   = np.array([       z,       z,     eps,       z,       z,       z,       z,       z,            z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps,       z ]) # Other std. load cases, inherited by HEXE8
    uu[:,1,21]   = np.array([       z,       z,       z,     eps,       z,       z,       z,       z,            z,       z,  .5*eps,  .5*eps,       z,       z,       z,       z,       z,       z,       z,  .5*eps ]) # Other std. load cases, inherited by HEXE8
    uu[:,1,22]   = np.array([       z,       z,       z,       z,     eps,       z,       z,       z,            z,       z,       z,       z,  .5*eps,       z,       z,  .5*eps,  .5*eps,       z,       z,       z ]) # Other std. load cases, inherited by HEXE8

    uu[:,1,23]   = np.array([       z,       z,       z,       z,       z,     eps,    -eps,       z,            z,       z,       z,       z,  .5*eps, -.5*eps,       z,       z,       z,  .5*eps, -.5*eps,       z ]) # Other std. load cases, inherited by HEXE8

    # Node 9                                                                                                                                                                                                             # independent Load Case | total Load Case
    uu[:,1,24]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,          eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  1 | tLC 25
    uu[:,2,24]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,          eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  1 | tLC 25
    # Node 10
    uu[:,0,25]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  2 | tLC 26
    uu[:,2,25]   = np.array([       z,    -eps,     eps,       z,       z,       z,       z,       z,            z,     eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  2 | tLC 26
    # Node 11
    uu[:,1,26]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  3 | tLC 27
    uu[:,2,26]   = np.array([       z,    -eps,     eps,       z,       z,       z,       z,       z,            z,       z,     eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  3 | tLC 27
    # Node 12
    uu[:,0,27]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,     eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  4 | tLC 28
    uu[:,2,27]   = np.array([       z,     eps,   3*eps,       z,       z,       z,       z,       z,            z,       z,       z,     eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  4 | tLC 28
    # Node 13
    uu[:,1,28]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,     eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  5 | tLC 29
    uu[:,2,28]   = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  5 | tLC 29
    # Node 14
    uu[:,0,29]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  6 | tLC 30
    uu[:,2,29]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC  6 | tLC 30
    # Node 15
    uu[:,1,30]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC  7 | tLC 31
    uu[:,2,30]   = np.array([       z,       z,       z,       z,       z,       z,    -eps,   2*eps,            z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC  7 | tLC 31
    # Node 16
    uu[:,0,31]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z,       z ]) # HEXE20 iLC  8 | tLC 32
    uu[:,2,31]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC  8 | tLC 32
    # Node 17
    uu[:,0,32]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ]) # HEXE20 iLC  9 | tLC 33
    uu[:,1,32]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ]) # HEXE20 iLC  9 | tLC 33
    # Node 18
    uu[:,0,33]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z ]) # HEXE20 iLC 10 | tLC 34
    uu[:,1,33]   = np.array([       z,       z,       z,       z,  -2*eps,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z ]) # HEXE20 iLC 10 | tLC 34
    # Node 19
    uu[:,0,34]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z ]) # HEXE20 iLC 11 | tLC 35
    uu[:,1,34]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z ]) # HEXE20 iLC 11 | tLC 35
    # Node 20
    uu[:,0,35]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,     eps ]) # HEXE20 iLC 12 | tLC 36
    uu[:,1,35]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps ]) # HEXE20 iLC 12 | tLC 36

    uu[:,0,36]   = np.array([     eps,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 13 | tLC 37 | Single Vertex
    uu[:,1,37]   = np.array([       z,     eps,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 14 | tLC 38 | Single Vertex
    uu[:,0,38]   = np.array([       z,       z,    -eps,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 15 | tLC 39 | Single Vertex
    uu[:,1,39]   = np.array([       z,       z,       z,    -eps,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 16 | tLC 40 | Single Vertex
    uu[:,0,40]   = np.array([       z,       z,       z,       z,     eps,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 17 | tLC 41 | Single Vertex
    uu[:,1,41]   = np.array([       z,       z,       z,       z,       z,     eps,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 18 | tLC 42 | Single Vertex
    uu[:,0,42]   = np.array([       z,       z,       z,       z,       z,       z,    -eps,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 19 | tLC 43 | Single Vertex
    uu[:,1,43]   = np.array([       z,       z,       z,       z,       z,       z,       z,    -eps,            z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 20 | tLC 44 | Single Vertex

    uu[:,1,44]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,     eps,    -eps,    -eps ]) # HEXE20 iLC 21 | tLC 45 | Compress Center
    uu[:,2,45]   = np.array([       z,       z,       z,  -3*eps,       z,       z,       z,       z,            z,     eps,       z,     eps,       z,    -eps,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC 22 | tLC 46 | Compress Center
    uu[:,0,46]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,    -eps,    -eps,     eps ]) # HEXE20 iLC 23 | tLC 47 | Compress Center
    uu[:,2,47]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,          eps,       z,     eps,       z,    -eps,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 24 | tLC 48 | Compress Center
    uu[:,0,48]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,    -eps,       z,     eps,       z,    -eps,       z,     eps,       z,       z,       z,       z ]) # HEXE20 iLC 25 | tLC 49 | Compress Center
    uu[:,1,49]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,          eps,       z,    -eps,       z,     eps,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 26 | tLC 50 | Compress Center

    uu[:,0,50]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,          eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 27 | tLC 51 | Move Center
    uu[:,1,51]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,     eps,       z,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 28 | tLC 52 | Move Center
    uu[:,0,52]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 29 | tLC 53 | Move Center
    uu[:,1,53]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,    -eps,       z,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 30 | tLC 54 | Move Center
    uu[:,0,54]   = np.array([       z,       z,       z,   2*eps,       z,       z,       z,       z,            z,       z,       z,       z,     eps,       z,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 31 | tLC 55 | Move Center
    uu[:,1,55]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,     eps,       z,       z,       z,       z,       z,       z ]) # HEXE20 iLC 32 | tLC 56 | Move Center
    uu[:,0,56]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z,       z ]) # HEXE20 iLC 33 | tLC 57 | Move Center
    uu[:,0,57]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z,       z,       z ]) # HEXE20 iLC 34 | tLC 58 | Move Center
    uu[:,2,58]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,     eps,       z,       z,       z ]) # HEXE20 iLC 35 | tLC 59 | Move Center
    uu[:,2,59]   = np.array([       z,       z,       z,       z,       z,       z,       z,       z,            z,       z,       z,       z,       z,       z,       z,       z,       z,    -eps,       z,       z ]) # HEXE20 iLC 36 | tLC 60 | Move Center
