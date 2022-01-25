#include "stdlib.h"
#include <stdio.h>
#include "metis.h"

/*************************************************************************
**************************************************************************/
void F_METIS_PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, 
          idx_t *vwgt, idx_t *vsize, idx_t *ncommon, idx_t *nparts, 
          real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, 
          idx_t *npart) 
{
 
  int rstatus = METIS_OK; 

  vwgt    = NULL;
  vsize   = NULL;
  tpwgts  = NULL;

  printf("## %d ##\n",METIS_NOPTIONS);

  rstatus = METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0; 

  printf("## %d ## %d ## %d ## %d ## %d ##\n",*ne,*nn, *nparts, *ncommon, sizeof(idx_t)); 

// Solid Cube
// ## 40 ##
// ## 8 ## 64 ## 2 ## 4 ## 4 ##


// FH01
// ## 40 ##
// ## 13997521 ## 111980168 ## 2 ## 4 ## 8 ##


  rstatus = METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize,
  			       ncommon, nparts, tpwgts, options,
  			       objval, epart, npart);

  if ( rstatus == METIS_OK ) {
    printf("METIS OK\n");  
  }
  else {
    printf("METIS NOT OK: %d \n", rstatus);
  }
  
}

/*************************************************************************
**************************************************************************/
void F_METIS_PartMeshNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, 
          idx_t *vwgt, idx_t *vsize, idx_t *nparts, 
          real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, 
          idx_t *npart) 
{

  int rstatus = METIS_OK;

  vwgt    = NULL;
  vsize   = NULL;
  tpwgts  = NULL;

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 1;

  printf("## %d ## %d ## %d ##\n",*ne,*nn, *nparts);

  printf("## %d ## %d ## %d ##\n",eptr[0],eptr[*ne], *nparts);

  rstatus = METIS_PartMeshNodal(ne, nn, eptr, eind, vwgt, vsize,
  			       nparts, tpwgts, options,
  			       objval, epart, npart);

  if ( rstatus == METIS_OK ) {
    printf("METIS OK\n"); 
  }
  else {
    printf("METIS NOT OK\n"); 
  }
  
}
