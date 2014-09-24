#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "feffpath.h"

int main()
{
  int i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret = initialize_path(path);

  /* path 4 from copper metal */
  /*  0.0000    0.0000    0.0000  0  29 Cu       absorbing atom */
  /*  0.0000    0.0000   -3.6100  1  29 Cu */
  /* -1.8050    0.0000   -1.8050  1  29 Cu */

  path->index = 4;
  path->nleg  = 3;
  path->deg   = 48.0;
  path->nnnn  = 1;
  path->json  = 0;

  path->ipot[1] = 1;
  path->ipot[2] = 1;

  path->rat[1][0] = 0;
  path->rat[1][1] = 0;
  path->rat[1][2] = -3.610000;

  path->rat[2][0] = -1.805000;
  path->rat[2][1] = 0;
  path->rat[2][2] = -1.805000;

  /* ret = add_scatterer(path, 0, 0, -3.61, 1); */

  ret = make_path(path);
  

  for (i = 0; i < path->ne; i++) {
    printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n",
    	   path->k[i], path->real_phc[i], path->mag_feff[i], path->pha_feff[i], path->red_fact[i], path->lam[i], path->rep[i]);
  }



}
