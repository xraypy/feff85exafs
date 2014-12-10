#include <stdio.h>
#include <stdlib.h>

#include "feffpath.h"

long main()
{
  long i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret = create_path(path);

  strcpy(path->phbin, "../fortran/phase.bin");
  path->nnnn    = 1;
  path->json    = 0;
  path->verbose = 1;

  /* first path in copper */
  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 1); /* sole scatterer in an SS path */
  ret = make_path(path);

  clear_path(path);

  /* fourth path in copper */
  path->index   = 4;
  path->deg     = 48.0;
  ret = add_scatterer(path,  0,     0, -3.61,  1); /* first atom after absorber */
  ret = add_scatterer(path, -1.805, 0, -1.805, 1); /* last atom before absorber */
  if (ret > 0) {  /* check for errors */
    printf("%s", path->errormessage);
    exit(ret);
  };
  ret = make_path(path);
  if (ret > 0) {  /* check for errors */
    printf("%s", path->errormessage);
    exit(ret);
  };
  
  /* for (i = 0; i < path->ne; i++) { */
  /*   printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n", */
  /*   	   path->k[i], path->real_phc[i], path->mag_feff[i], path->pha_feff[i], path->red_fact[i], path->lam[i], path->realp[i]); */
  /* } */

  cleanup(path);
  return 0;
}
