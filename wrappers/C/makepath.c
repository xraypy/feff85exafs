#include <stdio.h>
#include <stdlib.h>

#include "feffpath.h"

int main()
{
  int i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret = create_path(path);

  /* strcpy(path->phpad, "../fortran/phase_orig.pad"); */
  strcpy(path->phpad, "phase.pad");
  path->nnnn    = 1;
  path->xdi     = 1;
  path->verbose = 1;

  /* first path in copper */
  path->index   = 1;
  path->degen   = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 1); /* sole scatterer in an SS path */
  ret = make_path(path);

  /* for (i=0; i<nphx; i++) { */
  /*   printf(" %d  %d\n", i, path->iz[i]); */
  /* }; */

  clear_path(path);

  /* fourth path in copper */
  path->index   = 4;
  path->degen   = 48.0;
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
  /*   	   path->k[i], path->real_phc[i], path->mag_feff[i], path->pha_feff[i], path->red_fact[i], path->lam[i], path->rep[i]); */
  /* } */

  /* printf("%s = %31s \n", "version", path->version); */
  /* printf("%s = %9s  \n", "exch",    path->exch); */
  /* printf("%s = %.5f \n", "edge",    path->edge); */
  /* printf("%s = %.5f \n", "gam_ch",  path->gam_ch); */
  /* printf("%s = %.5f \n", "kf",      path->kf); */
  /* printf("%s = %.5f \n", "mu",      path->mu); */
  /* printf("%s = %.5f \n", "rnorman", path->rnorman); */
  /* printf("%s = %.5f \n", "rs_int",  path->rs_int); */
  /* printf("%s = %.5f \n", "vint",    path->vint); */

  cleanup(path);
  return 0;
}
