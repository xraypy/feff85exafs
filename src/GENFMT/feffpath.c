#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "feffpath.h"

_EXPORT(long) create_path(FEFFPATH *path) {
  /* Instantiate and initialize a FEFFPATH struct */
  /* Return an error code -- currently returns 0 in all cases. */
  long i;

  path->index   = 9999;
  path->iorder  = 2;
  path->nleg    = 0;
  path->deg     = 0.0;
  path->nnnn    = false;
  path->json    = false;
  path->verbose = false;
  path->ipol    = false;
  path->elpty   = 0.0;
  path->ne      = 0;

  /* --------------------------------------------------- */
  /* allocate array sizes consistent with Feff           */
  /* path         = calloc(1, sizeof *path); */
  assert (path != NULL);

  path->rat  = calloc(legtot+2, sizeof(double *));
  for (i = 0; i < legtot+2; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->ipot   = calloc(legtot+1, sizeof(long));
  path->ri     = calloc(legtot,   sizeof(double));
  path->beta   = calloc(legtot+1, sizeof(double));
  path->eta    = calloc(legtot+2, sizeof(double));

  path->evec   = calloc(3, sizeof(double));
  path->xivec  = calloc(3, sizeof(double));

  path->k        = calloc(nex, sizeof(double));
  path->real_phc = calloc(nex, sizeof(double));
  path->mag_feff = calloc(nex, sizeof(double));
  path->pha_feff = calloc(nex, sizeof(double));
  path->red_fact = calloc(nex, sizeof(double));
  path->lam      = calloc(nex, sizeof(double));
  path->rep      = calloc(nex, sizeof(double));
  /* --------------------------------------------------- */

  return 0;
}

_EXPORT(void) clear_path(FEFFPATH *path) {
  /* Reinitialize a FEFFPATH struct, returning everything to default except the three boolean members. */
  long i,j;

  path->index  = 9999;
  path->iorder = 2;
  path->nleg   = 0;
  path->deg    = 0.0;
  /* path->nnnn    = 0; */
  /* path->json    = 0; */
  /* path->verbose = 0; */
  path->ipol   = 0;
  path->elpty  = 0.0;
  path->reff   = 0.0;
  path->ne     = 0;
  for (i = 0; i < 3; i++) {
    path->evec[i]  = 0;
    path->xivec[i] = 0;
  }
  for (i = 0; i < legtot; i++) {
    path->ipot[i] = 0;
    path->ri[i]   = 0;
    path->beta[i] = 0;
    path->eta[i]  = 0;
  }
  path->ipot[legtot+1] = 0;
  path->beta[legtot+1] = 0;
  path->eta[legtot+1]  = 0;
  path->eta[legtot+2]  = 0;
  for (i = 0; i < legtot+2; i++) {
    for (j = 0; j < 3; j++) {
      path->rat[i][j] = 0;
    }
  }
  for (i = 0; i < nex; i++) {
    path->k[i]        = 0;
    path->real_phc[i] = 0;
    path->mag_feff[i] = 0;
    path->pha_feff[i] = 0;
    path->red_fact[i] = 0;
    path->lam[i]      = 0;
    path->rep[i]      = 0;
  }
}


_EXPORT(long) make_path(FEFFPATH *path) {
  /* Compute the path using the content of the FEFFPATH struct.*/
  /* Fill the struct with geometry and F_eff information. */
  /* Return an error code -- currently returns 0 in all cases. */
  long i, j;

  /* scattering and path geometry */
  long index;
  long iorder;
  long nleg;
  double deg;
  long ipot[legtot+1];
  double rat[legtot+2][3];
  double ri[legtot], beta[legtot+1], eta[legtot+2];

  /* onepath.f output */
  long nnnn, json, verbose;

  /* feffNNNN.dat columns */
  long ne;
  double k[nex], real_phc[nex], mag_feff[nex], pha_feff[nex], red_fact[nex], lam[nex], rep[nex];

  /* polarization and ellipticity */
  long ipol;
  double elpty;
  double evec[3];
  double xivec[3];

  iorder = path->iorder;
  index = path->index;
  nleg = path->nleg;
  deg = path->deg;
  for (i = 0; i < legtot+1; i++) {
    ipot[i] = path->ipot[i];
  }
  for (i = 0; i < legtot+2; i++) {
    for (j = 0; j < 3; j++) {
      rat[i][j] = path->rat[i][j]/bohr;
    }
  };
  ipol = path->ipol;
  elpty = path->elpty;
  for (i = 0; i < 3; i++) {
    evec[i]  = path->evec[i];
    xivec[i] = path->xivec[i];
  }
  nnnn    = path->nnnn;
  json    = path->json;
  verbose = path->verbose;
  
  onepath_(&index, &nleg, &deg, &iorder, &ipot, &rat,
	   &ipol, &evec, &elpty, &xivec,
	   &nnnn, &json, &verbose, &ri, &beta, &eta,
	   &ne, &k, &real_phc, &mag_feff, &pha_feff, &red_fact, &lam, &rep);


  /* --------------------------------------------------- */
  /* transfer everything into the struct                 */

  /* path geometry */
  for (i = 0; i < legtot; i++) {
    path->ri[i]   = ri[i]   * bohr;
    path->beta[i] = beta[i] * 180 / pi;
    path->eta[i]  = eta[i]  * 180 / pi;
  }
  path->beta[legtot+1] = beta[legtot+1];
  path->eta[legtot+1]  = eta[legtot+1];
  path->eta[legtot+2]  = eta[legtot+2];

  /* compute Reff for this path */
  path->reff = 0;
  for (i = 0; i < path->nleg; i++) {
    path->reff = path->reff + path->ri[i];
  }
  path->reff = path->reff / 2;

  /* array of F_eff */
  path->ne = ne;
  for (i = 0; i < path->ne; i++) {
    path->k[i]        = k[i];
    path->real_phc[i] = real_phc[i];
    path->mag_feff[i] = mag_feff[i];
    path->pha_feff[i] = pha_feff[i];
    path->red_fact[i] = red_fact[i];
    path->lam[i]      = lam[i];
    path->rep[i]      = rep[i];
    /* printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n", */
    /* 	   path->kgrid[i], path->caps[i], path->amff[i], path->phff[i], path->redfac[i], path->lambda[i], path->rep[i]); */
  }
  /* --------------------------------------------------- */

  return 0;
}

_EXPORT(long) add_scatterer(FEFFPATH *path, double x, double y, double z, long ip) {
  /* Add a scattering atom to the path. */
  /* Return the current value of nleg. */
  long nleg = path->nleg;
  if (nleg == 0) {nleg = 1;}
  nleg = nleg + 1;
  
  /* printf(">>> %.3f  %.3f  %.3f  %li   %li\n", x, y, z, ip, nleg); */

  path->rat[nleg-1][0] = x;
  path->rat[nleg-1][1] = y;
  path->rat[nleg-1][2] = z;
  path->ipot[nleg-1]   = ip;
  path->nleg           = nleg;

  return nleg;
}

_EXPORT(void) cleanup(FEFFPATH *path) {
  free(path);
}
