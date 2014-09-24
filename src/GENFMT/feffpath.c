#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "feffpath.h"

_EXPORT(int) initialize_path(FEFFPATH *path) {

  int i;

  path->index  = 9999;
  path->iorder = 2;
  path->nleg   = 0;
  path->deg    = 0.0;
  path->nnnn   = 0;
  path->json   = 0;
  path->ipol   = 0;
  path->elpty  = 0.0;
  path->ne     = 0;

  /* --------------------------------------------------- */
  /* allocate array sizes consistent with Feff           */
  /* path         = calloc(1, sizeof *path); */
  /* assert (path != NULL); */

  path->rat  = calloc(legtot+2, sizeof(double *));
  for (i = 0; i < legtot+2; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->ipot   = calloc(legtot+1, sizeof(int));
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

_EXPORT(int) make_path(FEFFPATH *path) {

  int i, j;

  /* scattering and path geometry */
  int index;
  int iorder;
  int nleg;
  double deg;
  int ipot[legtot+1];
  double rat[legtot+2][3];
  double ri[legtot], beta[legtot+1], eta[legtot+2];

  /* onepath.f output */
  int nnnn, json;

  /* feffNNNN.dat columns */
  int ne;
  double k[nex], real_phc[nex], mag_feff[nex], pha_feff[nex], red_fact[nex], lam[nex], rep[nex];

  /* polarization and ellipticity */
  int ipol;
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
  nnnn = path->nnnn;
  json = path->json;
  

  void onepath_(int *,                   /* path index */
		int *,                   /* nlegs */
		double *,                /* degeneracy */
		int *,                   /* iorder */
		/* scattering geometry */
		int (*)[legtot+1],       /* list of unique potentials */
		double (*)[legtot+2][3], /* list of cartesian coordinates */
		/* polarization and ellipticity */
		int *,                   /* flag to compute polarization */
		double (*)[3],           /* polarization vector */
		double *,                /* ellipticity */
		double (*)[3],           /* direction of travel */
		/* ouyput flags */
		int *,                   /* integer flag for writing feffNNNN.dat */
		int *,                   /* integer flag for writing feffNNNN.json */
		/* path geometry */
		double (*)[legtot],      /* Ri   */
		double (*)[legtot+1],    /* beta */
		double (*)[legtot+2],    /* eta  */
		int *,                   /* number of points in kgrid */
		/* seven columns of feffNNNN.dat file */
		double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex]);


  onepath_(&index, &nleg, &deg, &iorder, &ipot, &rat,
	   &ipol, &evec, &elpty, &xivec,
	   &nnnn, &json, &ri, &beta, &eta,
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


/* methods */
/* */
/*   1.  initialization, up to line 110 */
/*   2.  add atom: adds a position/ipot to rat and ipot, returns new nleg */
/*   3.  compute: call onepath_, filling the output arrays, lines 112-169 */
/* */
