#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "feffpath.h"

int main()
{

  int i, j;
  FEFFPATH *path;

  /* scattering and path geometry */
  int index = 9999;
  int iorder = 2;
  int nleg;
  double deg;
  int ipot[legtot+1];
  double rat[legtot+2][3];
  double ri[legtot], beta[legtot+1], eta[legtot+2];

  /* onepath.f output */
  int nnnn, json;

  /* feffNNNN.dat columns */
  int ne;
  double kgrid[nex], caps[nex], amff[nex], phff[nex], redfac[nex], lambda[nex], rep[nex];

  /* polarization and ellipticity */
  int ipol;
  double elpty;
  double evec[3]  = {0,0,0};
  double xivec[3] = {0,0,0};

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

  /* temporarily hardwire in this information */
  index  = 4;
  nleg   = 3;
  deg    = 48.0;
  nnnn   = 1;
  json   = 0;
  ipol   = 0;
  elpty  = 0;

  /* initialize scattering geometry */
  for (i = 0; i < legtot+1; i++) {
    ipot[i] = 0;
  }
  for (i = 0; i < legtot+2; i++) {
    for (j = 0; j < 3; j++) {
      rat[i][j] = 0;
    }
  };
 

  /* --------------------------------------------------- */
  /* allocate array sizes consistent with Feff           */
  path         = calloc(1, sizeof *path);
  assert (path != NULL);

  path->rat  = calloc(legtot+2, sizeof(double *));
  for (i = 0; i < legtot+2; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->ipot   = calloc(legtot+1,   sizeof(int));
  path->ri     = calloc(legtot,   sizeof(double));
  path->beta   = calloc(legtot+1, sizeof(double));
  path->eta    = calloc(legtot+2, sizeof(double));

  path->nnnn   = (nnnn ? true : false);
  path->json   = (json ? true : false);

  path->evec   = calloc(3, sizeof(double));
  path->xivec  = calloc(3, sizeof(double));

  path->kgrid  = calloc(nex, sizeof(double));
  path->caps   = calloc(nex, sizeof(double));
  path->amff   = calloc(nex, sizeof(double));
  path->phff   = calloc(nex, sizeof(double));
  path->redfac = calloc(nex, sizeof(double));
  path->lambda = calloc(nex, sizeof(double));
  path->rep    = calloc(nex, sizeof(double));
  /* --------------------------------------------------- */


  onepath_(&index, &nleg, &deg, &iorder, &ipot, &rat,
	   &ipol, &evec, &elpty, &xivec,
	   &nnnn, &json, &ri, &beta, &eta,
	   &ne, &kgrid, &caps, &amff, &phff, &redfac, &lambda, &rep);


  /* --------------------------------------------------- */
  /* transfer everything into the struct                 */

  /* scattering geometry */
  path->index  = index;
  path->nleg   = nleg;
  path->deg    = deg;
  path->iorder = iorder;
  for (i = 0; i < legtot+1; i++) {
    path->ipot[i] = ipot[i];
  }
  for (i = 0; i < legtot+2; i++) {
    for (j = 0; j < 3; j++) {
      path->rat[i][j] = rat[i][j];
    }
  };

  /* polarization and ellipticity */
  path->ipol   = ipol;
  path->elpty  = elpty;
  for (i = 0; i < 3; i++) {
    path->evec[i]  = evec[0];
    path->xivec[i] = xivec[0];
  }

  /* path geometry */
  for (i = 0; i < legtot; i++) {
    path->ri[i]   = ri[i]   * bohr;
    path->beta[i] = beta[i] * 180 / pi;
    path->eta[i]  = eta[i]  * 180 / pi;
  }

  /* compute Reff for this path */
  path->reff = 0;
  for (i = 0; i < path->nleg; i++) {
    path->reff = path->reff + path->ri[i];
  }
  path->reff = path->reff / 2;

  /* array of F_eff */
  path->ne = ne;
  for (i = 0; i < path->ne; i++) {
    path->kgrid[i]  = kgrid[i];
    path->caps[i]   = caps[i];
    path->amff[i]   = amff[i];
    path->phff[i]   = phff[i];
    path->redfac[i] = redfac[i];
    path->lambda[i] = lambda[i];
    path->rep[i]    = rep[i];
    /* printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n", */
    /* 	   path->kgrid[i], path->caps[i], path->amff[i], path->phff[i], path->redfac[i], path->lambda[i], path->rep[i]); */
  }
  /* --------------------------------------------------- */

  for (i = 0; i <= path->nleg; i++) {
    printf("%8.5f  %8.5f  %8.5f  %d\n", path->rat[i][0]*bohr, path->rat[i][1]*bohr, path->rat[i][2]*bohr, path->ipot[i]);
  }

  return 0;
}


/* methods */
/* */
/*   1.  initialization, up to line 110 */
/*   2.  add atom: adds a position/ipot to rat and ipot, returns new nleg */
/*   3.  compute: call onepath_, filling the output arrays, lines 112-169 */
/* */
