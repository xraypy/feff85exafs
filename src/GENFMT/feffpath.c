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

  int i;
  FEFFPATH *path;

  /* path geometry */
  int index, nleg;
  double deg;
  double ri[legtot], beta[legtot+1], eta[legtot+2];
  /* onepath.f output */
  int nnnn, json;
  /* feffNNNN.dat columns */
  int ne;
  double kgrid[nex], caps[nex], amff[nex], phff[nex], redfac[nex], lambda[nex], rep[nex];


  void onepath_(int *,      /* path index */
		int *,      /* nlegs */
		double *,   /* degeneracy */
		int *,      /* integer flag for writing feffNNNN.dat */
		int *,      /* integer flag for writing feffNNNN.json */
		/* path geometry */
		double (*)[legtot],   /* Ri   */
		double (*)[legtot+1], /* beta */
		double (*)[legtot+2], /* eta  */
		int *,      /* number of points in kgrid */
		/* seven columns of feffNNNN.dat file */
		double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex]);

  index = 4;
  nleg  = 3;
  deg   = 48.0;
  nnnn  = 1;
  json  = 0;

  /* allocate array sizes consistent with Feff */
  path         = calloc(1, sizeof *path);
  assert (path != NULL);

  path->rat  = calloc(legtot, sizeof(double *));
  for (i = 0; i < legtot; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->ipot   = calloc(legtot,   sizeof(int));
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

  onepath_(&index, &nleg, &deg, &nnnn, &json, &ri, &beta, &eta,
	   &ne, &kgrid, &caps, &amff, &phff, &redfac, &lambda, &rep);

  path->ne    = ne;
  path->index = index;
  path->nleg  = nleg;
  path->deg   = deg;

  for (i = 0; i < legtot; i++) {
    path->ri[i]   = ri[i] * bohr;
    path->beta[i] = beta[i] * 180 / pi;
    path->eta[i]  = eta[i] * 180 / pi;
  }

  path->reff = 0;
  for (i = 0; i < path->nleg; i++) {
    path->reff = path->reff + path->ri[i];
  }
  path->reff = bohr * path->reff / 2;

  /* transfer array values into the struct */
  for (i = 0; i < path->ne; i++) {
    path->kgrid[i]  = kgrid[i];
    path->caps[i]   = caps[i];
    path->amff[i]   = amff[i];
    path->phff[i]   = phff[i];
    path->redfac[i] = redfac[i];
    path->lambda[i] = lambda[i];
    path->rep[i]    = rep[i];

    printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n",
	   path->kgrid[i], path->caps[i], path->amff[i], path->phff[i], path->redfac[i], path->lambda[i], path->rep[i]);
  }
  printf(" %.5f %.5f %.5f\n", path->ri[0], path->ri[1], path->ri[2]);
  printf(" %.5f %.5f %.5f\n", path->beta[0], path->beta[1], path->beta[2]);
  printf(" %.5f %.5f %.5f\n", path->eta[0], path->eta[1], path->eta[2]);

  return 0;
}
