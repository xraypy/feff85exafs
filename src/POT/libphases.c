/* ********************************************************************** */
/* LICENSE AND COPYRIGHT                                                  */
/*                                                                        */
/* To the extent possible, the authors have waived all rights granted by  */
/* copyright law and related laws for the code and documentation that     */
/* make up the C Interface to the potph library.  While information       */
/* about Authorship may be retained in some files for historical reasons, */
/* this work is hereby placed in the Public Domain.  This work is         */
/* published from: United States.                                         */
/*                                                                        */
/* Note that the potph library itself is NOT public domain, nor is the    */
/* Fortran source code for Feff that it relies upon.                      */
/*                                                                        */
/* Author: Bruce Ravel (bravel AT bnl DOT gov).                           */
/* Created: 22 May, 2015                                                  */
/* ********************************************************************** */


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "libphases.h"

_EXPORT(long) create_phases(FEFFPHASES *phases) {
  /* Instantiate and initialize a FEFFPHASES struct */
  /* Return an error code -- currently hardwired to return 0. */
  long i;

  char titl[81] = {'\0'};
  char ptlb[7]  = {'\0'};

  phases->ntitle = 0;
  phases->nat    = 0;
  phases->nph    = 0;
  phases->ihole  = 1;
  phases->rscf   = 0.0;
  phases->lscf   = 0;
  phases->nscmt  = 0;
  phases->ca     = 0.0;
  phases->nmix   = 0;
  phases->ecv    = 0.0;
  phases->icoul  = 0;
  phases->ipol   = 0;
  phases->elpty  = 0.0;
  phases->ispin  = 0;
  phases->angks  = 0.0;
  phases->gamach = 0.0;
  phases->iafolp = 0;
  phases->rgrd   = 0.0;
  phases->iunf   = 0;
  phases->inters = 0;
  phases->totvol = 0.0;
  phases->jumprm = 0;
  phases->nohole = 0;

  /* allocate string arrays for titles and potlbl */


  phases->rat  = calloc(natx, sizeof(double *));
  for (i = 0; i < natx; i++) {
    phases->rat[i] = calloc(3, sizeof(double));
  }

  phases->iz     = calloc(nphx+1, sizeof(long));
  phases->lmaxsc = calloc(nphx+1, sizeof(long));
  phases->lmaxph = calloc(nphx+1, sizeof(long));
  phases->xnatph = calloc(nphx+1, sizeof(double));
  phases->spinph = calloc(nphx+1, sizeof(double));

  phases->evec   = calloc(3, sizeof(double));
  phases->xivec  = calloc(3, sizeof(double));
  phases->spvec  = calloc(3, sizeof(double));

  phases->folp   = calloc(nphx+1, sizeof(double));
  phases->xion   = calloc(nphx+1, sizeof(double));


  phases->ptz  = calloc(3, sizeof(double *));
  for (i = 0; i < 3; i++) {
    phases->ptz[i] = calloc(3, sizeof(double));
  }

  phases->titles = malloc(sizeof(char *) *  nheadx * 7);
  for (i = 0; i < nheadx; i++) {
    strcpy(phases->titles[i], titl);
  }
  phases->potlbl = malloc(sizeof(char *) * (nphx+1) * 81);
  for (i = 0; i < nphx+1; i++) {
    strcpy(phases->potlbl[i], ptlb);
  }

  return 0;
}


_EXPORT(void) clear_phases(FEFFPHASES *phases) {
  /* Reinitialize a FEFFPHASES struct, returning everything to default */
  long i,j;

  phases->ntitle = 0;
  phases->nat    = 0;
  phases->nph    = 0;
  phases->ihole  = 1;
  phases->rscf   = 0.0;
  phases->lscf   = 0;
  phases->nscmt  = 0;
  phases->ca     = 0.0;
  phases->nmix   = 0;
  phases->ecv    = 0.0;
  phases->icoul  = 0;
  phases->ipol   = 0;
  phases->elpty  = 0.0;
  phases->ispin  = 0;
  phases->angks  = 0.0;
  phases->gamach = 0.0;
  phases->iafolp = 0;
  phases->rgrd   = 0.0;
  phases->iunf   = 0;
  phases->inters = 0;
  phases->totvol = 0.0;
  phases->jumprm = 0;
  phases->nohole = 0;


  for (i = 0; i <= nheadx; i++) {
    strcpy(phases->titles[i], "");
  }

  for (i = 0; i <= nphx; i++) {
    phases->iz[i]     = 0;
    phases->lmaxsc[i] = 0;
    phases->lmaxph[i] = 0;
    phases->xnatph[i] = 0.0;
    phases->spinph[i] = 0.0;
    phases->folp[i]   = 0.0;
    phases->xion[i]   = 0.0;
    strcpy(phases->potlbl[i], "");
  }

  for (i = 0; i < natx; i++) {
    for (j = 0; j < 3; j++) {
      phases->rat[i][j] = 0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      phases->ptz[i][j] = 0;
    }
  }

  for (i = 0; i < 3; i++) {
    phases->evec[i]  = 0.0;
    phases->xivec[i] = 0.0;
    phases->spvec[i] = 0.0;
  }

}


_EXPORT(long) make_phases(FEFFPHASES *phases) {
  /* Instantiate and initialize a FEFFPHASES struct */
  /* Return an error code -- currently hardwired to return 0. */
  long i, j;

  long ntitle, nat, nph, ihole, lscf, nscmt, nmix, icoul, ipol, ispin, iafolp, iunf, inters, jumprm, nohole;
  double rscf, ca, ecv, elpty, angks, gamach, rgrd, totvol;
  long iz[nphx+1], lmaxsc[nphx+1], lmaxph[nphx+1], iphat[natx];
  double xnatph[nphx+1], spinph[nphx+1], folp[nphx+1], xion[nphx+1];
  double rat[natx][3], ptz[3][3];
  double evec[3], xivec[3], spvec[3];
  char potlbl[nphx+1][7];
  char titles[nheadx][81];


  for (i = 0; i < ntitle; i++) {
    strcpy(titles[i], phases->titles[i]);
  };

  for (i = 0; i < nph; i++) {
    iz[i]     = phases->iz[i];
    lmaxsc[i] = phases->lmaxsc[i];
    lmaxph[i] = phases->lmaxph[i];
    xnatph[i] = phases->xnatph[i];
    spinph[i] = phases->spinph[i];
    folp[i]   = phases->folp[i];
    xion[i]   = phases->xion[i];
    strcpy(potlbl[i], phases->potlbl[i]);
  };

  for (i = 0; i < nat; i++) {
    for (j = 0; j < 3; j++) {
      rat[i][j] = phases->rat[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      ptz[i][j] = phases->ptz[i][j];
    }
  }


  for (i = 0; i < 3; i++) {
    evec[i]  = phases->evec[i];
    xivec[i] = phases->xivec[i];
    spvec[i] = phases->spvec[i];
  }


  libpotph_(&ntitle, &titles, &nat, &rat, &iphat,
	    &nph, &iz, &potlbl, &lmaxsc, &lmaxph, &xnatph, &spinph,
	    &ihole, &rscf, &lscf, &nscmt, &ca, &nmix, &ecv, &icoul,
	    &ipol, &evec, &elpty, &xivec, &ispin, &spvec, &angks,
	    &ptz, &gamach, 
	    &iafolp, &folp, &xion, &rgrd, &iunf, &inters, &totvol, &jumprm, &nohole );

  return 0;
}


_EXPORT(void) cleanup(FEFFPHASES *phases) {

  long i;

  for (i = 0; i < natx; i++) {
    free(phases->rat[i]);
  }
  free(phases->rat);

  free(phases->iz);
  free(phases->lmaxsc);
  free(phases->lmaxph);
  free(phases->xnatph);
  free(phases->spinph);

  free(phases->evec);
  free(phases->xivec);
  free(phases->spvec);

  free(phases->folp);
  free(phases->xion);


  for (i = 0; i < 3; i++) {
    free(phases->ptz[i]);
  }
  free(phases->ptz);

  free(phases->titles);
  free(phases->potlbl);

};
