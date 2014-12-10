/* ********************************************************************** */
/* LICENSE AND COPYRIGHT                                                  */
/*                                                                        */
/* To the extent possible, the authors have waived all rights granted by  */
/* copyright law and related laws for the code and documentation that     */
/* make up the C Interface to the feffpath library.  While information    */
/* about Authorship may be retained in some files for historical reasons, */
/* this work is hereby placed in the Public Domain.  This work is         */
/* published from: United States.                                         */
/*                                                                        */
/* Note that the onepath library itself is NOT public domain, nor is the  */
/* Fortran source code for Feff that it relies upon.                      */
/*                                                                        */
/* Author: Bruce Ravel (bravel AT bnl DOT gov).                           */
/* Last update: 5 December, 2014                                          */
/* ********************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "feffpath.h"

_EXPORT(long) create_path(FEFFPATH *path) {
  /* Instantiate and initialize a FEFFPATH struct */
  /* Return an error code -- currently hardwired to return 0. */
  long i;
  char phbin[256] = {'\0'};
  sprintf(phbin, "%-256s", " ");

  path->index     = 9999;
  path->iorder    = 2;
  path->nleg      = 0;
  path->deg       = 1.0;
  path->nnnn      = false;
  path->json      = false;
  path->verbose   = false;
  path->ipol      = false;
  path->elpty     = 0.0;
  path->ne        = 0;
  path->errorcode = 0;

  /* --------------------------------------------------- */
  /* allocate array sizes consistent with Feff           */
  /* path         = calloc(1, sizeof *path); */
  assert (path != NULL);

  path->rat  = calloc(legtot+2, sizeof(double *));
  for (i = 0; i < legtot+2; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->ipot     = calloc(legtot+1, sizeof(long));
  path->ri       = calloc(legtot,   sizeof(double));
  path->beta     = calloc(legtot+1, sizeof(double));
  path->eta      = calloc(legtot+2, sizeof(double));

  path->evec     = calloc(3, sizeof(double));
  path->xivec    = calloc(3, sizeof(double));

  path->k        = calloc(nex, sizeof(double));
  path->real_phc = calloc(nex, sizeof(double));
  path->mag_feff = calloc(nex, sizeof(double));
  path->pha_feff = calloc(nex, sizeof(double));
  path->red_fact = calloc(nex, sizeof(double));
  path->lam      = calloc(nex, sizeof(double));
  path->realp    = calloc(nex, sizeof(double));
  /* --------------------------------------------------- */

  COPY_STRING(path->errormessage, "");
  COPY_STRING(path->phbin, phbin);

  return 0;
}

_EXPORT(void) clear_path(FEFFPATH *path) {
  /* Reinitialize a FEFFPATH struct, returning everything to default */
  /* except phbin, verbose, nnnn, and json. */
  long i,j;
  /* char phbin[256] = {'\0'}; */
  /* sprintf(phbin, "%-256s", " "); */

  path->index     = 9999;
  path->iorder    = 2;
  path->nleg      = 0;
  path->deg       = 0.0;
  /* path->nnnn       = 0; */
  /* path->json       = 0; */
  /* path->verbose    = 0; */
  path->ipol      = 0;
  path->elpty     = 0.0;
  path->reff      = 0.0;
  path->ne        = 0;
  path->errorcode = 0;

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
    path->realp[i]    = 0;
  }

  COPY_STRING(path->errormessage, "");
  /* COPY_STRING(path->phbin, "phase.bin"); */
}


_EXPORT(long) make_path(FEFFPATH *path) {
  /* Compute the path using the content of the FEFFPATH struct.*/
  /* Fill the struct with geometry and F_eff information. */
  /* Return an error code. */
  long i, j;
  long error = 0;

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
  double k[nex], real_phc[nex], mag_feff[nex], pha_feff[nex], red_fact[nex], lam[nex], realp[nex];

  /* polarization and ellipticity */
  long ipol;
  double elpty;
  double evec[3];
  double xivec[3];

  char phbin[256] = {'\0'};

  /* printf("entering make_path\n"); */
  /* fflush(stdout); */
  sprintf(phbin, "%-256s", path->phbin);
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

  if ( (rat[1][0] == 0) && (rat[1][1] == 0) && (rat[1][2] == 0) ) {
    error = error + ERR_FIRSTISABS;
  };
  if ( (rat[nleg-1][0] == 0) && (rat[nleg-1][1] == 0) && (rat[nleg-1][2] == 0) ) {
    error = error + ERR_NLEGISABS;
  };
  if (deg < 0) { /* degeneracy must be non-negative */
    error = error + ERR_DEGNEG;
  };
  if ( (index < 0) || (index > 9999) ) { /* 0 <= index <= 9999 */
    error = error + ERR_BADINDEX;
  };
  if ( (elpty < 0) || (elpty > 1) ) { /* 0 <= elpty <= 1 */
    error = error + ERR_BADELPTY;
  };
  if ( (iorder < 0) || (iorder > 10) ) { /* 0 <= iorder <= 10 */
    error = error + ERR_BADIORDER;
  };
  /* if( fopen(phbin, "r") == NULL ) { */
  /*   error = error + ERR_NOPHBIN; */
  /* } */
  path->errorcode = error;
  make_path_errorstring(path);
  /* printf("%ld\n", path->errorcode); */
  /* fflush(stdout); */
  if (error > 0) {
    return error;
  };

  /* printf(">%s<\n", phbin); */
  /* fflush(stdout); */
  onepath_(phbin, &index, &nleg, &deg, &iorder, &ipot, &rat,
	   &ipol, &evec, &elpty, &xivec,
	   &nnnn, &json, &verbose, &ri, &beta, &eta,
	   &ne, &k, &real_phc, &mag_feff, &pha_feff, &red_fact, &lam, &realp);
  /* printf("after onepath_\n"); */
  /* fflush(stdout); */

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
  /* printf("after geometry\n"); */
  /* fflush(stdout); */

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
    path->realp[i]    = realp[i];
    /* printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n", */
    /* 	   path->kgrid[i], path->caps[i], path->amff[i], path->phff[i], path->redfac[i], path->lambda[i], path->realp[i]); */
  }
  /* --------------------------------------------------- */
  /* printf("after arrays\n"); */
  /* fflush(stdout); */

  return error;
}

_EXPORT(long) add_scatterer(FEFFPATH *path, double x, double y, double z, long ip) {
  /* Add a scattering atom to the path. Fill the nleg member. */
  /* Return an error code. */
  long error;
  double length;
  long nleg = path->nleg;
  if (nleg == 0) {nleg = 1;}
  nleg = nleg + 1;
  
  path->rat[nleg-1][0] = x;
  path->rat[nleg-1][1] = y;
  path->rat[nleg-1][2] = z;
  path->ipot[nleg-1]   = ip;
  path->nleg           = nleg;
  
  error = 0;
  if (ip < 0) {
    error = error + ERR_NEGIPOT;
  };
  if (ip > 7) {
    error = error + ERR_BIGIPOT;
  };
  length = leglength(path);
  if (length < 0.5) {
    error = error + ERR_TOOCLOSE;
  };
  if (nleg > legtot) {
    error = error + ERR_TOOMANYLEGS;
  };

  path->errorcode = error;
  make_scatterer_errorstring(path);
  return error;
}

_EXPORT(void) cleanup(FEFFPATH *path) {
  /* does one need to explicitly free each part of the struct? */
  /* Google for "c free struct memory" */
  /* http://stackoverflow.com/questions/5324355/how-do-i-release-a-struct-from-memory-and-arrays-within-them */
  /* http://stackoverflow.com/questions/6720781/c-free-and-struct */
  free(path);
}

leglength(FEFFPATH *path) {
  double x1, y1, z1, x2, y2, z2;
  long nleg = path->nleg;
  x1 = path->rat[nleg-2][0];
  y1 = path->rat[nleg-2][1];
  z1 = path->rat[nleg-2][2];
  x2 = path->rat[nleg-1][0];
  y2 = path->rat[nleg-1][1];
  z2 = path->rat[nleg-1][2];
  return sqrt( pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2) );
}


/* error string interpretation */
make_scatterer_errorstring(FEFFPATH *path) {
  double x, y, z;
  long ip;
  char message[500] = {'\0'};
  char pos[100];
  long errcode = path->errorcode;
  long nleg = path->nleg;
  x  = path->rat[nleg-1][0];
  y  = path->rat[nleg-1][1];
  z  = path->rat[nleg-1][2];
  ip = path->ipot[nleg-1];

  if (errcode == 0) { return; }
  sprintf(pos, "Error in add_scatterer at atom (%.5f, %.5f, %.5f, %ld):\n", x, y, z, ip);
  strcat(message, pos);
  if (errcode & ERR_NEGIPOT) {
    strcat(message, "\t(code 1) ipot argument to add_scatterer is less than 0\n");
  };
  if (errcode & ERR_BIGIPOT) {
    strcat(message, "\t(code 2) ipot argument to add_scatterer is greater than 7\n");
  } ;
  if (errcode & ERR_TOOCLOSE) {
    strcat(message, "\t(code 4) coordinates are for an atom too close to the previous atom in the path\n");
  };
  if (errcode & ERR_TOOMANYLEGS) {
    strcat(message, "\t(code 8) nlegs greater than legtot\n");
  };
  COPY_STRING(path->errormessage, message);
}

make_path_errorstring(FEFFPATH *path) {
  double deg, elpty;
  long index, iorder;
  char message[240] = {'\0'};
  char str[100];
  long errcode = path->errorcode;
  long nleg = path->nleg;
  char phbin[256] = {'\0'};
  deg    = path->deg;
  index  = path->index;
  iorder = path->iorder;
  elpty  = path->elpty;
  sprintf(phbin, "%-256s", path->phbin);

  if (errcode == 0) { return; }
  strcat(message, "Error in make_path\n");
  if (errcode & ERR_FIRSTISABS) {
    strcat(message, "\t(code 1) the first atom specified is the absorber\n");
  };
  if (errcode & ERR_NLEGISABS) {
    strcat(message, "\t(code 2) the last atom specified is the absorber\n");
  };
  if (errcode & ERR_DEGNEG) {
    sprintf(str, "\t(code 4) path degeneracy (%.2f) is negative\n", deg);
    strcat(message, str);
  };
  if (errcode & ERR_BADINDEX) {
    sprintf(str, "\t(code 8) path index (%ld) not between 0 and 9999\n", index);
    strcat(message, str);
  };
  if (errcode & ERR_BADELPTY) {
    sprintf(str, "\t(code 16) ellipticity (%.2f) not between 0 and 1\n", elpty);
    strcat(message, str);
  };
  if (errcode & ERR_BADIORDER) {
    sprintf(str, "\t(code 32) iorder (%ld) not between 0 and 10\n", iorder);
    strcat(message, str);
  };
  if (errcode & ERR_NOPHBIN) {
    sprintf(str, "\t(code 64) phase.bin file (%s) does not exist\n", phbin);
    strcat(message, str);
  };
  if (errcode & ERR_PHBINNOREAD) {
    sprintf(str, "\t(code 128) phase.bin file (%s) cannot be read\n", phbin);
    strcat(message, str);
  };
  COPY_STRING(path->errormessage, message);
}

