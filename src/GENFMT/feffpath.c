/* ********************************************************************** */
/* LICENSE AND COPYRIGHT                                                  */
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

_EXPORT(int) create_path(FEFFPATH *path) {
  /* Instantiate and initialize a FEFFPATH struct */
  /* Return an error code -- currently hardwired to return 0. */
  int i;
  char message[500] = {'\0'};
  char phpad[257] = {'\0'};
  char exch[9] = {'\0'};
  char version[31] = {'\0'};

  strcpy(message, " ");
  strcpy(phpad, " ");
  strcpy(exch, " ");
  strcpy(version, " ");

  path->index     = 9999;
  path->iorder    = 2;
  path->nleg      = 0;
  path->degen     = 1.0;
  path->nnnn      = 0;
  path->xdi       = 0;
  path->verbose   = 0;
  path->ipol      = 0;
  path->elpty     = 0.0;
  path->ne        = 0;
  path->errorcode = 0;

  path->edge      = 0.0;
  path->gam_ch    = 0.0;
  path->kf        = 0.0;
  path->mu        = 0.0;
  path->rnorman   = 0.0;
  path->rs_int    = 0.0;
  path->vint      = 0.0;
  path->exch = calloc(9, sizeof(char));
  strcpy(path->exch, exch);
  path->version = calloc(31, sizeof(char));
  strcpy(path->version, version);


  /* --------------------------------------------------- */
  /* allocate array sizes consistent with Feff           */
  /* path         = calloc(1, sizeof *path); */
  assert (path != NULL);

  path->rat  = calloc(legtot+2, sizeof(double *));
  for (i = 0; i < legtot+2; i++) {
    path->rat[i] = calloc(3, sizeof(double));
  }
  path->iz       = calloc(nphx+1,   sizeof(int));
  path->ipot     = calloc(legtot+1, sizeof(int));
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
  path->rep      = calloc(nex, sizeof(double));
  /* --------------------------------------------------- */

  /* COPY_STRING(path->errormessage, message); */
  /* COPY_STRING(path->phpad, phpad); */
  path->errormessage = calloc(500, sizeof(char));
  strcpy(path->errormessage, message);
  path->phpad = calloc(257, sizeof(char));
  strcpy(path->phpad, phpad);

  return 0;
}

_EXPORT(void) clear_path(FEFFPATH *path) {
  /* Reinitialize a FEFFPATH struct, returning everything to default */
  /* except phpad, verbose, nnnn, xdi. */
  int i,j;
  /* char phpad[256] = {'\0'}; */
  /* sprintf(phpad, "%-256s", " "); */

  path->index     = 9999;
  path->iorder    = 2;
  path->nleg      = 0;
  path->degen     = 1.0;
  /* path->nnnn       = 0; */
  /* path->xdi        = 0; */
  /* path->verbose    = 0; */
  path->ipol      = 0;
  path->elpty     = 0.0;
  path->reff      = 0.0;
  path->ne        = 0;
  path->errorcode = 0;

  path->edge      = 0;
  path->gam_ch    = 0;
  path->kf        = 0;
  path->mu        = 0;
  path->rnorman   = 0;
  path->rs_int    = 0;
  path->vint      = 0;
  strcpy(path->exch,    "");
  strcpy(path->version, "");

  for (i = 0; i < 3; i++) {
    path->evec[i]  = 0;
    path->xivec[i] = 0;
  }
  for (i = 0; i <= nphx; i++) {
    path->iz[i] = 0;
  }
  for (i = 0; i < legtot; i++) {
    path->ipot[i] = 0;
    path->ri[i]   = 0;
    path->beta[i] = 0;
    path->eta[i]  = 0;
  }
  path->ipot[legtot] = 0;
  path->beta[legtot] = 0;
  path->eta[legtot]  = 0;
  path->eta[legtot+1]  = 0;
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

  strcpy(path->errormessage, "");
  /* strcpy(path->phpad, "phase.pad"); */
}


_EXPORT(int) make_path(FEFFPATH *path) {
  /* Compute the path using the content of the FEFFPATH struct.*/
  /* Fill the struct with geometry and F_eff information. */
  /* Return an error code. */
  int i, j;
  int error = 0;

  /* scattering and path geometry */
  int index, iorder, nleg;
  double degen;
  // int ipot[legtot+1], iz[nphx+1];
  // double rat[legtot+2][3];
  // double ri[legtot], beta[legtot+1], eta[legtot+2];
  int *ipot, *iz;
  double **rat;
  double *ri, *beta, *eta;

  /* potentials parameters */
  int ixc;
  double rs, vint, mu, edge, kf, rnrmav, gamach;

  /* onepath.f output */
  int nnnn, xdi, verbose;

  /* feffNNNN.dat columns */
  int ne;
  // double k[nex], real_phc[nex], mag_feff[nex], pha_feff[nex], red_fact[nex], lam[nex], rep[nex];
  double *k, *real_phc, *mag_feff, *pha_feff, *red_fact, *lam, *rep;
  // double k[nex], real_phc[nex], mag_feff[nex], pha_feff[nex], red_fact[nex], lam[nex], rep[nex];

  /* polarization and ellipticity */
  int ipol;
  double elpty;
  // double evec[3], xivec[3];
  double *evec, *xivec;

  char phpad[257] = {'\0'};
  char exch[9] = {'\0'};
  char version[31] = {'\0'};
  FILE *ifp;

  /* printf("entering make_path\n"); */
  /* fflush(stdout); */
  strcpy(phpad, path->phpad);
  iorder = path->iorder;
  index = path->index;
  nleg = path->nleg;
  degen = path->degen;
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
  xdi     = path->xdi;
  verbose = path->verbose;

  if ( (rat[1][0] == 0) && (rat[1][1] == 0) && (rat[1][2] == 0) ) {
    error = error + ERR_FIRSTISABS;
  };
  if ( (rat[nleg-1][0] == 0) && (rat[nleg-1][1] == 0) && (rat[nleg-1][2] == 0) ) {
    error = error + ERR_NLEGISABS;
  };
  if (degen < 0) { /* degeneracy must be non-negative */
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
  ifp = fopen(phpad, "r");  /* cannot find or read phase.pad */
  if( ifp == NULL ) {
    error = error + ERR_NOPHPAD;
  } else {
    fclose(ifp);
  };
  path->errorcode = error;
  make_path_errorstring(path);
  /* printf("%d\n", path->errorcode); */
  /* fflush(stdout); */
  if (error > 0) {
    return error;
  };

  /* printf(">%s<\n", phpad); */
  /* fflush(stdout); */
  onepath_(phpad, &index, &nleg, &degen, &iorder,
           exch, &rs, &vint, &mu, &edge, &kf, &rnrmav, &gamach,
           version, &ipot, rat, &iz, &ipol, &evec, &elpty, &xivec,
           &nnnn, &xdi, &verbose, &ri, &beta, &eta,
           &ne, &k, &real_phc, &mag_feff, &pha_feff, &red_fact, &lam, &rep);
  /* printf("after onepath_\n"); */
  /* fflush(stdout); */

  /* --------------------------------------------------- */
  /* transfer everything into the struct                 */

  /* potentials parameters */
  path->rs_int = rs;
  path->vint = vint * hart;
  path->gam_ch = gamach;
  path->kf = kf / bohr;
  path->edge = edge * hart;
  path->mu = mu * hart;
  path->rnorman = rnrmav;
  strncpy(path->exch, exch, 8);
  strncpy(path->version, version, 30);

  /* path geometry */
  for (i = 0; i <= nphx; i++) {
    path->iz[i] = iz[i];
  }
  for (i = 0; i < legtot; i++) {
    path->ri[i]   = ri[i]   * bohr;
    path->beta[i] = beta[i] * 180 / pi;
    path->eta[i]  = eta[i]  * 180 / pi;
  }
  path->beta[legtot]  = beta[legtot];
  path->eta[legtot]   = eta[legtot];
  path->eta[legtot+1] = eta[legtot+1];
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
    path->rep[i]      = rep[i];
    /* printf(" %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e\n", */
    /* 	   path->kgrid[i], path->caps[i], path->amff[i], path->phff[i], path->redfac[i], path->lambda[i], path->rep[i]); */
  }
  /* --------------------------------------------------- */
  /* printf("after arrays\n"); */
  /* fflush(stdout); */

  return error;
}

_EXPORT(int) add_scatterer(FEFFPATH *path, double x, double y, double z, int ip) {
  /* Add a scattering atom to the path. Fill the nleg member. */
  /* Return an error code. */
  int error;
  double length;
  int nleg = path->nleg;
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
  /* one needs to explicitly free each part of the struct */
  int i;

  for (i = 0; i < legtot+2; i++) {
    free(path->rat[i]);
  }
  free(path->rat);
  free(path->ipot);
  free(path->iz);
  free(path->ri);
  free(path->beta);
  free(path->eta);

  free(path->evec);
  free(path->xivec);

  free(path->k);
  free(path->real_phc);
  free(path->mag_feff);
  free(path->pha_feff);
  free(path->red_fact);
  free(path->lam);
  free(path->rep);

  free(path->errormessage);
  free(path->phpad);
  free(path->exch);
  free(path->version);

  free(path);
}

_EXPORT(double) leglength(FEFFPATH *path) {
  double x1, y1, z1, x2, y2, z2;
  int nleg = path->nleg;
  x1 = path->rat[nleg-2][0];
  y1 = path->rat[nleg-2][1];
  z1 = path->rat[nleg-2][2];
  x2 = path->rat[nleg-1][0];
  y2 = path->rat[nleg-1][1];
  z2 = path->rat[nleg-1][2];
  return sqrt( pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2) );
}


/* error string interpretation */
_EXPORT(void) make_scatterer_errorstring(FEFFPATH *path) {
  double x, y, z;
  int ip;
  char message[500];
  char error[100];
  int errcode = path->errorcode;
  int nleg = path->nleg;
  x  = path->rat[nleg-1][0];
  y  = path->rat[nleg-1][1];
  z  = path->rat[nleg-1][2];
  ip = path->ipot[nleg-1];

  if (errcode == 0) { return; }
  sprintf(message, "Error in add_scatterer at atom (%.5f, %.5f, %.5f, %d):\n", x, y, z, ip);
  if (errcode & ERR_NEGIPOT) {
    strcpy(error, "\t(code 1) ipot argument to add_scatterer is less than 0\n");
    strcat(message, error);
  };
  if (errcode & ERR_BIGIPOT) {
    strcpy(error, "\t(code 2) ipot argument to add_scatterer is greater than 7\n");
    strcat(message, error);
  } ;
  if (errcode & ERR_TOOCLOSE) {
    strcpy(error, "\t(code 4) coordinates are for an atom too close to the previous atom in the path\n");
    strcat(message, error);
  };
  if (errcode & ERR_TOOMANYLEGS) {
    strcpy(error, "\t(code 8) nlegs greater than legtot\n");
    strcat(message, error);
  };
  strcpy(path->errormessage, message);
}

_EXPORT(void) make_path_errorstring(FEFFPATH *path) {
  double degen, elpty;
  int index, iorder;
  char message[500];
  char error[100];
  int errcode = path->errorcode;
  int nleg = path->nleg;
  char phpad[256] = {'\0'};
  degen  = path->degen;
  index  = path->index;
  iorder = path->iorder;
  elpty  = path->elpty;
  strcpy(phpad, path->phpad);

  if (errcode == 0) { return; }
  sprintf(message, "Error in make_path\n");
  if (errcode & ERR_FIRSTISABS) {
    strcpy(error, "\t(code 1) the first atom specified is the absorber\n");
    strcat(message, error);
  };
  if (errcode & ERR_NLEGISABS) {
    strcpy(error, "\t(code 2) the last atom specified is the absorber\n");
    strcat(message, error);
  };
  if (errcode & ERR_DEGNEG) {
    sprintf(error, "\t(code 4) path degeneracy (%.2f) is negative\n", degen);
    strcat(message, error);
  };
  if (errcode & ERR_BADINDEX) {
    sprintf(error, "\t(code 8) path index (%d) not between 0 and 9999\n", index);
    strcat(message, error);
  };
  if (errcode & ERR_BADELPTY) {
    sprintf(error, "\t(code 16) ellipticity (%.2f) not between 0 and 1\n", elpty);
    strcat(message, error);
  };
  if (errcode & ERR_BADIORDER) {
    sprintf(error, "\t(code 32) iorder (%d) not between 0 and 10\n", iorder);
    strcat(message, error);
  };
  if (errcode & ERR_NOPHPAD) {
    sprintf(error, "\t(code 64) phase.pad file (%s) does not exist or cannot be read\n", phpad);
    strcat(message, error);
  };
  strcpy(path->errormessage, message);
}

/* simple wrapper of Fortran onepath for better portability */
_EXPORT(void) calc_onepath(char *phpad, int *index, int *nlegs, double *degen, int *iorder,
                           char *exch, double *rs, double *vint, double *mu, double *edge,
                           double *kf, double *rnorman, double *gamach, char *version,
                           int **ipot, double **rat, int **iz,
                           int *ipol, double **evec, double *elpty, double **xivec,
                           int *nnnn_out, int *xdi_out, int *verbose,
                           double **ri, double **beta, double **eta, int *nout,
                           double **a1, double **a2, double **a3,
                           double **a4, double **a5, double **a6, double **a7) {
  return onepath_(phpad, index, nlegs, degen, iorder, exch, rs, vint, mu, edge, kf,
                  rnorman, gamach, version, ipot, rat, iz, ipol, evec, elpty, xivec,
                  nnnn_out, xdi_out, verbose,
                  ri, beta, eta, nout, a1, a2, a3, a4, a5, a6, a7);
}
