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

#include "feffphases.h"
#include "nxjson.c"

_EXPORT(int) create_phases(FEFFPHASES *phases) {
  /* Instantiate and initialize a FEFFPHASES struct */
  /* Return an error code -- currently hardwired to return 0. */
  int i;

  char message[500] = {'\0'};
  char jsonfl[257]  = {'\0'};
  char titl[81]     = {'\0'};
  char ptlb[7]      = {'\0'};

  strcpy(message, " ");
  strcpy(jsonfl,  " ");
  strcpy(titl,    " ");
  strcpy(ptlb,    " ");

  /* ints and doubles */
  phases->errorcode	= 0;
  phases->ntitle	= 0;
  phases->nat		= 0;
  phases->nph		= 0;
  phases->ihole		= 1;
  phases->rscf		= 0.0;
  phases->lscf		= 0;
  phases->nscmt		= 0;
  phases->ca		= 0.0;
  phases->nmix		= 0;
  phases->ecv		= 0.0;
  phases->icoul		= 0;
  phases->ipol		= 0;
  phases->elpty		= 0.0;
  phases->ispin		= 0;
  phases->angks		= 0.0;
  phases->gamach	= 0.0;
  phases->ixc		= 0;
  phases->vr0		= 0.0;
  phases->vi0		= 0.0;
  phases->ixc0		= 0;
  phases->iafolp	= 0;
  phases->rgrd		= 0.0;
  phases->iunf		= 0;
  phases->inters	= 0;
  phases->totvol	= 0.0;
  phases->jumprm	= 0;
  phases->nohole	= 0;

  /* atom cluster */
  phases->rat  = calloc(natx, sizeof(double *));
  for (i = 0; i < natx; i++) {
    phases->rat[i] = calloc(3, sizeof(double));
  }
  phases->iphat  = calloc(natx, sizeof(int));

  /* properties of unique potentials */
  phases->iz     = calloc(nphx+1, sizeof(int));
  phases->lmaxsc = calloc(nphx+1, sizeof(int));
  phases->lmaxph = calloc(nphx+1, sizeof(int));
  phases->xnatph = calloc(nphx+1, sizeof(double));
  phases->spinph = calloc(nphx+1, sizeof(double));
  phases->folp   = calloc(nphx+1, sizeof(double));
  phases->xion   = calloc(nphx+1, sizeof(double));

  /* 3 vectors */
  phases->evec   = calloc(3, sizeof(double));
  phases->xivec  = calloc(3, sizeof(double));
  phases->spvec  = calloc(3, sizeof(double));

  /* polarization tensor */
  phases->ptz  = calloc(3, sizeof(double complex *));
  for (i = 0; i < 3; i++) {
    phases->ptz[i] = calloc(3, sizeof(double complex));
  }

  /* strings */
  phases->titles = malloc(sizeof(char *) *  nheadx);
  for (i = 0; i < nheadx; i++) {
    phases->titles[i] = malloc(sizeof(char *) *  81);
    strcpy(phases->titles[i], titl);
  }
  phases->potlbl = malloc(sizeof(char *) * (nphx+1));
  for (i = 0; i < nphx+1; i++) {
    phases->potlbl[i] = malloc(sizeof(char *) * 7);
    strcpy(phases->potlbl[i], ptlb);
  }

  phases->jsonfile = calloc(257, sizeof(char));
  strcpy(phases->jsonfile, jsonfl);

  phases->errormessage = calloc(257, sizeof(char));
  strcpy(phases->errormessage, message);


  return 0;
}


_EXPORT(void) clear_phases(FEFFPHASES *phases) {
  /* Reinitialize a FEFFPHASES struct, returning everything to default */
  int i,j;

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
  phases->ixc    = 0;
  phases->vr0    = 0.0;
  phases->vi0    = 0.0;
  phases->ixc0   = 0;
  phases->iafolp = 0;
  phases->rgrd   = 0.0;
  phases->iunf   = 0;
  phases->inters = 0;
  phases->totvol = 0.0;
  phases->jumprm = 0;
  phases->nohole = 0;


  for (i = 0; i <= nheadx; i++) {
    strcpy(phases->titles[i], " ");
  }

  for (i = 0; i <= nphx+1; i++) {
    phases->iz[i]     = 0;
    phases->lmaxsc[i] = 0;
    phases->lmaxph[i] = 0;
    phases->xnatph[i] = 0.0;
    phases->spinph[i] = 0.0;
    phases->folp[i]   = 0.0;
    phases->xion[i]   = 0.0;
    strcpy(phases->potlbl[i], " ");
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

  phases->errorcode = 0;
  strcpy(phases->errormessage, " ");
  strcpy(phases->jsonfile, " ");

}


_EXPORT(int) read_libpotph_json(FEFFPHASES *phases) {
  /* read the libpotph.json file written by rdinp */

  int i;
  char string[81] = {'\0'};

  const nx_json* json;
  const nx_json* arr;
  const nx_json* item;
  const nx_json* rr;
  const nx_json* ii;

  char *buffer;
  FILE *fh = fopen(phases->jsonfile, "rb");

  if ( fh != NULL )  {
    fseek(fh, 0L, SEEK_END);
    int s = ftell(fh);
    rewind(fh);
    buffer = malloc(s);
    if ( buffer != NULL ) {
      fread(buffer, s, 1, fh);
      // we can now close the file
      fclose(fh); fh = NULL;
 
      // do something, e.g.
      // printf("%s", buffer);
      json = nx_json_parse(buffer, 0);
      if (json) {
	/* integers */
	phases->ntitle = nx_json_get(json, "ntitle")->int_value;
	phases->nat    = nx_json_get(json, "natt"  )->int_value;
	phases->nph    = nx_json_get(json, "nph"   )->int_value;
	phases->ihole  = nx_json_get(json, "ihole" )->int_value;
	phases->lscf   = nx_json_get(json, "lfms1" )->int_value;
	phases->nscmt  = nx_json_get(json, "nscmt" )->int_value;
	phases->nmix   = nx_json_get(json, "nmix"  )->int_value;
	phases->icoul  = nx_json_get(json, "icoul" )->int_value;
	phases->ipol   = nx_json_get(json, "ipol"  )->int_value;
	phases->ispin  = nx_json_get(json, "ispin" )->int_value;
	phases->ixc    = nx_json_get(json, "ixc"   )->int_value;
	phases->ixc0   = nx_json_get(json, "ixc0"  )->int_value;
	phases->iafolp = nx_json_get(json, "iafolp")->int_value;
	phases->iunf   = nx_json_get(json, "iunf"  )->int_value;
	phases->inters = nx_json_get(json, "inters")->int_value;
	phases->jumprm = nx_json_get(json, "jumprm")->int_value;
	phases->nohole = nx_json_get(json, "nohole")->int_value;

	/* doubles */
	phases->rscf   = nx_json_get(json, "rfms1" )->dbl_value;
	phases->ca     = nx_json_get(json, "ca1"   )->dbl_value;
	phases->ecv    = nx_json_get(json, "ecv"   )->dbl_value;
	phases->elpty  = nx_json_get(json, "elpty" )->dbl_value;
	phases->angks  = nx_json_get(json, "angks" )->dbl_value;
	phases->gamach = nx_json_get(json, "gamach")->dbl_value;
	phases->vr0    = nx_json_get(json, "vro"   )->dbl_value;
	phases->vi0    = nx_json_get(json, "vio"   )->dbl_value;
	phases->rgrd   = nx_json_get(json, "rgrd"  )->dbl_value;
	phases->totvol = nx_json_get(json, "totvol")->dbl_value;

	/* load matrix of cartesian coordinates */
	arr = nx_json_get(json, "x");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->rat[i][0] = item->dbl_value;
	}
	arr = nx_json_get(json, "y");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->rat[i][1] = item->dbl_value;
	}
	arr = nx_json_get(json, "z");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->rat[i][2] = item->dbl_value;
	}
	/* and the corresponding potential indeces */
	arr = nx_json_get(json, "iphatx");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->iphat[i] = item->int_value;
	}
    

	/* all the arrays that are of length nphx+1 */
	arr = nx_json_get(json, "iz");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->iz[i] = item->int_value;
	}
	arr = nx_json_get(json, "lmaxsc");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->lmaxsc[i] = item->int_value;
	}
	arr = nx_json_get(json, "lmaxph");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->lmaxph[i] = item->int_value;
	}
	arr = nx_json_get(json, "xnatph");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->xnatph[i] = item->dbl_value;
	}
	arr = nx_json_get(json, "spinph");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->spinph[i] = item->dbl_value;
	}
	arr = nx_json_get(json, "folp");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->folp[i] = item->dbl_value;
	}
	arr = nx_json_get(json, "xion");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->xion[i] = item->dbl_value;
	}


	/* 3 vectors */
	arr = nx_json_get(json, "evec");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->evec[i] = item->dbl_value;
	}
	arr = nx_json_get(json, "xivec");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->xivec[i] = item->dbl_value;
	}
	arr = nx_json_get(json, "spvec");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  phases->spvec[i] = item->dbl_value;
	}

	/* load polarization tensor */
	arr = nx_json_get(json, "ptz1");
	rr=nx_json_item(arr, 0);
	ii=nx_json_item(arr, 1);
	phases->ptz[0][0] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 2);
	ii=nx_json_item(arr, 3);
	phases->ptz[0][1] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 4);
	ii=nx_json_item(arr, 5);
	phases->ptz[0][2] = rr->dbl_value + ii->dbl_value * I;

	arr = nx_json_get(json, "ptz2");
	rr=nx_json_item(arr, 0);
	ii=nx_json_item(arr, 1);
	phases->ptz[1][0] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 2);
	ii=nx_json_item(arr, 3);
	phases->ptz[1][1] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 4);
	ii=nx_json_item(arr, 5);
	phases->ptz[1][2] = rr->dbl_value + ii->dbl_value * I;

	arr = nx_json_get(json, "ptz3");
	rr=nx_json_item(arr, 0);
	ii=nx_json_item(arr, 1);
	phases->ptz[2][0] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 2);
	ii=nx_json_item(arr, 3);
	phases->ptz[2][1] = rr->dbl_value + ii->dbl_value * I;
	rr=nx_json_item(arr, 4);
	ii=nx_json_item(arr, 5);
	phases->ptz[2][2] = rr->dbl_value + ii->dbl_value * I;


	/* title strings */
	/***********************************************************************************/
        /* these need to be written to a string of length 80 (with trailing blanks)	   */
	/* in order to avoid a valgrind warning about "Conditional jump or move		   */
        /* depends on uninitialised value(s)" coming from istrln in COMMON/str.f	   */
        /***********************************************************************************/
	arr = nx_json_get(json, "titles");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  sprintf(string, "%-80s", item->text_value);
	  strcpy(phases->titles[i], string);
	}

	/* potential labels */
	arr = nx_json_get(json, "potlbl");
	for (i=0; i<arr->length; i++) {
	  item=nx_json_item(arr, i);
	  strcpy(phases->potlbl[i], item->text_value);
	}

	nx_json_free(json);
      };
      free(buffer);
    }
    if (fh != NULL) fclose(fh);
  } else {
    phases->errorcode = JSN_NOFILE;
    sprintf(phases->errormessage, "Error reading JSON file \"%s\"", phases->jsonfile);
    return JSN_NOFILE;
  }

  phases->errorcode = 0;
  strcpy(phases->errormessage, "");
  return 0;
};


_EXPORT(int) make_phases(FEFFPHASES *phases) {
  /************************************************************/
  /* Instantiate and initialize a FEFFPHASES struct	      */
  /* Conversion to code units happens here!		      */
  /* Return an error code -- currently hardwired to return 0. */
  /************************************************************/
  int i, j;

  int ntitle, nat, nph, ihole, lscf, nscmt, nmix, icoul, ipol, ispin, ixc, ixc0, iafolp, iunf, inters, jumprm, nohole;
  double rscf, ca, ecv, elpty, angks, gamach, vr0, vi0, rgrd, totvol;
  int iz[nphx+1], lmaxsc[nphx+1], lmaxph[nphx+1], iphat[natx];
  double xnatph[nphx+1], spinph[nphx+1], folp[nphx+1], xion[nphx+1];
  double rat[natx][3];
  double complex ptz[3][3];
  double evec[3], xivec[3], spvec[3];
  char potlbl[nphx+1][7];
  char titles[nheadx][81];

  ntitle	= phases->ntitle;
  nat		= phases->nat;
  nph		= phases->nph;
  ihole		= phases->ihole;
  lscf		= phases->lscf;
  nscmt		= phases->nscmt;
  nmix		= phases->nmix;
  icoul		= phases->icoul;
  ipol		= phases->ipol;
  ispin		= phases->ispin;
  ixc		= phases->ixc;
  ixc0		= phases->ixc0;
  iafolp	= phases->iafolp;
  iunf		= phases->iunf;
  inters	= phases->inters;
  jumprm	= phases->jumprm;
  nohole	= phases->nohole;
  
  rscf		= (phases->rscf) / (bohr); /* code units! */
  ca		= phases->ca;
  ecv		= (phases->ecv) / (hart); /* code units! */
  elpty		= phases->elpty;
  angks		= phases->angks;
  gamach	= (phases->gamach) / (hart); /* code units! */
  vr0		= phases->vr0;
  vi0		= phases->vi0;
  rgrd		= phases->rgrd;
  totvol	= (phases->totvol) / (bohr*bohr*bohr); /* code units! */
  /************************************************************/
  /* i don't quite understand why the parens are necessary,   */
  /* but the divisions by hart were not evaluating correctly  */
  /* without them....	 				      */
  /************************************************************/

  for (i = 0; i < ntitle; i++) {
    strcpy(titles[i], phases->titles[i]);
  };

  for (i = 0; i <= nph; i++) {	/* <= is used because the fortran arrays start at 0 not 1! -- [0:nphx] */
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
    iphat[i] = phases->iphat[i];
    for (j = 0; j < 3; j++) {
      rat[i][j] = phases->rat[i][j] / bohr; /* code units! */
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
	    &ptz, &gamach, &ixc, &vr0, &vi0, &ixc0,
	    &iafolp, &folp, &xion, &rgrd, &iunf, &inters, &totvol, &jumprm, &nohole );

  return 0;
}


_EXPORT(void) cleanup(FEFFPHASES *phases) {

  int i;

  /* atoms list */
  for (i = 0; i < natx; i++) {
    free(phases->rat[i]);
  }
  free(phases->rat);
  free(phases->iphat);

  /* info about potentials */
  free(phases->iz);
  free(phases->lmaxsc);
  free(phases->lmaxph);
  free(phases->xnatph);
  free(phases->spinph);
  free(phases->folp);
  free(phases->xion);

  /* 3 vectors */
  free(phases->evec);
  free(phases->xivec);
  free(phases->spvec);

  /* polarization tensor */
  for (i = 0; i < 3; i++) {
    free(phases->ptz[i]);
  }
  free(phases->ptz);

  /* strings: titles, potential labels, json file name */
  for (i = 0; i < nheadx; i++) {
    free(phases->titles[i]);
  }
  free(phases->titles);
  for (i = 0; i < nphx+1; i++) {
    free(phases->potlbl[i]);
  }
  free(phases->potlbl);
  free(phases->jsonfile);
  free(phases->errormessage);

};
