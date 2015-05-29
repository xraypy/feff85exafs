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

#include <stdbool.h>
#include <string.h>
#include <complex.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#define _EXPORT(a) __declspec(dllexport) a _stdcall
#else
#define _EXPORT(a) a
#endif


#define natx   1000           /* from feff, see src/HEADERS/dim.h */
#define nphx   11             /* from feff */
#define nheadx 30             /* from feff */
#define bohr   0.529177249    /* bohr, a unit of length, in Ångström */
#define pi     3.1415926535897932384626433  /* π */
#define ryd    13.605698      /* rydberg, a unit of energy, in eV */
#define hart   2.0 * ryd      /* hartree, a unit of energy, in eV */


typedef struct {
  /* error handling*/
  int errorcode;                /* error code                                   */
  char *errormessage;           /* error message                                */
  /* json file name from rdinp */
  char   *jsonfile;
  /* TITLE */
  int    ntitle;		/* number of header lines                       */
  char   **titles;		/* (nheadx) array of header string              */
  /* ATOMS */
  int    nat;			/* number of atoms in cluster                   */
  double **rat; 		/* (3,natx) cartesian coordinates of atoms in cluster  */
  int    *iphat;		/* (natx) unique potential indeces of atoms in cluster */
  /* POTENTIALS */
  int    nph;			/* number of unique potentials                  */
  int    *iz;			/* (0:nphx) Z numbers of unique potentials      */
  char   **potlbl;		/* (0:nphx) labels of unique potentials         */
  int    *lmaxsc;		/* (0:nphx) l max for SCF for each potential    */
  int    *lmaxph;		/* (0:nphx) l max for FMS for each potential    */
  double *xnatph;		/* (0:nphx) stoichiometry of each potential     */
  double *spinph;		/* (0:nphx) spin on each unique potential       */
  /* HOLE/EDGE */
  int    ihole;			/* edge index, 1=K, 4=L3, etc                   */
  /* SCF */
  double rscf;			/* cluster radius for self-consistent calc.     */
  int    lscf;			/* 0=solid, 1=molecule                          */
  int    nscmt;			/* max number of self-consistency iterations    */
  double ca;			/* self-consistency convergence accelerator     */
  int    nmix;			/* number of mixing iterations before Broyden   */
  double ecv;			/* core/valence separation energy               */
  int    icoul;			/* obsolete param. for handling Coulomb pot.    */
  /* POLARIZATION and ELLIPTICITY */
  int    ipol;			/* 1=do polarization calculation                */
  double *evec;			/* (3) polarization array                       */
  double elpty;			/* eccentricity of ellilptical light            */
  double *xivec;		/* (3) ellipticity array                        */
  /* SPIN */
  int    ispin;			/* 1=do spin calculation                        */
  double *spvec;		/* (3) spin array                               */
  double angks;			/* angle between spin and incidient beam        */
  /* return */
  double complex **ptz;	        /* (-1:1,-1:1) polarization tensor              */
  double gamach;		/* tabulated core-hole lifetime                 */
  /* EXCHANGE */
  int    ixc;	                /* exchange index                               */
  double vr0;                   /* Fermi level offset                           */
  double vi0;                   /* constant broadening                          */
  int    ixc0;                  /*                                              */
  /* AFOLP and FOLP */
  int    iafolp;		/* 1=do automated overlapping                   */
  double *folp;			/* (0:nphx) overlapping fractions               */
  /* ION */
  double *xion;			/* (0:nphx) potential ionizations               */
  /* RGRID */
  double rgrd;			/* radial grid used for the potentials/phases   */
  /* UNFREEZEF */
  int    iunf;			/* 1=unfreeze f electrons                       */
  /* INTERSTITIAL */
  int    inters;
  double totvol;
  /* JUMPRM */
  int    jumprm;		/* 1=remove potential jumps at muffin tin radii */
  /* NOHOLE */
  int    nohole;		/* 1=compute without core-hole                  */
} FEFFPHASES;

int create_phases(FEFFPHASES*);
void clear_phases(FEFFPHASES*);
int make_phases(FEFFPHASES*);
int read_libpotph_json(FEFFPHASES*);

void libpotph_(int *,		     /* ntitle */
	       char (*)[nheadx][81], /* titles */
	       int *,		     /* nat    */
	       double (*)[natx][3],  /* rat    */
	       int (*)[natx],	     /* iphat  */
	       int *,		     /* nph    */
	       int (*)[nphx+1],      /* iz     */
	       char (*)[nphx+1][7],  /* potlbl */
	       int (*)[nphx+1],      /* lmaxsc */
	       int (*)[nphx+1],      /* lmaxph */
	       double (*)[nphx+1],   /* xnatph */
	       double (*)[nphx+1],   /* spinph */
	       int *,		     /* ihole  */
	       double *,	     /* rscf   */
	       int *,		     /* lscf   */
	       int *,		     /* nscmt  */
	       double *,	     /* ca     */
	       int *,		     /* nmix   */
	       double *,	     /* ecv    */
	       int *,		     /* icoul  */
	       int *,		     /* ipol   */
	       double (*)[3],	     /* evec   */
	       double *,	     /* elpty  */
	       double (*)[3],	     /* xivec  */
	       int *,		     /* ispin  */
	       double (*)[3],	     /* spvec  */
	       double *,	     /* angks  */
	       double complex (*)[3][3], /* ptz    */
	       double *,	     /* gamach */
	       int *,		     /* ixc    */
	       double *,	     /* vr0    */
	       double *,	     /* vi0    */
	       int *,		     /* ixc0   */
	       int *,		     /* iafolp */
	       double (*)[nphx+1],   /* folp   */
	       double (*)[nphx+1],   /* xion   */
	       double *,	     /* rgrd   */
	       int *,		     /* iunf   */
	       int *,		     /* inters */
	       double *,	     /* totvol */
	       int *,		     /* jumprm */
	       int *		     /* nohole */
	       );
	       
/* json file reader error codes */
#define JSN_NOFILE           1  /* json file does not exist */
