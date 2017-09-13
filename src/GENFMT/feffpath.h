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

#include <stdbool.h>
#include <string.h>


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#define _EXPORT(a) __declspec(dllexport) a _stdcall
#else
#define _EXPORT(a) a
#endif


#define nex    150            /* from feff, see src/HEADERS/dim.h */
#define nphx   11             /* from feff */
#define npatx  8              /* from feff */
#define legtot npatx+1        /* from feff */
#define bohr   0.529177249    /* bohr, a unit of length, in Ångström */
#define pi     3.1415926535897932384626433  /* π */
#define ryd    13.605698      /* rydberg, a unit of energy, in eV */
#define hart   2.0 * ryd      /* hartree, a unit of energy, in eV */

typedef struct {
  /* INPUT: path to phase.pad file                                                   */
  char *phpad;

  /* INPUT: structure of path                                                       */
  int index;         /* path index                                default = 9999    */
  int nleg;          /* number of legs in path                    use add_scatterer */
  double degen;      /* path degeneracy                           must be supplied  */
  double **rat;      /* cartesian positions of atoms in path      use add_scatterer */
  int *ipot;         /* unique potentials of atoms in path        use add_scatterer */
  int iorder;        /* order of approximation in genfmt          default = 2       */

  /* INPUT: output flags for saving F_eff to a file                                 */
  int nnnn;          /* flag to write feffNNNN.dat file           default = false   */
  int xdi;           /* flag to write feffNNNN.xdi file           default = false   */
  int verbose;       /* flag to write screen messages             default = false   */

  /* INPUT: parameters controlling polarization                                      */
  int ipol;           /* flag to do polarization calculation       default = false   */
  double *evec;       /* polarization vector                       default = (0,0,0) */
  double elpty;       /* ellipticity                               default = 0       */
  double *xivec;      /* direction of X-ray propagation            default = (0,0,0) */

  /* OUTPUT: various strings and physical constants from feffNNNN.dat header         */
  double edge;        /* energy threshold relative to atomic value (a poor estimate) */
  double gam_ch;      /* core level energy width                                     */
  double kf;          /* k value at Fermi level                                      */
  double mu;          /* Fermi level, eV                                             */
  double rnorman;     /* Norman radius                                               */
  double rs_int;      /* interstitial radius                                         */
  double vint;        /* interstitial potential                                      */
  char *exch;         /* string describing electronic exchange model                 */
  char *version;      /* Feff version                                                */

  /* OUTPUT: geometry information (leg length, beta, eta, Z)                         */
  int *iz;            /* atomic numbers of atoms in path     obtained from phase.pad */
  double *ri;         /* leg lengths                                                 */
  double *beta;       /* beta angles                                                 */
  double *eta;        /* eta angles                                                  */
  double reff;        /* half path length                          computed from ri  */

  /* OUTPUT: columns of feffNNNN.dat                                                 */
  int ne;             /* number of energy points actually used by Feff               */
  double *k;          /* k grid for feff path calculation   column 1 in feffNNNN.dat */
  double *real_phc;   /* central atom phase shifts          column 2 in feffNNNN.dat */
  double *mag_feff;   /* magnitude of F_eff                 column 3 in feffNNNN.dat */
  double *pha_feff;   /* phase of F_eff                     column 4 in feffNNNN.dat */
  double *red_fact;   /* reduction factor                   column 5 in feffNNNN.dat */
  double *lam;        /* mean free path                     column 6 in feffNNNN.dat */
  double *rep;        /* real part of complex momentum      column 7 in feffNNNN.dat */

  /* OUTPUT: error handling                                                          */
  int errorcode;      /* error code from add_scatterer or make_path                  */
  char *errormessage; /* error code from add_scatterer or make_path                  */
} FEFFPATH;

/* --------------------------------------------------------------------------------------------------------------- */
/* still need to capture the following items in Larch's _feffdat group:                                            */
/* (see http://xraypy.github.io/xraylarch/xafs/feffpaths.html#the-feffdat-group-full-details-of-the-feff-dat-file) */
/*    - edge          energy threshold relative to atomic valu (a poor estimate)                                   */
/*    - exch          string describing electronic exchange model                                                  */
/*    - gam_ch        core level energy width                                                                      */
/*    - kf            k value at Fermi level                                                                       */
/*    - mu            Fermi level, eV                                                                              */
/*    - potentials    path potentials: list of (ipot, z, r_MuffinTin, r_Norman)                                    */
/*    - rnorman       Norman radius                                                                                */
/*    - rs_int        interstitial radius                                                                          */
/*    - title         user title                                                                                   */
/*    - version       Feff version                                                                                 */
/*    - vint          interstitial potential                                                                       */
/* --------------------------------------------------------------------------------------------------------------- */

_EXPORT(int) add_scatterer(FEFFPATH*, double, double, double, int);
_EXPORT(int) create_path(FEFFPATH*);
_EXPORT(void) clear_path(FEFFPATH*);
_EXPORT(int) make_path(FEFFPATH*);
_EXPORT(void) cleanup(FEFFPATH*);
_EXPORT(void) make_path_errorstring(FEFFPATH*);
_EXPORT(void) make_scatterer_errorstring(FEFFPATH*);
_EXPORT(double) leglength(FEFFPATH*);

_EXPORT(void) calc_onepath(char *,    /* phpad char[257],  path to phase.pad file */
			   int *,     /* index,            path index */
			   int *,     /* nlegs,            number of path legs */
			   double *,  /* degen,            degeneracy */
			   int *,     /* iorder,           exchange order */
			   char *,    /* exch char[9],     potential model description */
			   double *,  /* rs,               interstitial radius estimate */
			   double *,  /* vint,             interstitial potential energy */
			   double *,  /* mu,               Fermi level in eV */
			   double *,  /* edge,             estimate of energy threshold */
			   double *,  /* kf                k value of Fermi level */
			   double *,  /* rnorman,          average R_norman */
			   double *,  /* gamach,           core hole lifetime in eV */
			   char *,    /* version char[31], versioning string */

			   /* scattering geometry */
			   int **,    /* iz [nphx+1],      list of unique potentials */
			   double **, /* rat [legtot+2,3], list of cartesian coordinates */
			   int **,    /* ipot [legtot+1],  list of atomic numbers */

			   /* polarization and ellipticity */
			   int *,     /* ipol,             flag to compute polarization */
			   double **, /* evec [3],         polarization vector */
			   double *,  /* elpty,            ellipticity */
			   double **, /* xivec [3].        direction of travel */

			   /* flags controlling output */
			   int *,     /* nnnn_out,         flag for writing feffNNNN.dat */
			   int *,     /* xdi_out,          flag for writing feffNNNN.xdi */
			   int *,     /* verbose,          flag for writing screen messages */

			   /* path geometry */
			   double **, /* ri [legtot],      ri   */
			   double **, /* beta [legtot+1],  beta*/
			   double **, /* eta [legtot+2],   eta */
			   int *,     /* ne,               number of energy/k points */

			   /* output arrays corresponding to columns of feffNNNN.dat file */
			   double **, /* k [nex],          k grid for other arrays */
			   double **, /* real_phc [nex],   central atom phase shifts  */
			   double **, /* mag_feff [nex],   magnitude of F_eff  */
			   double **, /* pha_feff [nex],   phase of F_eff  */
			   double **, /* red_fact [nex],   reduction factor */
			   double **, /* lam [nex],        mean free path */
			   double **  /* rep [nex],        real part of complex momentum */
			   );

/* see calc_onepath for details of arg list */
void onepath_(char *, int *, int *, double *, int *, char *, double *,
	      double *, double *, double *, double *, double *, double *,
	      char *, int **, double **, int **, int *, double **, double *,
	      double **, int *, int *, int *, double **, double **,
	      double **, int *, double **, double **, double **, double **,
	      double **, double **, double **);

/* add_scatterer error codes */
#define ERR_NEGIPOT           1  /* ipot argument to add_scatterer lt 0 */
#define ERR_BIGIPOT           2  /* ipot argument to add_scatterer gt 7 */
#define ERR_TOOCLOSE          4  /* coordinates are for an atom too close to the previous atom in the path */
#define ERR_TOOMANYLEGS       8  /* nlegs gt legtot */

/* make_path error codes */
#define ERR_FIRSTISABS        1  /* the first atom specified is the absorber */
#define ERR_NLEGISABS         2  /* the last atom specified is the absorber */
#define ERR_DEGNEG            4  /* degeneracy is negative */
#define ERR_BADINDEX          8  /* index lt 0 or gt 9999 */
#define ERR_BADELPTY         16  /* elpty lt 0 or gt 1 */
#define ERR_BADIORDER        32  /* iorder lt 0 or gt ? */
#define ERR_NOPHPAD          64  /* phase.pad file cannot be found or cannot be read */

#define COPY_STRING(dest,src)  dest=calloc(strlen(src)+1, sizeof(char));\
  strcpy(dest, src);
