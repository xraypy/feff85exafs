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


#define nex    150            /* from feff */
#define nphx   11             /* from feff */
#define npatx  8              /* from feff */
#define legtot npatx+1        /* from feff */
#define bohr   0.529177249    /* bohr, a unit of length, in Ångström */
#define pi     3.1415926535897932384626433  /* π */
#define ryd    13.605698      /* rydberg, a unit of energy, in eV */
#define hart   2.0 * ryd      /* hartree, a unit of energy, in eV */

typedef struct {
  /* INPUT: path to phase.bin file                                                   */
  char *phbin;

  /* INPUT: structure of path                                                        */
  long index;         /* path index                                default = 9999    */
  long nleg;          /* number of legs in path                    use add_scatterer */
  double degen;       /* path degeneracy                           must be supplied  */
  double **rat;       /* cartesian positions of atoms in path      use add_scatterer */
  long *ipot;         /* unique potentials of atoms in path        use add_scatterer */
  long iorder;        /* order of approximation in genfmt          default = 2       */

  /* INPUT: output flags for saving F_eff to a file                                  */
  bool nnnn;          /* flag to write feffNNNN.dat file           default = false   */
  bool json;          /* flag to write feffNNNN.json file          default = false   */
  bool verbose;       /* flag to write screen messages             default = false   */

  /* INPUT: parameters controlling polarization                                      */
  bool ipol;          /* flag to do polarization calculation       default = false   */
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
  long *iz;           /* atomic numbers of atoms in path     obtained from phase.bin */
  double *ri;         /* leg lengths                                                 */
  double *beta;       /* beta angles                                                 */
  double *eta;        /* eta angles                                                  */
  double reff;        /* half path length                          computed from ri  */

  /* OUTPUT: columns of feffNNNN.dat                                                 */ 
  long ne;            /* number of energy points actually used by Feff               */
  double *k;          /* k grid for feff path calculation   column 1 in feffNNNN.dat */
  double *real_phc;   /* central atom phase shifts          column 2 in feffNNNN.dat */ 
  double *mag_feff;   /* magnitude of F_eff                 column 3 in feffNNNN.dat */ 
  double *pha_feff;   /* phase of F_eff                     column 4 in feffNNNN.dat */ 
  double *red_fact;   /* reduction factor                   column 5 in feffNNNN.dat */ 
  double *lam;        /* mean free path                     column 6 in feffNNNN.dat */
  double *rep;        /* real part of complex momentum      column 7 in feffNNNN.dat */

  /* OUTPUT: error handling                                                          */
  long errorcode;     /* error code from add_scatterer or make_path                  */
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

long add_scatterer(FEFFPATH*, double, double, double, long);
long create_path(FEFFPATH*);
void clear_path(FEFFPATH*);
long make_path(FEFFPATH*);
void cleanup(FEFFPATH*);

void onepath_(char *,
	      long *,                   /* path index */
	      long *,                   /* nlegs */
	      double *,                 /* degeneracy */
	      long *,                   /* iorder */
	      long *,                   /* ixc, potential model index */
	      double *,                 /* rs, interstitial radius estimate */
	      double *,                 /* vint, interstitial potential energy */
	      double *,                 /* mu */
	      double *,                 /* edge */
	      double *,                 /* kf */
	      double *,                 /* rnrmav, average R_norman */
	      double *,                 /* gamach, chore hole lifetime in eV */
	      /* scattering geometry */
	      long (*)[legtot+1],       /* list of unique potentials */
	      double (*)[legtot+2][3],  /* list of cartesian coordinates */
	      long (*)[nphx+1],         /* list of atomic numbers */
	      /* polarization and ellipticity */
	      long *,                   /* flag to compute polarization */
	      double (*)[3],            /* polarization vector */
	      double *,                 /* ellipticity */
	      double (*)[3],            /* direction of travel */
	      /* output flags */
	      long *,                   /* integer flag for writing feffNNNN.dat */
	      long *,                   /* integer flag for writing feffNNNN.json */
	      long *,                   /* integer flag for writing screen messages */
	      /* path geometry */
	      double (*)[legtot],       /* Ri   */
	      double (*)[legtot+1],     /* beta */
	      double (*)[legtot+2],     /* eta  */
	      long *,                   /* number of points in kgrid */
	      /* seven columns of feffNNNN.dat file */
	      double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex]);


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
#define ERR_NOPHBIN          64  /* phase.bin file cannot be found or cannot be read */

#define COPY_STRING(dest,src)  dest=calloc(strlen(src)+1, sizeof(char));\
  strcpy(dest, src);
