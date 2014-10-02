#include <stdbool.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#define _EXPORT(a) __declspec(dllexport) a _stdcall
#else
#define _EXPORT(a) a
#endif


#define nex    150            /* from feff */
#define npatx  8              /* from feff */
#define legtot npatx+1        /* from feff */
#define bohr   0.529177249    /* bohr, a unit of length, in Ångström */
#define pi     3.1415926535897932384626433  /* π */


typedef struct {
  /* INPUT: structure of path                                                       */
  long index;        /* path index                                default = 9999    */
  long nleg;         /* number of legs in path                    use add_scatterer */
  double deg;        /* path degeneracy                           must be supplied  */
  double **rat;      /* cartesian positions of atoms in path      use add_scatterer */
  long *ipot;        /* unique potentials of atoms in path        use add_scatterer */
  long iorder;       /* order of approximation in genfmt          default = 2       */

  /* INPUT: output flags for saving F_eff to a file                                 */
  bool nnnn;         /* flag to write feffNNNN.dat file           default = false   */
  bool json;         /* flag to write feffNNNN.json file          default = false   */
  bool verbose;      /* flag to write screen messages             default = false   */

  /* INPUT: parameters controlling polarization                                     */
  bool ipol;         /* flag to do polarization calculation       default = false   */
  double *evec;      /* polarization vector                       default = (0,0,0) */
  double elpty;      /* ellipticity                               default = 0       */
  double *xivec;     /* direction of X-ray propagation            default = (0,0,0) */

  /* OUTPUT: geometry information (leg length, beta, eta)                           */
  double *ri;        /* leg lengths                                                 */
  double *beta;      /* beta angles                                                 */
  double *eta;       /* eta angles                                                  */
  double reff;       /* half path length                          computed from ri  */

  /* OUTPUT: columns of feffNNNN.dat                                                */ 
  long ne;           /* number of energy points actually used by Feff               */
  double *k;         /* k grid for feff path calculation   column 1 in feffNNNN.dat */
  double *real_phc;  /* central atom phase shifts          column 2 in feffNNNN.dat */ 
  double *mag_feff;  /* magnitude of F_eff                 column 3 in feffNNNN.dat */ 
  double *pha_feff;  /* phase of F_eff                     column 4 in feffNNNN.dat */ 
  double *red_fact;  /* reduction factor                   column 5 in feffNNNN.dat */ 
  double *lam;       /* mean free path                     column 6 in feffNNNN.dat */
  double *rep;       /* real part of complex momentum      column 7 in feffNNNN.dat */
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

void onepath_(long *,                   /* path index */
	      long *,                   /* nlegs */
	      double *,                /* degeneracy */
	      long *,                   /* iorder */
	      /* scattering geometry */
	      long (*)[legtot+1],       /* list of unique potentials */
	      double (*)[legtot+2][3], /* list of cartesian coordinates */
	      /* polarization and ellipticity */
	      long *,                   /* flag to compute polarization */
	      double (*)[3],           /* polarization vector */
	      double *,                /* ellipticity */
	      double (*)[3],           /* direction of travel */
	      /* output flags */
	      long *,                   /* integer flag for writing feffNNNN.dat */
	      long *,                   /* integer flag for writing feffNNNN.json */
	      long *,                   /* integer flag for writing screen messages */
	      /* path geometry */
	      double (*)[legtot],      /* Ri   */
	      double (*)[legtot+1],    /* beta */
	      double (*)[legtot+2],    /* eta  */
	      long *,                   /* number of points in kgrid */
	      /* seven columns of feffNNNN.dat file */
	      double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex], double (*)[nex]);


