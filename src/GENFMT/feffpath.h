
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#define _EXPORT(a) __declspec(dllexport) a _stdcall
#else
#define _EXPORT(a) a
#endif


#define nex    150            /* from feff */
#define npatx  8              /* from feff */
#define legtot npatx+1        /* from feff */
#define bohr   0.529177249    /* bohr, a unit of length, in Angstrom */
#define pi     3.1415926535897932384626433  /* π */


typedef struct {
  /* INPUT: structure of path                                                      */
  int index;        /* path index                                default = 9999    */
  int nleg;         /* number of legs in path                    must be supplied  */
  double reff;      /* half path length                          must be supplied  */
  double deg;       /* path degeneracy                           must be supplied  */
  double **rat;     /* cartesian positions of atoms in path      must be supplied  */
  int *ipot;        /* unique potentials of atoms in path        must be supplied  */
  int iorder;       /* order of approximation in genfmt          default = 2       */

  /* INPUT: output flags for saving F_eff to a file                                */
  bool nnnn;        /* flag to write feffNNNN.dat file           default = false   */
  bool json;        /* flag to write feffNNNN.json file          default = false   */

  /* INPUT: parameters controlling polarization                                    */
  bool ipol;        /* flag to do polarization calculation       default = false   */
  double *evec;     /* polarization vector                       default = (0,0,0) */
  double elpty;     /* ellipticity                               default = 0       */
  double *xivec;    /* direction of X-ray propagation            default = (0,0,0) */

  /* OUTPUT: geometry information (leg length, beta, eta)                          */
  double *ri;       /* leg lengths                                                 */
  double *beta;     /* beta angles                                                 */
  double *eta;      /* eta angles                                                  */

  /* OUTPUT: columns of feffNNNN.dat                                               */ 
  int ne;           /* number of energy points actually used by Feff               */
  double *kgrid;    /* k grid for feff path calculation   column 1 in feffNNNN.dat */
  double *caps;     /* central atom phase shifts          column 2 in feffNNNN.dat */ 
  double *amff;     /* magnitude of F_eff                 column 3 in feffNNNN.dat */ 
  double *phff;     /* phase of F_eff                     column 4 in feffNNNN.dat */ 
  double *redfac;   /* reduction factor                   column 5 in feffNNNN.dat */ 
  double *lambda;   /* mean free path                     column 6 in feffNNNN.dat */
  double *rep;      /* real part of complex momentum      column 7 in feffNNNN.dat */
} FEFFPATH;

