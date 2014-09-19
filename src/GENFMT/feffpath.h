#define nex    150
#define npatx  8
#define legtot npatx+1
#define bohr   0.529177249
#define pi     3.1415926535897932384626433

typedef struct {
  /* structure of path */
  int index;        /* path index */
  int nleg;         /* number of legs in path */
  double reff;      /* half path length */
  double deg;       /* path degeneracy */
  double **rat;     /* cartesian positions of atoms in path */
  int *ipot;        /* unique potentials of atoms in path */
  int iorder;       /* order of approximation in genfmt */

  /* output flags for saving F_eff to a file */
  bool nnnn;        /* flag to write feffNNNN.dat file */
  bool json;        /* flag to write feffNNNN.json file */

  /* input parameters controlling polarization */
  bool ipol;        /* flag to do polarization calculation */
  double *evec;     /* polarization vector */
  double elpty;     /* ellipticity */
  double *xivec;    /* direction of X-ray propagation */

  /* output with geometry information (leg length, beta, eta) */
  double *ri;       /* leg lengths */
  double *beta;     /* beta angles */
  double *eta;      /* eta angles  */

  /* output containing F_eff */ 
  int ne;           /* number of energy points actually used by Feff */
  double *kgrid;    /* k grid for feff path calculation, column 1 in feffNNNN.dat  */
  double *caps;     /* central atom phase shifts, column 2 in feffNNNN.dat */ 
  double *amff;     /* magnitude of F_eff, column 3 in feffNNNN.dat*/ 
  double *phff;     /* phase of F_eff, column 4 in feffNNNN.dat*/ 
  double *redfac;   /* reduction factor, column 5 in feffNNNN.dat*/ 
  double *lambda;   /* mean free path, column 6 in feffNNNN.dat*/ 
  double *rep;      /* real part of complex momentum, column 7 in feffNNNN.dat*/ 
} FEFFPATH;

