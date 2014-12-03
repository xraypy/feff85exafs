#if defined(SWIGPERL)
%module "Xray::FeffPathWrapper"
#else
%module feffpath
#endif

%{
/* Includes the header in the wrapper code */
#include "feffpath.h"
%}
 
/* Parse the header file to generate wrappers */
%include "feffpath.h"


%inline %{
  void set_evec(FEFFPATH *path, int index, double val) {
    if ((index < 0) || (index >= 3)) return;
    else path->evec[index] = val;
  }
  double get_evec(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= 3)) return -1;
    else return path->evec[index];
  }
  void set_xivec(FEFFPATH *path, int index, double val) {
    if ((index < 0) || (index >= 3)) return;
    else path->xivec[index] = val;
  }
  double get_xivec(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= 3)) return -1;
    else return path->xivec[index];
  }

  double get_ri(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= legtot)) return -1;
    else return path->ri[index];
  }
  double get_beta(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= legtot)) return -1;
    else return path->beta[index];
  }
  double get_eta(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= legtot)) return -1;
    else return path->eta[index];
  }

  double get_k(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->k[index];
  }
  double get_real_phc(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->real_phc[index];
  }
  double get_mag_feff(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->mag_feff[index];
  }
  double get_pha_feff(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->pha_feff[index];
  }
  double get_red_fact(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->red_fact[index];
  }
  double get_lam(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->lam[index];
  }
  double get_realp(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->realp[index];
  }
%}
