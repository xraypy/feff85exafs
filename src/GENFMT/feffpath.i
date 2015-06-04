/* ********************************************************************** */
/* LICENSE AND COPYRIGHT                                                  */
/*                                                                        */
/* To the extent possible, the authors have waived all rights granted by  */
/* copyright law and related laws for the code and documentation that     */
/* make up the SWIG Interface to the feffpath library.  While information */
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

#if defined(SWIGPERL)
%module "Xray::FeffPathWrapper"
#else
%module feffpathwrapper
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

  int get_iz(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nphx)) return -1;
    else return path->iz[index];
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
  double get_rep(FEFFPATH *path, int index) {
    if ((index < 0) || (index >= nex)) return -1;
    else return path->rep[index];
  }
%}
