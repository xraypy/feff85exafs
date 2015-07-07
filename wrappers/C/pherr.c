#include <stdio.h>
#include <stdlib.h>

#include "feffphases.h"

int main()
{
  int i, ret;
  FEFFPHASES *phases;

  phases = malloc(sizeof(FEFFPHASES));
  ret = create_phases(phases);

  strcpy(phases->jsonfile, "../fortran/libpotphXX.json");
  ret = read_libpotph_json(phases);
  if (ret > 0) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize problem with %s\n\n", phases->jsonfile);
  }

  strcpy(phases->jsonfile, "../fortran/libpotph.json");
  ret = read_libpotph_json(phases);

  phases->nph = 12;
  ret = make_phases(phases);
  if (ret & ERR_NPHX) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize nph>nphx  %d\n\n", phases->nph);
  }
  phases->nph = 2;

  phases->nat = 1500;
  ret = make_phases(phases);
  if (ret & ERR_NATX) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize nat>natx  %d\n\n", phases->nat);
  }
  phases->nat = 177;

  phases->ihole = 15;
  ret = make_phases(phases);
  if (ret & ERR_IHOLE) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad ihole  %d\n\n", phases->ihole);
  }
  phases->ihole = 1;

  phases->iz[1] = 105;
  ret = make_phases(phases);
  if (ret & ERR_Z) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad Z number  %d\n\n", phases->iz[1]);
  }
  phases->iz[1] = 2;

  phases->lmaxsc[1] = 8;
  ret = make_phases(phases);
  if (ret & ERR_L) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad angular momentum  %d\n\n", phases->lmaxsc[1]);
  }
  phases->lmaxsc[1] = 2;


  phases->xnatph[1] = -7;
  ret = make_phases(phases);
  if (ret & ERR_STOI) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad stoichiometry  %lf\n\n", phases->xnatph[1]);
  }
  phases->xnatph[1] = 100;



  phases->folp[1] = -2.0;
  ret = make_phases(phases);
  if (ret & ERR_FOLP) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad overlap  %lf\n\n", phases->folp[1]);
  }
  phases->folp[1] = 1.15;






  phases->ca = 1.0;
  ret = make_phases(phases);
  if (ret & ERR_CA) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad ca  %lf\n\n", phases->ca);
  }
  phases->ca = 0.0;


  phases->ecv = 1.0;
  ret = make_phases(phases);
  if (ret & ERR_ECV) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad ecv  %lf\n\n", phases->ecv);
  }
  phases->ecv = -40.0;



  /* phases->rscf = -0.5; */
  /* ret = make_phases(phases); */
  /* if (ret & ERR_RSCF) { */
  /*   printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode); */
  /* } else { */
  /*   printf("*** failed to recognize bad rscf  %lf\n\n", phases->rscf); */
  /* } */
  /* phases->rscf = -1.0; */



  phases->ixc = 7;
  ret = make_phases(phases);
  if (ret & ERR_IXC) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad rgrd  %d\n\n", phases->ixc);
  }
  phases->ixc = 0;



  phases->rgrd = -0.05;
  ret = make_phases(phases);
  if (ret & ERR_RGRD) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad rgrd  %lf\n\n", phases->rgrd);
  }
  phases->rgrd = 0.05;



  phases->iphat[17] = 3;
  ret = make_phases(phases);
  if (ret & ERR_IPOT) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize ipot too large  %d\n\n", phases->iphat[17]);
  }
  phases->iphat[17] = 0;
  ret = make_phases(phases);
  if (ret & ERR_IPOT) {
    printf("%s (error code %d)\n\n", phases->errormessage, phases->errorcode);
  } else {
    printf("*** failed to recognize bad ipot  %d\n\n", phases->iphat[17]);
  }
  phases->iphat[17] = 1;





  cleanup(phases);
  free(phases);
  return 0;
}
