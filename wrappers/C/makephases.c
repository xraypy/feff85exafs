#include <stdio.h>
#include <stdlib.h>

#include "feffphases.h"
#include "nxjson.c"
#include <complex.h>

long main()
{
  long i, ret;

  long ntitle, nat, nph, ihole, lfms, nscmt, nmix, icoul, ipol, ispin, ixc, ixc0, iafolp, iunf, inters, jumprm, nohole;
  double rfms, ca, ecv, elpty, angks, gamach, vro, vri, rgrd, totvol;
  char titles[nheadx][81], potlbl[nphx+1][7];
  double rat[natx][3];
  long iphat[natx];
  long iz[nphx+1], lmaxsc[nphx+1], lmaxph[nphx+1];
  double xnatph[nphx+1], spinph[nphx+1], folp[nphx+1], xion[nphx+1];
  double evec[3], xivec[3], spvec[3];
  double ptz1[6], ptz2[6], ptz3[6];
  double complex ptz[3][3];

  FEFFPHASES *phase;

  const nx_json* json;
  const nx_json* arr;
  const nx_json* item;
  const nx_json* rr;
  const nx_json* ii;

  char *buffer;
  FILE *fh = fopen("../fortran/libpotph.json", "rb");

  if ( fh != NULL )  {
    fseek(fh, 0L, SEEK_END);
    long s = ftell(fh);
    rewind(fh);
    buffer = malloc(s);
    if ( buffer != NULL ) {
      fread(buffer, s, 1, fh);
      // we can now close the file
      fclose(fh); fh = NULL;
 
      // do something, e.g.
      // printf("%s", buffer);
      json = nx_json_parse(buffer, 0);
 
      free(buffer);
    }
    if (fh != NULL) fclose(fh);
  }

  if (json) {
    /* integers */
    ntitle = nx_json_get(json, "ntitle")->int_value;
    nat    = nx_json_get(json, "natt"  )->int_value;
    nph    = nx_json_get(json, "nph"   )->int_value;
    ihole  = nx_json_get(json, "ihole" )->int_value;
    lfms   = nx_json_get(json, "lfms1" )->int_value;
    nscmt  = nx_json_get(json, "nscmt" )->int_value;
    nmix   = nx_json_get(json, "nmix"  )->int_value;
    icoul  = nx_json_get(json, "icoul" )->int_value;
    ipol   = nx_json_get(json, "ipol"  )->int_value;
    ispin  = nx_json_get(json, "ispin" )->int_value;
    ixc    = nx_json_get(json, "ixc"   )->int_value;
    ixc0   = nx_json_get(json, "ixc0"  )->int_value;
    iafolp = nx_json_get(json, "iafolp")->int_value;
    iunf   = nx_json_get(json, "iunf"  )->int_value;
    inters = nx_json_get(json, "inters")->int_value;
    jumprm = nx_json_get(json, "jumprm")->int_value;
    nohole = nx_json_get(json, "nohole")->int_value;

    /* doubles */
    rfms   = nx_json_get(json, "rfms1" )->dbl_value;
    ca     = nx_json_get(json, "ca1"   )->dbl_value;
    ecv    = nx_json_get(json, "ecv"   )->dbl_value;
    elpty  = nx_json_get(json, "elpty" )->dbl_value;
    angks  = nx_json_get(json, "angks" )->dbl_value;
    gamach = nx_json_get(json, "gamach")->dbl_value;
    vro    = nx_json_get(json, "vro"   )->dbl_value;
    vri    = nx_json_get(json, "vri"   )->dbl_value;
    rgrd   = nx_json_get(json, "rgrd"  )->dbl_value;
    totvol = nx_json_get(json, "totvol")->dbl_value;

    /* load matrix of cartesian coordinates */
    arr = nx_json_get(json, "x");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      rat[i][0] = item->dbl_value;
    }
    arr = nx_json_get(json, "y");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      rat[i][1] = item->dbl_value;
    }
    arr = nx_json_get(json, "z");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      rat[i][2] = item->dbl_value;
    }
    /* and the corresponding potential indeces */
    arr = nx_json_get(json, "iphatx");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      iphat[i] = item->int_value;
    }
    

    /* all the arrays that are nphx+1 long */
    arr = nx_json_get(json, "iz");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      iz[i] = item->int_value;
    }
    arr = nx_json_get(json, "lmaxsc");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      lmaxsc[i] = item->int_value;
    }
    arr = nx_json_get(json, "lmaxph");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      lmaxph[i] = item->int_value;
    }
    arr = nx_json_get(json, "xnatph");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      xnatph[i] = item->dbl_value;
    }
    arr = nx_json_get(json, "spinph");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      spinph[i] = item->dbl_value;
    }
    arr = nx_json_get(json, "folp");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      folp[i] = item->dbl_value;
    }
    arr = nx_json_get(json, "xion");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      xion[i] = item->dbl_value;
    }


    /* 3 vectors */
    arr = nx_json_get(json, "evec");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      evec[i] = item->dbl_value;
    }
    arr = nx_json_get(json, "xivec");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      xivec[i] = item->dbl_value;
    }
    arr = nx_json_get(json, "spvec");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      spvec[i] = item->dbl_value;
    }

    /* load polarization tensor */
    arr = nx_json_get(json, "ptz1");
    rr=nx_json_item(arr, 0);
    ii=nx_json_item(arr, 1);
    ptz[0][0] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 2);
    ii=nx_json_item(arr, 3);
    ptz[0][1] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 4);
    ii=nx_json_item(arr, 5);
    ptz[0][2] = rr->dbl_value + ii->dbl_value * I;

    arr = nx_json_get(json, "ptz2");
    rr=nx_json_item(arr, 0);
    ii=nx_json_item(arr, 1);
    ptz[1][0] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 2);
    ii=nx_json_item(arr, 3);
    ptz[1][1] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 4);
    ii=nx_json_item(arr, 5);
    ptz[1][2] = rr->dbl_value + ii->dbl_value * I;

    arr = nx_json_get(json, "ptz3");
    rr=nx_json_item(arr, 0);
    ii=nx_json_item(arr, 1);
    ptz[2][0] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 2);
    ii=nx_json_item(arr, 3);
    ptz[2][1] = rr->dbl_value + ii->dbl_value * I;
    rr=nx_json_item(arr, 4);
    ii=nx_json_item(arr, 5);
    ptz[2][2] = rr->dbl_value + ii->dbl_value * I;


    /* title strings */
    arr = nx_json_get(json, "titles");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      strcpy(titles[i], item->text_value);
    }

    /* potential labels */
    arr = nx_json_get(json, "potlbl");
    for (i=0; i<arr->length; i++) {
      item=nx_json_item(arr, i);
      strcpy(potlbl[i], item->text_value);
    }


    printf("nat=%ld\n", nat);
    printf("gamach=%f\n", nx_json_get(json, "gamach")->dbl_value);
    printf("vfeff=\"%s\"\n", nx_json_get(json, "vfeff")->text_value);
    printf("%i  %f  %f  %f  %ld\n", 0, rat[0][0], rat[0][1], rat[0][2], iphat[0]);
    printf("%i  %f  %f  %f  %ld\n", 1, rat[1][0], rat[1][1], rat[1][2], iphat[1]);
    printf("%i  %f  %f  %f  %ld\n", 2, rat[2][0], rat[2][1], rat[2][2], iphat[2]);
    printf("%i  %f  %f  %f  %ld\n", 23, rat[23][0], rat[23][1], rat[23][2], iphat[23]);
    printf("%i  %f  %f  %f  %ld\n", 176, rat[176][0], rat[176][1], rat[176][2], iphat[176]);

    printf(">%s<   >%s<\n", potlbl[0], potlbl[1]);
    printf(">%s<\n>%s<\n", titles[0], titles[1]);

    const nx_json* arr=nx_json_get(json, "xnatph");
    int i;
    for (i=0; i<arr->length; i++) {
      const nx_json* item=nx_json_item(arr, i);
      printf("arr[%d] = %lf\n", i, item->dbl_value);
    }

    nx_json_free(json);
  };

  phase = malloc(sizeof(FEFFPHASES));
  ret = create_phases(phase);


  cleanup(phase);
  free(phase);
  return 0;
}
