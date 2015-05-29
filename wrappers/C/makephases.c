#include <stdio.h>
#include <stdlib.h>

#include "feffphases.h"

int main()
{
  int i, ret;
  FEFFPHASES *phases;

  phases = malloc(sizeof(FEFFPHASES));
  ret = create_phases(phases);

  strcpy(phases->jsonfile, "../fortran/libpotph.json");
  read_libpotph_json(phases);

  /* printf("reading from >%s<\n", phases->jsonfile); */
  /* printf("nat=%d\n", phases->nat); */
  /* printf("gamach=%f\n", phases->gamach); */
  /* printf("%i  %f  %f  %f  %d\n",   0, phases->rat[0][0],   phases->rat[0][1],   phases->rat[0][2],   phases->iphat[0]); */
  /* printf("%i  %f  %f  %f  %d\n",   1, phases->rat[1][0],   phases->rat[1][1],   phases->rat[1][2],   phases->iphat[1]); */
  /* printf("%i  %f  %f  %f  %d\n",   2, phases->rat[2][0],   phases->rat[2][1],   phases->rat[2][2],   phases->iphat[2]); */
  /* printf("%i  %f  %f  %f  %d\n",  23, phases->rat[23][0],  phases->rat[23][1],  phases->rat[23][2],  phases->iphat[23]); */
  /* printf("%i  %f  %f  %f  %d\n", 176, phases->rat[176][0], phases->rat[176][1], phases->rat[176][2], phases->iphat[176]); */

  /* printf(">%f<   >%f<\n", phases->folp[0], phases->folp[1]); */
  /* printf(">%s<   >%s<\n", phases->potlbl[0], phases->potlbl[1]); */
  /* printf(">%s<\n>%s<\n", phases->titles[0], phases->titles[1]); */

  make_phases(phases);



  cleanup(phases);
  free(phases);
  return 0;
}
