#include <stdio.h>
#include <stdlib.h>

#include "feffphases.h"

int main()
{
  int ret = 0;
  FEFFPHASES *phases = malloc(sizeof(FEFFPHASES));
  ret = create_phases(phases);

  strcpy(phases->jsonfile, "../fortran/libpotph.json");
  /* strcpy(phases->jsonfile, "libpotph.json"); */
  /* strcpy(phases->phpad, "foo.pad"); */
  ret = read_libpotph_json(phases);
  phases->verbose = true;
  if (ret > 0) {
    printf("%s (error code %d)\n", phases->errormessage, phases->errorcode);
  } else {
    ret = make_phases(phases);
    if (ret > 0) {
      printf("%s (error code %d)\n", phases->errormessage, phases->errorcode);
    }
  }


  cleanup(phases);
  free(phases);
  return 0;
}
