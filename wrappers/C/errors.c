#include <stdio.h>
#include <stdlib.h>

#include "feffpath.h"

long main()
{
  long i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret  = create_path(path);

  strcpy(path->phpad, "../fortran/phase.pad");
  path->nnnn    = 1;
  /* path->json    = 0; */
  path->verbose = 1;


  /* --------- Errors in add_scatterer ----------------------------------------------------- */

  path->index   = 1;
  path->degen   = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* ipot negative, add_scatterer error 1 */
  if (path->errorcode == 0) {
    ret = make_path(path);
  } else {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->degen   = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 9); /* ipot too big, add_scatterer error 2 */
  if (path->errorcode == 0) {
    ret = make_path(path);
  } else {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->degen   = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 1); /* atoms too close, add_scatterer error 4 */
  ret = add_scatterer(path,  1.805, 0,  1.905, 1);
  if (path->errorcode == 0) {
    ret = make_path(path);
  } else {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  /* --------- Errors in make_path ----------------------------------------------------- */


  path->index   = 1;
  path->degen   = 12.0;
  ret = add_scatterer(path,  0,     0,  0,     1); /* first atom absorber, make_path error 1 */
  ret = add_scatterer(path,  1.805, 0,  1.805, 1);
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  /* recognize >1 error in the one call to make_path */
  path->index   = 1;
  path->degen   = -12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* negative degeneracy, make_path error 4 */
  ret = add_scatterer(path,  0,     0,  0,     1); /* last atom absorber, make_path error 2 */
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);



  path->index   = 40000;  /* bad index, make_path error 8 */
  path->degen   = 12.0;
  path->elpty   = -0.5;
  path->iorder  = -1;    /* bad iorder, make_path error 32 */
  strcpy(path->phpad, "foo.bar");  /* bad phpad, make_path error 64 */
  ret = add_scatterer(path,  1.805, 0,  1.805, -1);
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };


  cleanup(path);
  return 0;
}
