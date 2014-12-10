#include <stdio.h>
#include <stdlib.h>

#include "feffpath.h"

long main()
{
  long i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret = create_path(path);

  strcpy(path->phbin, "../fortran/phase.bin");
  path->nnnn    = 1;
  path->json    = 0;
  path->verbose = 1;


  /* --------- Errors in add_scatterer ----------------------------------------------------- */

  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 9); /* ipot negative, add_scatterer error 1 */
  if (path->errorcode == 0) {
    ret = make_path(path);
  } else {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* ipot too big, add_scatterer error 2 */
  if (path->errorcode == 0) {
    ret = make_path(path);
  } else {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->deg     = 12.0;
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
  path->deg     = 12.0;
  ret = add_scatterer(path,  0,     0,  0,     1); /* first atom absorber, make_path error 1 */
  ret = add_scatterer(path,  1.805, 0,  1.805, 1);
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 1);
  ret = add_scatterer(path,  0,     0,  0,     1); /* last atom absorber, make_path error 2 */
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  path->index   = 1;
  path->deg     = -12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* negative degeneracy, make_path error 4 */
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);



  path->index   = 40000;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* bad index, make_path error 8 */
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);



  path->index   = 1;
  path->deg     = 12.0;
  path->iorder   = -1;
  ret = add_scatterer(path,  1.805, 0,  1.805, -1); /* bad iorder, make_path error 32 */
  ret = make_path(path);
  if (path->errorcode != 0) {
    printf("%s\n", path->errormessage);
  };

  clear_path(path);


  /* this is not working yet.... */
  /* strcpy(path->phbin, "foo.bar"); */
  /* path->index   = 1; */
  /* path->deg     = 12.0; */
  /* ret = add_scatterer(path,  1.805, 0,  1.805, -1); /\* bad phbin, make_path error 64 *\/ */
  /* ret = make_path(path); */
  /* if (path->errorcode != 0) { */
  /*   printf("%s\n", path->errormessage); */
  /* }; */

  /* clear_path(path); */



  cleanup(path);
  return 0;
}
