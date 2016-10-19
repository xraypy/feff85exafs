# Using the feffpath C wrapper

The program `makepath.c` will read from `phase.pad` (this example is
calculated from copper metal) and write the files `feff0001.dat` and
`feff0004.dat`.

The `makepath.c` program likely will not compile before you build and
install feff85exafs.  The build script will look for the `feffpath`
library and the `feffpath.h` header file in their installation
locations

The `errors.c` program exercises most of the error conditions
resulting from improper input to `add\_scatterer` or `make_path`.

To compile the `makepath` and `errors` sample programs:

	~> make

That's it!

## Sample program

Here is the simplest program using the C wrapper:

```C
#include <stdio.h>
#include <stdlib.h>
#include "feffpath.h"

long main()
{
  long i, ret;
  FEFFPATH *path;

  path = malloc(sizeof(FEFFPATH));
  ret = create_path(path);

  strcpy(path->phpad, "../fortran/phase.pad");
  path->nnnn    = 1;
  path->xdi     = 0;
  path->verbose = 1;

  /* first path in copper */
  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path, 1.805, 0, 1.805, 1);
  if (ret > 0) {
    printf("%s", path->errormessage);
    exit(ret);
  };
  ret = make_path(path);
  if (ret > 0) {
    printf("%s", path->errormessage);
    exit(ret);
  };

  clear_path(path);

  free(path);
  return 0;
}

```

1. The `feffpath.h` header file is included.

2. A FEFFPATH struct is created and called `path`.  Memory is
   allocated for it.

3. The call to `create_path` allocates memory for all elements of the
   struct and initializes everything.  The `phpad` attribute, which
   specifies the location of the `phase.pad` file, is set to
   `phase.pad` in the current working directory.

4. `strcpy` is used to set the `phpad` attribute to the actual location
   of the `phase.pad` file.  Note that the fortran parameter in the
   `onepath` library has this string dimensioned to be 256 characters.
   And it's fortran, so it probably won't deal gracefully with unicode
   characters.

5. Setting `nnnn` and `verbose` to a true value tells the program to
   write `feffNNNN.dat` files and to write a short message to the
   screen when it is written.  Setting `xdi` to a true value will
   cause `feffNNNN.xdi` to be written.  This is equivalent to
   `feffNNNN.dat`, but in XDI format.

6. The path index is set to 1.  This means the output file will be
   called `feff0001.dat`.

7. The degeneracy of the path is specified.

8. The `add_scatterer` method is used to build the scattering
   geometry.  It's arguments are the Cartesian coordinates (referenced
   to the absorber *at the origin*), the unique potential index of the
   atom.  The return value is an error code (see `feffpath.h`).  0
   means no error.

   In this case, a single scattering path is calculated, so
   `add_scatterer` is called only once.  For a multiple scattering
   path, `add_scatterer` would be called repeatedly, once for each
   atom (or leg) in the path.  `add_scatterer` fills the `rat` and
   `ipot` members of the struct and also sets `nleg`.

9. The call to `make_path` computes the parts of F-effective and
   stores them in the `path` struct.  The return value is an error
   code (see `feffpath.h`).  0 means no error.

10. The call to `clear_path` reinitializes the struct.  This is not
	strictly necessary in this case, but would be were the program to
	go on to compute another path.  Note that clear_path does not
	reset `phpad`, `nnnn`, or `xdi`, which are assumed to be
	constant from time to time calling `make_path`.

11. Finally, the memory for `path` is deallocated and the program
    terminates.


## Contents of the FEFFPATH struct

Here are all the elements of the struct, their data types,
descriptions, and default values.  Input parameters are marked with
"I" and output parameters, set by the `make_path` method, are marked
with "O".  Most of the names are chosen to be consistent with the
naming conventions in Feff.  The output arrays for the columns of
`feffNNNN.dat` are those
[used by Larch](http://xraypy.github.io/xraylarch/xafs/feffpaths.html#the-feffdat-group-full-details-of-the-feff-dat-file).

| attribute  | type       | I/O | description                             | default              |
| ---------- | --------   | --- |---------------------------------------- | -------------------- |
|  phpad     | \*char     | I   | path to `phase.pad`                     |  `phase.pad`         |
|  index     | long       | I   | path index                              |  9999                |
|  deg       | double     | I   | path degeneracy                         |  required input      |
|  nleg      | long       | I   | number of legs in path                  |  use add\_scatterer  |
|  rat       | \*\*double | I   | cartesian positions of atoms in path    |  use add\_scatterer  |
|  ipot      | \*long     | I   | unique potentials of atoms in path      |  use add\_scatterer  |
|  iorder    | long       | I   | order of approximation in genfmt        |  2                   |
|  nnnn      | bool       | I   | flag to write `feffNNNN.dat` file       |  false               |
|  xdi       | bool       | I   | flag to write `feffNNNN.xdi` file       |  false               |
|  verbose   | bool       | I   | flag to write screen messages           |  false               |
|  ipol      | bool       | I   | flag to do polarization calculation     |  false               |
|  evec      | \*double   | I   | polarization vector                     |  (0,0,0)             |
|  elpty     | double     | I   | ellipticity                             |  0                   |
|  xivec     | \*double   | I   | direction of X-ray propagation          |  (0,0,0)             |
|  ri        | \*double   | O   | leg lengths                             |                      |
|  beta      | \*double   | O   | beta angles                             |                      |
|  eta       | \*double   | O   | eta angles                              |                      |
|  reff      | double     | O   | half path length                        |  computed from ri    |
|  ne        | long       | O   | number of energy points actually used by Feff                | |
|  k         | \*double   | O   | k grid for feff path calculation, column 1 in `feffNNNN.dat` | |
|  real\_phc | \*double   | O   | central atom phase shifts. column 2 in `feffNNNN.dat`        | |
|  mag\_feff | \*double   | O   | magnitude of F\_eff, column 3 in `feffNNNN.dat`              | |
|  pha\_feff | \*double   | O   | phase of F\_eff, column 4 in `feffNNNN.dat`                  | |
|  red\_fact | \*double   | O   | reduction factor, column 5 in `feffNNNN.dat`                 | |
|  lam       | \*double   | O   | mean free path, column 6 in `feffNNNN.dat`                   | |
|  rep       | \*double   | O   | real part of complex momentum, column 7 in `feffNNNN.dat`    | |
|  errorcode | long       | O   | error code from `add_scatterer` or `make_path`               | |
|  errormessage | \*char  | O   | error message from `add_scatterer` or `make_path`            | |

Additionally, this struct captures several bits of information about
the potential model that Feff writes to the header of the `feffNNNN.dat` file:

| attribute  | type    | I/O | description                             |
| ---------- | ------- | --- |---------------------------------------- |
|  edge      | float   |  O  | energy threshold relative to atomic value (a poor estimate) |
|  gam\_ch   | float   |  O  | core level energy width |
|  kf        | float   |  O  | k value at Fermi level |
|  mu        | float   |  O  | Fermi level, eV |
|  rnorman   | float   |  O  | average Norman radius |
|  version   | string  |  O  | Feff version |
|  exch      | string  |  O  | brief description of the electronic exchange model |
|  rs\_int   | float   |  O  | interstitial radius |
|  vint      | float   |  O  | interstitial potential |
|  version   | string  |  O  | the version of feff and the feffpath revision |

A polarization calculation is enabled by setting the `ipol` element to
a true value.  `evec` has 3 elements and represents the polarization
vector of the incident beam.  `elpty` and `xivec` are the ellipticity
and Poynting vector of the incident beam for a calculation with
ellipticity.

When `make_path` is called, the `ri`, `beta`, and `eta` arrays are
filled to be `nleg` elements long and contain the geometry of the
scattering path.

Also when `make_path` is called, the arrays containing the columns of
a traditional `feffNNNN.dat` file are filled to be `ne` elements long.
These arrays are the same (besides precision) as the corresponding
columns.  While a `feffNNNN.dat` file can be exported (using the
`nnnn` flag), the need to write/read that file is obviated.

Direct access to `rat` and `ipot` is inconvenient and discouraged.
Use `add_scatterer`.

## A comment on the add_scatterer method

The way of building the geometry of the scattering path was
purposefully kept completely general.  While a program might use the
output of Feff's pathfinder, we also want to suport any other way of
generating path geometries.

To use the output of Feff's pathfinder, one would parse the
`paths.dat` file, making a path for each "paragraph" in the file and
calling `add_scatterer` for each atom in the path.

As a more complicated example, a reverse Monte Carlo approach to a fit
will move atoms to minimize difference between a model and the data.
As the atoms move about, the RMC program will generate new SS and MS
paths based on the current atom positions.  In that case, the current
coordinates can be used as the input to `add_scatterer` (with the
caveat that the `feffpath` library expects the absorber to be at the
origin, so the arguments to `add_scatterer` should have the absorber
position subtracted off).
