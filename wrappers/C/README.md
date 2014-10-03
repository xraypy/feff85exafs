# Using the feffpath C wrapper

The wrapper will be installed when feff85exafs is built.  This folder
contains an example of its use in a C program.  The program
`makepath.c` will read from `phase.bin` (calculated from copper metal)
and write the files `feff0001.dat` and `feff0004.dat`.

The `makepath.c` program likely will not compile before you build and
install feff85exafs.  The build script will look for various Feff
libraries and the `feffpath.h` header file in their installation
location.

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

  path->nnnn    = 1;
  path->verbose = 1;

  /* first path in copper */
  path->index   = 1;
  path->deg     = 12.0;
  ret = add_scatterer(path,  1.805, 0,  1.805, 1);
  ret = make_path(path);

  clear_path(path);

  free(path);
  return 0;
}

```

1. The `feffpath.h` header file is included.

2. A FEFFPATH struct is created and called `path`.  Memory is
   allocated for it.

3. The call to `create_path` allocates memory for all elements of the
   struct and initializes everything.

4. Setting `nnnn` and `verbose` to a true value tells the program to
   write `feffNNNN.dat` files and to write a short message to the
   screen when it is written

5. The path index is set to 1.  This means the output file will be
   called `feff0001.dat`.

6. The degeneracy of the path is specified.

7. The `add_scatterer` method is used to build the scattering
   geometry.  It's arguments are the Cartesian coordinates (referenced
   to the absorber at the origin) and the unique potential index of
   the atom.  In this case, a single scattering path is calculated, so
   `add_scatterer` is called only once.  For a multiple scattering
   path, `add_scatterer` would be called repeatedly, once for each
   atom (or leg) in the path.  `add_scatterer` fills the `rat` and
   `ipot` members of the struct and also sets `nleg`.

8. The call to `make_path` computes the parts of F-effective and
   stores them in the the `path` struct.

9. The call to `clear_path` reinitializes the struct.  This is not
   strictly necessary in this case, but would be were the program to
   go on to caompute another path.

10. Finally, the memory for `path` is deallocated and the program
    terminates.


## Contents of the FEFFPATH struct

Here are all the elements of the struct, their data types,
descriptions, and default values.  Input parameters are marked with
"I" and output parameters, set by the `make_path` method, are marked
with "O".

| element    | type       | I/O | description                             | default              |
| ---------- | --------   | --- |---------------------------------------- | -------------------- |
|  index     | long       | I   | path index                              |  9999                |
|  deg       | double     | I   | path degeneracy                         |  required input      |
|  nleg      | long       | I   | number of legs in path                  |  use add\_scatterer  |
|  rat       | \*\*double | I   | cartesian positions of atoms in path    |  use add\_scatterer  |
|  ipot      | \*long     | I   | unique potentials of atoms in path      |  use add\_scatterer  |
|  iorder    | long       | I   | order of approximation in genfmt        |  2                   |
|  nnnn      | bool       | I   | flag to write `feffNNNN.dat` file       |  false               |
|  json      | bool       | I   | flag to write `feffNNNN.json` file      |  false               |
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

A polarization calculation is enabled by setting the `ipol` element to
a true value.  `evec` contains a 3-vector with the polarization vector
of the incident beam.  `elpty` and `xivec` are the ellipticity and
Poynting vector of the incident beam for a calculation with
ellipticity.

When `make_path` is called, the `ri`, `beta`, and `eta` arrays are
filled to be `nleg` elements long and contain the geometry of the
scattering path.

Also when `make_path` is called, the arrays containing the columns of
a traditional `feffNNNN.dat` file are filled to be `ne` elements long.
These arrays are the same (besides precision) as the corresponding
columns.  While a `feffNNNN.dat` file can be exported (using the
`nnnn` flag), the need to write/read that file is obviated.


## A comment on the add_scatterer method

The way of building the geometry of the scattering path was
purposefully kept completely general.  While a program might use the
output of Feff's pathfinder, we also want to suport other ways of
generating path geometries.

To use the output of Feff's pathfinder, one would parse the
`paths.dat` file, calling `add_scatterer` for each "paragraph" in the
file.

As a more complicated example, a reverse Monte Carlo approach to a fit
will move atoms to minimize difference between a model and the data.
As the atoms move about, the RMC program will generate new SS and MS
paths based on the current atom positions.  In that case, the current
coordinates can be used as the input to `add_scatterer` (with the
caveat that the `feffpath` library expects the absorber to be at the
origin, so the arguments to `add_scatterer` should have the absorber
position subtracted off).
