# Using the feffphases C wrapper

The program `makephases.c` will read from `libpotph.json` (this
example is made from a `feff.inp` file for copper metal) and write the
`phase.pad` file.

The `makephases.c` program likely will not compile before you build and
install feff85exafs.  The build script will look for the `feffphases`
library and the `feffphases.h` header file in their installation
locations

To compile the `makephases` sample programs:

	~> make

That's it!

## Sample program

Here is the simplest program using the C wrapper in `libfeffphases`:

```C
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
  strcpy(phases->phpad, "/path/to/phase.pad");
  ret = read_libpotph_json(phases);
  if (ret > 0) {
    printf("%s (error code %d)\n", phases->errormessage, phases->errorcode);
  } else {
	ret = make_phases(phases);
    if (ret > 0) {
      printf("%s (error code %d)\n", phases->errormessage, phases->errorcode);
    };
  }

  cleanup(phases);
  free(phases);
  return 0;
}

```

1. The `feffphases.h` header file is included.

2. A FEFFPHASES struct is created and is named `phases`.  Memory is
   allocated for it.

3. The call to `create_phases` allocates memory for all elements of the
   struct and initializes everything.

4. The `jsonfile` attribute of the struct is set to the location of
   the `libpotph.json` file, then it is read and imported into the
   struct.  (See note below.)

5. If no problems are found among the input data, the `phase.pad` file
   is calcuated and written to disk.

6. Finally, the memory for `phases` is deallocated and the program
   terminates.

## Contents of the FEFFPATH struct

| attribute    | type               | I/O | description                                         | feff.inp card |
| ------------ | ------------------ | --- |---------------------------------------------------- | ------------- |
| phpad        | \*char             | I   | path to output `phase.pad` file                     |               |
| jsonfile     | \*char             | I   | path to `libpotph.json`                             |               |
| verbose      | bool               | I   | flag to wrrite Feff screen messages, F=suppress     |               |
| errorcode    | int                | O   | integer error code                                  |               |
| errormessage | \*char             | O   | message associated with error code                  |               |
| ntitle       | int                | I   | number of header lines                              | TITLE         |
| titles       | \*\*char           | I   | (nheadx) array of header string                     | TITLE         |
| nat          | int                | I   | number of atoms in cluster                          | ATOMS         |
| rat          | \*\*double         | I   | (3,natx) cartesian coordinates of atoms in cluster  | ATOMS         |
| iphat        | \*int              | I   | (natx) unique potential indeces of atoms in cluster | ATOMS         |
| nph          | int                | I   | number of unique potentials                         | POTENTIALS    |
| iz           | \*int              | I   | (0:nphx) Z numbers of unique potentials             | POTENTIALS    |
| potlbl       | \*\*char           | I   | (0:nphx) labels of unique potentials                | POTENTIALS    |
| lmaxsc       | \*int              | I   | (0:nphx) l max for SCF for each potential           | POTENTIALS    |
| lmaxph       | \*int              | I   | (0:nphx) l max for FMS for each potential           | POTENTIALS    |
| xnatph       | \*double           | I   | (0:nphx) stoichiometry of each potential            | POTENTIALS    |
| spinph       | \*double           | I   | (0:nphx) spin on each unique potential              | POTENTIALS    |
| ihole        | int                | I   | edge index, 1=K, 4=L3, etc                          | HOLE/EDGE     |
| rscf         | float              | I   | cluster radius for self-consistent calculation      | SCF           |
| lscf         | int                | I   | 0=solid, 1=molecule                                 | SCF           |
| nscmt        | int                | I   | max number of self-consistency iterations           | SCF           |
| ca           | double             | I   | self-consistency convergence accelerator            | SCF           |
| nmix         | int                | I   | number of mixing iterations before Broyden          | SCF           |
| ecv          | double             | I   | core/valence separation energy                      | SCF           |
| icoul        | int                | I   | obsolete param. for handling Coulomb potential      | SCF           |
| ipol         | int                | I   | 1=do polarization calculation                       | POLARIZATION  |
| evec         | \*double           | I   | (3) polarization array                              | POLARIZATION  |
| elpty        | double             | I   | eccentricity of elliptical light                    | ELLIPTICITY   |
| xivec        | \*double           | I   | (3) ellipticity array                               | ELLIPTICITY   |
| ispin        | int                | I   | +/-2 = do spin calculation                          | SPIN          |
| spvec        | \*double           | I   | (3) spin array                                      | SPIN          |
| angks        | double             | I   | angle between spin and incidient beam               | SPIN          |
| ptz          | \*\*double complex | O   | (-1:1,-1:1) polarization tensor                     | return        |
| gamach       | double             | O   | tabulated core-hole lifetime                        | return        |
| ixc          | int                | I   | exchange index                                      | EXCHANGE      |
| vr0          | double             | I   | Fermi level offset                                  | EXCHANGE      |
| vi0          | double             | I   | constant broadening                                 | EXCHANGE      |
| ixc0         | int                | I   | exchange index for background function              | EXCHANGE      |
| iafolp       | int                | I   | 1=do automated overlapping                          | FOLP & AFOLP  |
| folp         | \*double           | I   | (0:nphx) overlapping fractions                      | FOLP & AFOLP  |
| xion         | \*double           | I   | (0:nphx) potential ionizations                      | ION           |
| rgrd         | double             | I   | radial grid used for the potentials/phases          | RGRID         |
| iunf         | int                | I   | 1=unfreeze f electrons                              | UNFREEZEF     |
| inters       | int                | I   |                                                     | INTERSTITIAL  |
| totvol       | double             | I   |                                                     | INTERSTITIAL  |
| jumprm       | int                | I   | 1=remove potential jumps at muffin tin radii        | JUMPRM        |
| nohole       | int                | I   | 1=compute without core-hole                         | NOHOLE        |
| iplsmn       | int                | I   | 1=compute with the MPSE approximation               | PLASMON       |


Note that `rscf` is a _single precision float_!

Note, where appropriate (`rscf`, `gamach`, `ecv`, `totvol`, `rat`),
these parameters are in *natural units* (Angstrom and eV).  They will
be converted to *code units* (Bohr and Hartree) by the library.

The error code and message will be 0 and an empty string when no
errors are found.

### A note about Feff screen messages

The default is to suppress Feff's many screen messages, ie. `verbose`
is false by default.  No effort is made to capture the path to
convergence when computing self-consistent potentials.  In order to
preserve that information, you must set `verbose` to true then somehow
capture and parse the screen output.  This is left as a chore for the
user interface.

### Attributes to consider trimming

1. `iafolp`: this may not be needed in the struct, the folp array gets filled
2. does the polarization and ellipticity info have any impact on the phases?
3. is spin to be supported in feff85exafs?
4. what does INTERSTITIAL do?  I've never quite understood that one....


## A note about the libpotph.json file

At this time (the date at the time of writing this note is 8 July
2015) the revamping of feff85exafs is a work in progress.  The
purpose is to facilitate the tight integration of feff into EXAFS data
analysis software.  Step one was to combine the roles of GENFMT and
FF2X into a stand-alone library, allowing user software to easily
generate the data table contained in a `feffNNNN.dat` file.  Step two
is to combine the roles of POT and XSPH into this stand-alone library,
alowing user software to easily generate the phase shifts needed for
path generation.

The full flowchart for an interaction with feff starts with gathering
data about the cluster and the details of the calculation.
Historically, this was done by reading a file called `feff.inp` using
the RDINP part of feff.  Armed with this information, the phases would
be calculated then the pathfinder would be used to enumerate the full
list of paths represented in the input cluster.  Finally,
`feffNNNN.dat` files would be generated.

The stand-alone libraries for phases and path are part of a plan to
disrupt this work flow.  Rather than relying on the quirky `feff.inp`
file, better user experiences could be created using modern UI and GUI
tools.  Until those new, better user experiences exist, we need a way
to get information into the phases library.

As a stop-gap measure, `rdinp` has been modified to write out a file
called `libpotph.json` which contains the content from a `feff.inp`
file required by the phases library in JSON format.  The
`read_libpotph_json` method used in the example above reads this JSON
file.

The eventual goal is that the UI somehow organizes the data needed by
the phases library, eliminating the need for either `feff.inp` or
`libpotph.json`.
