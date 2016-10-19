# Using the onepath Fortran entry point

The `libpotph` library will be installed when feff85exafs is built.
This folder contains an example of its use in a Fortran program.  The
program `makepotph` will read from `libpotph.json` (this example is
generated for copper metal) and write the file `phase.pad`.

The `makepotph.f` program likely will not compile before you build and
install feff85exafs.  The build script will look for the `libpotph`
library.

To compile the `makepotph` sample program:

	~> make

That's it!

## Sample program

Here is the simplest program using the Fortran entry point in
`libpotph`:

```fortran
  	  program makepotph
	  implicit double precision (a-h, o-z)

      integer nphx, natx, nphx
      parameter (nphx = 11, natx = 1000, nheadx = 30)

      character*80 title(nheadx)
      integer ntitle, nat, nph, iphat(natx), ipol, ispin, ihole
      integer lfms1, nmix, icoul, ixc, ixc0, iafolp, iunf, inters
      integer jumprm, nohole, iplsmn
      integer iz(0:nphx), lmaxsc(0:nphx), lmaxph(0:nphx)
      double precision evec(3), xivec(3), spvec(3), spinph(0:nphx)
      complex*16 ptz(-1:1, -1:1)
      double precision rat(3,natx)
      double precision xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
      character*6 potlbl(0:nphx)
      real rfms1
      double precision elpty, angks, gamach, ca1, ecv
      double precision vr0, vi0, rgrd, totvol

	  character*256 phpad
	  logical verbse

	  call inipotph(ntitle, title, nat, rat, iphat,
     1       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
     2       ihole, rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     3       ipol, evec, elpty, xivec, ispin, spvec, angks,
     4       ptz, gamach, ixc, vr0, vi0, ixc0,
     5       iafolp, folp, xion, rgrd, iunf,
     6       inters, totvol, jumprm, nohole, iplsmn)

	  call json_read_libpotph(ntitle, title, nat, rat, iphat,
     1       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
     2       ihole, rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     3       ipol, evec, elpty, xivec, ispin, spvec, angks,
     4       ptz, gamach, ixc, vr0, vi0, ixc0,
     5       iafolp, folp, xion, rgrd, iunf,
     6       inters, totvol, jumprm, nohole, iplsmn)

	  phpad = 'phase.pad'
	  verbse = .true.
	  call libpotph(phpad, verbse,
	 1       ntitle, title, nat, rat, iphat,
     2       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
     3       ihole, rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     4       ipol, evec, elpty, xivec, ispin, spvec, angks,
     5       ptz, gamach, ixc, vr0, vi0, ixc0,
     6       iafolp, folp, xion, rgrd, iunf,
     7       inters, totvol, jumprm, nohole, iplsmn)

	  stop
	  end
```

 1. All the necessary variables are typed, dimensioned, and initialized.

 2. Everything gets initialized.  Note that `rfms1` is single precision,
    not double.

 3. Read the values of parameters from `feff.inp`, which were written
    to a file called `libpotph.json` by RDINP.  (See note below.)

 4. `json_read_libpotph` converts to code units (bohr and hartree).
    `libpotph` expects code units, **not** natural units.  If you use
    some other mechanism to get data into `libpotph`, you **must**
    remember to convert `rfms1`, `gamach`, `ecv`, `totvol`, and all of
    `rat` to code units.  See, for example, lines 342 - 350 in
    `src/JSON/json_read_libpotph`.

 5. Specify the path and name of the output `phase.pad` file.

 6. Call the libpotph library, which sorts the input cluster, computes
    the muffin tin potentials, and writes phase shifts to a file
    as specified by the phpad argument.

## Arguments of the libpotph subroutine

Here are all the arguments to `libpotph`, their data types,
descriptions, and which `feff.inp` card they are associated with.
Input parameters are marked with "I" and output parameters, set by the
`make_path` method, are marked with "O".  Most of the names are chosen
to be consistent with the naming conventions in Feff.

| element    | type                 | I/O | description                                                  | Feff card    |
| ---------- | -------------------- | --- |------------------------------------------------------------- | ------------ |
|  phpad     | character\*256       | O   | path and name of output phase.pad file                       |              |
|  verbse    | logical              | I   | flag to write screen messages, .false.=suppress              |              |
|  ntitle    | integer              | I   | number of title lines                                        | TITLE        |
|  title     | character\*80        | I   | array(nheadx) of title lines                                 | TITLE        |
|  nat       | integer              | I   | number of atoms in cluster                                   | ATOMS        |
|  rat       | double precision     | I   | array(3,natx) of Cartesian coordinates                       | ATOMS        |
|  iphat     | integer              | I   | array(natx) of potential indeces                             | ATOMS        |
|  nph       | integer              | I   | number of potential indeces                                  | POTENTIALS   |
|  iz        | integer              | I   | array(0:nphx) of unique potential Z numbers                  | POTENTIALS   |
|  potlbl    | character\*6         | I   | array(0:nphx) of unique potential labels                     | POTENTIALS   |
|  lmaxsc    | integer              | I   | array(0:nphx) of unique potential l max for self consistency | POTENTIALS   |
|  lmaxph    | integer              | I   | array(0:nphx) of unique potential l max for FMS              | POTENTIALS   |
|  xnatph    | double precision     | I   | array(0:nphx) of unique potential stoichiometries            | POTENTIALS   |
|  spinph    | double precision     | I   | array(0:nphx) of unique potential spin values                | POTENTIALS   |
|  ihole     | integer              | I   | edge index (K=1, L:3=4, etc)                                 | HOLE/EDGE    |
|  rfms1     | double precision     | I   | radial size of cluster for self consistency                  | SCF          |
|  lfms1     | integer              | I   | for self consistency, 0=solid, 1=molecule                    | SCF          |
|  nscmt     | integer              | I   | max number of self-consistency iterations                    | SCF          |
|  ca1       | double precision     | I   | self-consistency convergence accelerator                     | SCF          |
|  nmix      | integer              | I   | number of mixing iterations before Broyden                   | SCF          |
|  ecv       | double precision     | I   | core valence separation energy                               | SCF          |
|  icoul     | integer              | I   | obsolete param. for handling Coulomb pot.                    | SCF          |
|  ipol      | integer              | I   | 1=do polarization calculation                                | POLARIZATION |
|  evec      | double precision     | I   | array(3) polarization vector                                 | POLARIZATION |
|  elpty     | double precision     | I   | eccentricity of elleptical light                             | ELLIPTICITY  |
|  xivec     | double precision     | I   | array(3) ellipticity vector                                  | ELLIPTICITY  |
|  ispin     | integer              | I   | 1=do spin calculation                                        | SPIN         |
|  spvec     | double precision     | I   | array(3) spin vector                                         | SPIN         |
|  angks     | double precision     | I   | angle between spin vector and incident beam                  | SPIN         |
|  ptz       | double precision     | O   | array(-1:1, -1:1) polarization tensor                        |              |
|  gamach    | double precision     | O   | tabulated core-hole lifetime                                 |              |
|  ixc       | integer              | I   | echange index                                                | EXCHANGE     |
|  vr0       | double precision     | I   | Fermi level shift                                            | EXCHANGE     |
|  vi0       | double precision     | I   | "optical" loss term                                          | EXCHANGE     |
|  ixc0      | integer              | I   | exchange index used for background function                  | EXCHANGE     |
|  iafolp    | integer              | I   | 1=do automated potential overlap                             | AFOLP        |
|  folp      | double precision     | I   | array(0:nphx) overlap fractions                              | FOLP/AFOLP   |
|  xion      | double precision     | I   | array(0:nphx) potential ionizations                          | ION          |
|  rgrd      | double precision     | I   | radial grid used for the potentials/phases                   | RGRID        |
|  iunf      | integer              | I   | 1=unfreeze f electrons                                       | UNFREEZEF    |
|  inters    | integer              | I   | interstitial potential and density param                     | INTERSTITIAL |
|  totvol    | double precision     | I   | interstitial potential and density param                     | INTERSTITIAL |
|  jumprm    | integer              | I   | 1=remove potential jumps at muffin tin radii                 | JUMPRM       |
|  nohole    | integer              | I   | 1=compute without core-hole                                  | NOHOLE       |
|  iplsmn    | integer              | I   | 1=compute with the MPSE approximation                        | PLASMON      |


### A note about Feff's screen messages

Feff traditionally spits a lot of text to the screen, particularly
when computing self-consistent potentials.  This screen output can be
suppressed by setting `verbse` to `.false.`, making the potentials and
phases calculation run completely silently (except for terminal
errors).  Setting `verbse` to `.true.` makes the library behave much
like the stand-alone `pot` and `xsph` programs.


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
