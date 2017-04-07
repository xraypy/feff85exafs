# Using the onepath Fortran entry point

The `onepath` library will be installed when feff85exafs is built.
This folder contains an example of its use in a Fortran program.  The
program `makepath` will read from `phase.pad` (this example is
calculated from copper metal) and write the files `feff0001.dat` and
`feff0004.dat`.

The `makepath.f` program likely will not compile before you build and
install feff85exafs.  The build script will look for the `onepath`
library.

To compile the `makepath` sample program:

	~> make

That's it!

## Sample program

Here is the simplest program using the Fortran entry point:

```fortran
  	  program makepath
	  implicit double precision (a-h, o-z)

c     taken from the feff HEADERS/dim.h
      integer nex, legtot, ixc
      parameter (nex = 150, legtot = 9)
      double precision evec(3), xivec(3)

	  integer iorder, innnn, ixdi, ivrbse
	  integer ipot(0:legtot), iz(0:nphx)
      double precision elpty, rat(3,0:legtot+1)

      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)

      character*256 phpad
      character*8 cxc

c     initialize everything
      phpad   = 'phase.pad'
      index   = 9999
      nleg    = 0
      deg     = 1.0
      iorder  = 2
      innnn   = 1 
      ixdi    = 1
      ivrbse  = 1
      ipol    = 0 
      elpty   = 0.0
      ne      = 0

      cxc    = ''
      versn  = ''
      rs     = 0.
      vint   = 0.
      xmu    = 0.
      edge   = 0.
      xkf    = 0.
      rnrmav = 0.
      gamach = 0.
      
      do 5  i=1,3
         evec(i)  = 0
         xivec(i) = 0
 5    continue
      do 10 i=0, legtot
         rat(1,i) = 0
         rat(2,i) = 0
         rat(3,i) = 0
         ipot(i)  = 0
         eta(i)   = 0
         if (i>0) then
            ri(i)   = 0
            beta(i) = 0
            eta(i)  = 0
         endif
 10   continue
      beta(legtot+1) = 0
      eta(legtot+1)  = 0

      do 15 i=0,nphx
         iz(i) = 0
 15   continue

	  do 20 i=1, nex
         col1(i) = 0
         col2(i) = 0
         col3(i) = 0
         col4(i) = 0
         col5(i) = 0
         col6(i) = 0
         col7(i) = 0
 20   continue
c     done initializing

c     compute first shell of Copper (SS, deg=12)
      index  = 1
      nleg   = 2
      deg    = 12

c                  leg   x	  y      z   ip  ipot and rat arrays
      call addatom(1, -1.805, 0., -1.805, 1, ipot, rat)
      call onepath(phpad, index, nleg, deg, iorder,
     &     cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &     versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &     innnn, ixdi, ivrbse, ri, beta, eta,
     &     ne,col1,col2,col3,col4,col5,col6,col7)

c
c     now you can do things with the colN arrays
c
      end

```

1. All the necessary variables are typed, dimensioned, and initialized.

2. `phpad` is a character*256 parameter containing the path to the
   `phase.pad` file.  If not specified (or if the file does not
   exist), `onepath` will try to read `phase.pad` in the current
   working directory.

3. Setting `nnnn` and `verbose` to a true value tells the program to
   write `feffNNNN.dat` files and to write a short message to the
   screen when it is written.  (Set `xdi` to true to get this file in
   XDI format.

4. The path index is set to 1.  This means the output file will be
   called `feff0001.dat`.

5. The degeneracy of the path is specified.

6. The `addatom` subroutine (scroll down for an explanation) is used
   to build the scattering geometry.  It's arguments are the leg
   index, the Cartesian coordinates (referenced to the absorber *at
   the origin*) and the unique potential index of the atom.  In this
   case, a single scattering path is calculated, so `addatom` is
   called only once.  For a multiple scattering path, `addatom` would
   be called repeatedly, once for each atom (or leg) in the path.
   `addatom` fills the `rat` and `ipot` arrays.

7. The call to `onepath` computes the parts of F-effective and
   stores them in the `colN` arrays.


## Arguments of the onepath subroutine

Here are all the elements of the struct, their data types,
descriptions, and default values.  Input parameters are marked with
"I" and output parameters, set by the `make_path` method, are marked
with "O".  Many of the names are chosen to be consistent with the
naming conventions in Feff.

| element    | type          | I/O | description                             | default              |
| ---------- | ------------  | --- |---------------------------------------- | -------------------- |
|  phpad     | character*256 | I   | path to `phase.pad` file                |  `phase.pad`         |
|  index     | integer       | I   | path index                              |  9999                |
|  deg       | double        | I   | path degeneracy                         |  required input      |
|  nleg      | integer       | I   | number of legs in path                  |  required input      |
|  rat       | double(3,0:legtot+1) | I   | cartesian positions of atoms in path |  use addatom     |
|  ipot      | integer(legtot)      | I   | unique potentials of atoms in path   |  use addatom     |
|  iorder    | integer       | I   | order of approximation in genfmt        |  2                   |
|  innnn     | integer       | I   | flag to write `feffNNNN.dat` file       |  0                   |
|  ixdi      | integer       | I   | flag to write `feffNNNN.xdi` file       |  0                   |
|  iverbose  | integer       | I   | flag to write screen messages           |  0                   |
|  ipol      | integer       | I   | flag to do polarization calculation     |  0                   |
|  evec      | double(3)     | I   | polarization vector                     |  (0,0,0)             |
|  elpty     | double        | I   | ellipticity                             |  0                   |
|  xivec     | double(3)     | I   | direction of X-ray propagation          |  (0,0,0)             |
|  ri        | double(legtot)     | O | leg lengths                          |                      |
|  beta      | double(legtot+1)   | O | beta angles                          |                      |
|  eta       | double(0:legtot+1) | O | eta angles                           |                      |
|  ne        | integer       | O   | number of energy points actually used by Feff                | |
|  col1      | double(nex)   | O   | k grid for calculation, column 1 in `feffNNNN.dat`           | |
|  col2      | double(nex)   | O   | central atom phase shifts. column 2 in `feffNNNN.dat`        | |
|  col3      | double(nex)   | O   | magnitude of F\_eff, column 3 in `feffNNNN.dat`              | |
|  col4      | double(nex)   | O   | phase of F\_eff, column 4 in `feffNNNN.dat`                  | |
|  col5      | double(nex)   | O   | reduction factor, column 5 in `feffNNNN.dat`                 | |
|  col6      | double(nex)   | O   | mean free path, column 6 in `feffNNNN.dat`                   | |
|  col7      | double(nex)   | O   | real part of complex momentum, column 7 in `feffNNNN.dat`    | |

Additionally, the entry point returns several bits of information
about the potential model that Feff traditionally writes to the header
of the `feffNNNN.dat` file:

| attribute  | type    | I/O | description                             |
| ---------- | ------- | --- |---------------------------------------- |
|  versn     | string  |  O  | the version of feff and the feffpath revision |
|  cxc       | string  |  O  | brief description of the electronic exchange model |
|  edge      | float   |  O  | energy threshold relative to atomic value (a poor estimate) |
|  gamach    | float   |  O  | core level energy width |
|  xkf       | float   |  O  | k value at Fermi level |
|  xmu       | float   |  O  | Fermi level, eV |
|  rnrmav    | float   |  O  | average Norman radius |
|  rs        | float   |  O  | interstitial radius |
|  vint      | float   |  O  | interstitial potential |

A polarization calculation is enabled by setting the `ipol` element to
a true value.  `evec` has 3 elements and represents the polarization
vector of the incident beam.  `elpty` and `xivec` are the ellipticity
and Poynting vector of the incident beam for a calculation with
ellipticity.

When `onepath` is called, the `ri`, `beta`, and `eta` arrays are
filled to be `nleg` elements long and contain the geometry of the
scattering path.

Also when `onepath` is called, the arrays containing the columns of
a traditional `feffNNNN.dat` file are filled to be `ne` elements long.
These arrays are the same (besides precision) as the corresponding
columns.  While a `feffNNNN.dat` file can be exported (using the
`nnnn` flag), the need to write/read that file is obviated.

While you can certainly fill `rat` and `ipot` directly, `addatom` is
more convenient.  In any case, be aware that the units in `rat` are
bohr, not Angstrom.

## A comment on the addatom subroutine

This is a helper subroutine to help organize the contents of `ipot`
and `rat` and to convert the Cartesian coordinates to code units.

```fortran
      subroutine addatom(leg, x, y, z, ip, ipot, rat)
      implicit double precision (a-h, o-z)
      integer legtot
      parameter (legtot=9)

      integer leg, ip
      real x, y, z
      double precision rat(3,0:legtot+1)
      integer ipot(0:legtot)
      double precision bohr
      parameter (bohr = 0.529 177 249d0)

      ipot(leg)  = ip
      rat(1,leg) = dble(x) / bohr
      rat(2,leg) = dble(y) / bohr
      rat(3,leg) = dble(z) / bohr
      return
      end
```

The way of building the geometry of the scattering path was
purposefully kept completely general.  While a program might use the
output of Feff's pathfinder, we also want to suport any other way of
generating path geometries.

To use the output of Feff's pathfinder, one would parse the
`paths.dat` file, making a path for each "paragraph" in the file and
calling `addatom` for each atom in the path.

As a more complicated example, a reverse Monte Carlo approach to a fit
will move atoms to minimize difference between a model and the data.
As the atoms move about, the RMC program will generate new SS and MS
paths based on the current atom positions.  In that case, the current
coordinates can be used as the input to `addatom` (with the
caveat that the `onepath` library expects the absorber to be at the
origin, so the arguments to `addatom` should have the absorber
position subtracted off).
