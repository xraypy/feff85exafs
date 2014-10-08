# Using the onepath fortran entry point

The onepath library will be installed when feff85exafs is built.  This
folder contains an example of its use in a fortran program.  The
program `makepath` will read from `phase.bin` (this example is
calculated from copper metal) and write the files `feff0001.dat` and
`feff0004.dat`.

The `makepath.f` program likely will not compile before you build and
install feff85exafs.  The build script will look for the `onepath`
libraries,

To compile the `makepath` sample program:

	~> scons

That's it!

## Sample program

Here is the simplest program using the C wrapper:

```fortran
  	  program makepath
	  implicit double precision (a-h, o-z)

c     taken from the feff HEADERS/dim.h
      integer nex, legtot
      parameter (nex = 150, legtot = 9)
      double precision evec(3), xivec(3)

	  integer  iorder
      double precision elpty
      integer innnn, ijson, ivrbse

      double precision rat(3,0:legtot+1)
      integer ipot(0:legtot)

      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)

c     initialize everything
      index   = 9999
      nleg    = 0
      deg     = 0
      iorder  = 2
      innnn   = 1 
      ijson   = 0 
      ivrbse  = 1
      ipol    = 0 
      elpty   = 0
      ne      = 0
      
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
      index    = 1
      nleg     = 2
      deg      = 12
      call addatom(1, -1.805, 0., -1.805, 1, ipot, rat)
      call onepath(index, nleg, deg, iorder,
     &       ipot, rat,
     &       ipol, evec, elpty, xivec,
     &       innnn, ijson, ivrbse, ri, beta, eta,
     &       ne,col1,col2,col3,col4,col5,col6,col7)

      end

```

1. All the necessary variables are typed, dimensioned, and initialized.

2. Setting `nnnn` and `verbose` to a true value tells the program to
   write `feffNNNN.dat` files and to write a short message to the
   screen when it is written.

3. The path index is set to 1.  This means the output file will be
   called `feff0001.dat`.

4. The degeneracy of the path is specified.

5. The `addatom` subroutine is used to build the scattering geometry.
   It's arguments are leg index, the Cartesian coordinates (referenced
   to the absorber *at the origin*) and the unique potential index of
   the atom.  In this case, a single scattering path is calculated, so
   `addatom` is called only once.  For a multiple scattering
   path, `addatom` would be called repeatedly, once for each
   atom (or leg) in the path.  `addatom` fills the `rat` and
   `ipot` arrays.

6. The call to `onepath` computes the parts of F-effective and
   stores them in the `colN` arrays.


## Arguments of the onepath subroutine

Here are all the elements of the struct, their data types,
descriptions, and default values.  Input parameters are marked with
"I" and output parameters, set by the `make_path` method, are marked
with "O".  Most of the names are chosen to be consistent with the
naming conventions in Feff.  The output arrays for the columns of
`feffNNNN.dat` are those
[used by Larch](http://xraypy.github.io/xraylarch/xafs/feffpaths.html#the-feffdat-group-full-details-of-the-feff-dat-file).

| element    | type        | I/O | description                             | default              |
| ---------- | ----------  | --- |---------------------------------------- | -------------------- |
|  index     | integer     | I   | path index                              |  9999                |
|  deg       | double      | I   | path degeneracy                         |  required input      |
|  nleg      | integer     | I   | number of legs in path                  |  required input      |
|  rat       | double(3,0:legtot+1 | I   | cartesian positions of atoms in path |  use addatom    |
|  ipot      | integer(legtot)  | I   | unique potentials of atoms in path      |  use addatom    |
|  iorder    | integer     | I   | order of approximation in genfmt        |  2                   |
|  innnn     | integer     | I   | flag to write `feffNNNN.dat` file       |  0                   |
|  ijson     | integer     | I   | flag to write `feffNNNN.json` file      |  0                   |
|  iverbose  | integer     | I   | flag to write screen messages           |  0                   |
|  ipol      | integer     | I   | flag to do polarization calculation     |  0                   |
|  evec      | double(3)   | I   | polarization vector                     |  (0,0,0)             |
|  elpty     | double      | I   | ellipticity                             |  0                   |
|  xivec     | double(3)   | I   | direction of X-ray propagation          |  (0,0,0)             |
|  ri        | double(legtot)     | O | leg lengths                        |                      |
|  beta      | double(legtot+1)   | O | beta angles                        |                      |
|  eta       | double(0:legtot+1) | O | eta angles                         |                      |
|  ne        | integer     | O   | number of energy points actually used by Feff                | |
|  k         | double(nex) | O   | k grid for feff path calculation, column 1 in `feffNNNN.dat` | |
|  real\_phc | double(nex) | O   | central atom phase shifts. column 2 in `feffNNNN.dat`        | |
|  mag\_feff | double(nex) | O   | magnitude of F\_eff, column 3 in `feffNNNN.dat`              | |
|  pha\_feff | double(nex) | O   | phase of F\_eff, column 4 in `feffNNNN.dat`                  | |
|  red\_fact | double(nex) | O   | reduction factor, column 5 in `feffNNNN.dat`                 | |
|  lam       | double(nex) | O   | mean free path, column 6 in `feffNNNN.dat`                   | |
|  rep       | double(nex) | O   | real part of complex momentum, column 7 in `feffNNNN.dat`    | |

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

Direct access to `rat` and `ipot` is inconvenient and discouraged.
Use `add_scatterer`.

## A comment on the addatom subroutine

This is a helper subroutine to help organize the contents of `ipot`
and `rat` and to convert the Cartesian coordinates to code units.

```fortran
      subroutine addatom(leg, x, y, z, ip, ipot, rat)
      implicit double precision (a-h, o-z)
      integer npatx
      parameter (npatx = 8)
      integer legtot
      parameter (legtot=npatx+1)

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
