
# Content of the OPCONSAT folder

This directory contains routines for themany-pole self-energy model
using the optical constants: atomic approximation model.  See
[the 2007 Kas, et al. paper](http://dx.doi.org/10.1103/PhysRevB.76.195116).

The original Fortran 77 and Fortran 90 routines in this directory are
covered by the [LICENSE](../HEADERS/license.h).

This directory also contains the `libfeffloss` library, which combines
the functionality of the `opconsat` and `eps2exc` programs into a
single, callable library.  This library is, in turn, called by the
[`libfeffphases` library](../POT/README.md) whenever the MPSE
calculation is turned on.

# Build and install

To build, type `scons`.  This will build:

 * `opconsat`: the first stand-alone program
 * `eps2exc`: the second stand-alone program
 * `libfeffloss.so`: the Fortran entry point for computing the MPSE

Once built, type `scons install` to install everything:

 * `libfeffloss.so`: installed to `/usr/local/lib`
 * `opconsat`, `eps2exc`: installed to `/usr/local/bin`

You **must** install before building the Perl or Python wrappers.
Other wrappers almost certainly require at least that `libfeffloss.so`
be installed.

# Sample program

Assuming that a
[`libpotph.json`](../../wrappers/fortran/libpotph.json) file is
present in the current directory, this program will read the
`libpotph.json` information and compute the `loss.dat` file.

```fortran
      program loss

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'

      integer epsmax
      parameter(epsmax = 700)
      
c      integer iz(0:nphx)
      integer npoles
      double precision rnrm(0:nphx), eps0, gamma
      logical write_loss, write_opcons, write_exc, verbose

      double precision wpcorr(MxPole), delta(MxPole), ampfac(MxPole)

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


      
      write_loss   = .true.
      write_opcons = .false.
      write_exc    = .false.
      verbose      = .true.
      npoles       = 100
      eps0         = -1.d0



      call inipotph(
c     TITLE
     1       ntitle, title,
c     ATOMS
     2       nat, rat, iphat,
c     POTENTIALS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
c     HOLE/EDGE
     4       ihole,
c     SCF
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c     POLARIZATION, ELLIPTICITY
     6       ipol, evec, elpty, xivec,
c     SPIN
     7       ispin, spvec, angks,
c     computed
     8       ptz, gamach,
c     EXCHANGE
     9       ixc, vr0, vi0, ixc0,
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
     _       iafolp, folp, xion, rgrd, iunf,
c     INTERSTITIAL, JUMPRM, NOHOLE
     1       inters, totvol, jumprm, nohole, iplsmn)


c*****************************************************************************
c     read the contents of a json file that includes all of global.json,
c     atoms.json, pot.json & xpsh.json (i.e. global.dat, atoms.dat, mod1.inp,
c     and mod2.inp)
c*****************************************************************************
      call json_read_libpotph(
c     TITLE
     1       ntitle, title,
c     ATOMS
     2       nat, rat, iphat,
c     POTENTIALS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
c     HOLE/EDGE
     4       ihole,
c     SCF
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c     POLARIZATION, ELLIPTICITY
     6       ipol, evec, elpty, xivec,
c     SPIN
     7       ispin, spvec, angks,
c     computed
     8       ptz, gamach,
c     EXCHANGE
     9       ixc, vr0, vi0, ixc0,
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
     _       iafolp, folp, xion, rgrd, iunf,
c     INTERSTITIAL, JUMPRM, NOHOLE
     1       inters, totvol, jumprm, nohole, iplsmn)

c  this is a cheat .. rather than calling rdpot and all that, just
c  hard wire the values for a Copper calculation
      rnrm(0) = 2.8384890981392523
      rnrm(1) = 2.6294479894989911

      call feffloss(nph, iz, xnatph, rnrm, npoles, eps0,
     1       write_opcons, write_loss, write_exc, verbose,
     2       wpcorr, gamma, ampfac, delta)
      

      end

```


## Arguments to the feffloss subroutine



| element        | type                 | I/O | description                                                  | Feff card    |
| -------------- | -------------------- | --- |------------------------------------------------------------- | ------------ |
|  nph           | integer              | I   | number of potential indeces                                  | POTENTIALS   |
|  iz            | integer              | I   | array(0:nphx) of unique potential Z numbers                  | POTENTIALS   |
|  xnatph        | double precision     | I   | array(0:nphx) of unique potential stoichiometries            | POTENTIALS   |
|  rnrm          | double precision     | I   | array(0:nphx) of norman radii                                | from rdpot   |
|  npoles        | integer              | I   | number of poles to use in the calculation                    |              |
|  eps0          | double precision     | I   | dielectric constant                                          |              |
|  write\_opcons | logical              | I   | write opconsEE.dat file                                      |              |
|  write\_loss   | logical              | I   | write loss.dat file                                          |              |
|  write\_exc    | logical              | I   | write exc.dat file                                           |              |
|  verbose       | logical              | I   | write screen messages                                        |              |
|  wpcorr        | double precision     | O   | array(MxPole) energy grid for MPSE calculation               |              |
|  gamma         | double precision     | O   | array(MxPole) .... for MPSE calculation                      |              |
|  ampfac        | double precision     | O   | array(MxPole) .... for MPSE calculation                      |              |
|  delta         | double precision     | O   | array(MxPole) .... for MPSE calculation                      |              |

