# Building feff85exafs

In an attempt to move beyond the inadequate collection of Makefiles
that this project inherited from the Feff Project, Bruce decided to
use the [scons software construction tool](http://www.scons.org/) as a
replacement.  scons seemed easier to understand than GNU Autotools
(which Matt used ages ago for feff6l and Ifeffit).  Because scons is
written in python and uses configuration files that contain python
code, it seemed like a good fit for something that would eventually be
integrated into [Larch](https://github.com/xraypy/xraylarch).

Compilation is as simple as typing `scons` after cd-ing into the
`src/` folder.

You can also cd into one of the subfolders and do `scons` to build
just that piece.

The file `FeffBuild.py` has a bit of logic for specifying compilation
flags, which may be attractive for optimization or other reasons, and
installation locations.  This is then imported by each of the
lower-level SConstruct files.  Except for gfortran, the initial
guesses for compilation flags were taken from the top level Makefile
that came from the Feff Project.  The gfortran flags were selected by
Bruce after playing around on his Ubuntu system.

---

# Installed files

Many of the libraires and stand-alone programs have been renamed from
the names chosen by the original build system for feff85exafs.

## Feff Libraries

The following libraries contain the various parts of Feff.

* `ATOM/libfeffatom.so`
* `COMMON/libfeffcom.so`
* `DEBYE/libfeffdw.so`
* `EXCH/libfeffexch.so`
* `FMS/libfefffms.so`
* `FOVRG/libfeffpha.so`
* `GENFMT/libfeffgenfmt.so`
* `JSON/libfeffjson.so`
* `MATH/libfeffmath.so`
* `PAR/libfeffpar.so`
* `POT/libfeffint.so`

Default installation locations:

* Linux: `/usr/local/lib`
* Windows: `C:\path\to\lib`  :FIXME:
* Mac: `/some/where/lib`     :FIXME:

This can be set from the command line:

	~> scons prefix="/other/location"

where the default value for "prefix" is `/usr/local` on Linux, etc.

This location **must** be in the linker/loader path.  With bash, for
example, you may need to do

	~> export LD_LIBRARY_PATH=/usr/local/lib

or, if LD\_LIBRARY\_PATH already has a value:

	~> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH ":/usr/local/lib"

You may want to put that in your `.bashrc` file.

## Stand-alone programs

* `RDINP/rdinp`: input file reader 
* `POT/pot`: module 1, potentials calculation
* `XSPH/xsph`: module 2, phase shifts calculation
* `PATH/pathfinder`: module 4, path finder
* `GENFMT/genfmt`: module 5, F-matrix calculation
* `FF2X/ff2x`: module 6, output files

Note that module 3, `fms`, the full multiple scattering XANES
calculator, is not a part of feff85exafs.

Default installation locations:

* Linux: `/usr/local/bin`
* Windows: `C:\path\to\bin`  :FIXME:
* Mac: `/some/where/bin`     :FIXME:

This can be set from the command line:

	~> scons prefix="/other/location"

where the default value for "prefix" is `/usr/local` on Linux, etc.

Note that the ultimate goal of the feff85exafs project is to do away
with the stand-alone programs.
 * `rdinp` is a chore better handled by a GUI or other user interface.
 * `genfmt` and `ff2x` are replaced by the feffpath library, which can
   be called directly by a program written in fortran, C, or some
   other language
 * eventually `pot` and `xsph` will be replaced by similarly callable
   library
 * finally the pathfinder is missing critical features (most
   prominently: caching geometry of degenerate paths and fuzzy
   degeneracy).  the pathfinder has alrady been rewritten in Perl for
   Demeter, for example.


## The feffpath library and its wrappers

These files are the entry points for various programming languages to
the "feffpath" calculation, that is, the calculation of a single
`feffNNNN.dat` file given an input scattering geometry.

This presumes that `pot` and `xsph` have already been run and that the
`phase.bin` file is accessible.

* `GENFMT/libonepath.so`: This is the Fortran entry point.
* `GENFMT/libfeffpath.so`: This is the C wrapper around the Fortran onepath
* `GENFMT/feffpath.h`: This is the header file, almost certainly required by any language wrapper
* `GENFMT/perl/feffpath_wrap.c`: The SWIG generated wrapper file for use with the Python wrapper
* `GENFMT/perl/FeffPathWrapper.pm`: This is the SWIG generated Perl wrapper
* `GENFMT/python/feffpath_wrap.c`: The SWIG generated wrapper file for use with the Perl wrapper
* `GENFMT/python/feffpathwrapper.py`: This is the SWIG generated Perl wrapper

`libonepath.so` and `libfeffpath.so` will be installed to the same
location (`/usr/local/lib`, etc) as the Feff libraries.

The other files are transferred into the proper place in the
`wrappers` folder of the feff85exafs distribution so that the
language-specific wrapper can be built using language-specific build
tools.


---

# Compiling against json-fortran

At the moment, feff85exafs has been modified to replace all
intermediate i/o files (`modN.inp`, `global.dat`, `geom.dat`, and so
on) with JSON file.  Reading and writing the JSON files is currently
done using
[json-fortran](https://github.com/jacobwilliams/json-fortran).

According to the json-fortran web page, it can be compiled with
Visual Studio 2010, Intel Fortran, and gfortran 4.9.  My only
experience is with gfortran 4.9.  I know that gfortran 4.8 will not
do.

Once I installed gfortran 4.9
([this helped on my Ubuntu machine](http://askubuntu.com/questions/428198/getting-installing-gcc-g-4-9-on-ubuntu)),
I was able to use the `build.sh` script that comes with json-fortran
to build the static library and the test program.  The static library
and corresponding module end up in the `lib/` folder.

The files in `lib/` need to end up someplace where the compiler
and linker can find them.  The simple build script doesn't do an
install.  I put the `libjsonfortran.a` file in `/usr/lib` and the
`json_module.mod` file in `src/JSON/`.  This way, the linker was able
to find the library and the `FortranCompilation.py` file, which
configures compilation flags for the scons build system, sets the `-I`
flag to find the module file.

Obviously this is a bit crufty, but it works well enough for now.

## Executive summary

  1. Build json-fortran by running the `build.sh` script
  2. Copy `libjsonfortran.a` to `/usr/lib`
  3. Copy `json_module.mod` to `src/JSON/`
  4. Run scons to build feff against json-fortran

## Matt's comment on json-fortran


Matt has these reasonable things to say about compiling against json-fortran:

    Do you know what versions of Fortran the json-fortran library
    requires?  I know it may not be easy, but I'd like to leave the
    numerical code to be F77 (long arg lists and no common blocks!) to
    be more easily called from C (and so easier to wrap).  I don't
    know how portable doing this really is for Fortran versions beyond
    F77.  My understanding is that many F90 "objects" are not that
    easy to pass between C and Fortran, and even some F90 structures
    are apparently challenging, though I don't know the details and
    have never delved into this.
 
	That's not to say I don't completely agree that json is better
    than FEFF.inp, and it might be that the json-fortran library is no
    problem.  And we might be able to simply not need to wrap this
    code. But I'd like to be a little careful about adding
    dependencies on F90 and later.

Clearly, if we stay with json-fortran, these concerns need to be
addressed and well tested.

---

# To be fixed

1. Many of the SConstruct files explicitly link to libgfortran and
   libm.  That needs to be fixed for other platforms/compilers.

2. C wrapper's SConstruct has /usr/local/lib hardwired


---

# Simple static analysis

The following is a very useful command:

	../src> ftnchek -mkhtml */*.f

When run in the `src/` folder, this will generate a small html file
for every fortran source file below `src/` in the repository.  These
html files explain function and subroutine arguments, common block
parameters, local variables, and several other useful tidbits for each
file.

This combined with the call graphs in the `tree/` folder under each
source code folder provide a fairly thorough mapping of information
through the many parts of Feff.
