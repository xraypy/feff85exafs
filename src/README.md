# Building feff85exafs

In an attempt to move beyond the inadequate collection of Makefiles
that this project inherited from the Feff Project, Bruce decided to
use the [scons software construction tool](http://www.scons.org/) as a
replacement.  scons seemed easier to understand than GNU Autotools
(which Matt used ages ago for feff6l and Ifeffit).  Because scons is
written in python and uses configuration files that contain python
code, it seemed like a good fit for something that would eventually be
integrated into [Larch](https://github.com/xraypy/xraylarch).

Once all these scons scripts are working correctly, it should be as
simple as typing `scons` after cd-ing into the `src/` folder.

You can also cd into one of the subfolders and do `scons` to build
just that piece.

The file `FortranCompilation` has a bit of logic for specifying
compilation flags, which may be attractive for optimization or other
reasons.  This is then imported by each of the lower-level SConstruct
files.  Except for gfortran, the initial guesses for compilation flags
were taken from the top level Makefile that came from the Feff
Project.  The gfortran flags were selected by Bruce after playing
around on his Ubuntu system.


# Compiling against json-fortran

At the moment, feff85exafs has been modified to replace all
intermediate i/o files (`modN.inp`, `global.dat`, `geom.dat`, and so
on) with JSON file.  Reading and writing the JSON files is currently
done using
[json-fortran](https://github.com/jacobwilliams/json-fortran).

According to the json-fortran web page, it can be compiled with
Visualk Studio 2010, Intel Fortran, and gfortran 4.9.  My only
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


Matt has these reasonmable things to say about compiling against json-fortran:

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
