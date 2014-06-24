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

# Build Issues

 1. What is the purpose of RDINP?
 
