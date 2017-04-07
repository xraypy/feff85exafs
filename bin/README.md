# Executable scripts for feff8l

## f85e

This is a bash script which emulates the behavior of a monolithic Feff
by running the "modules" one after another given a `feff.inp` file.

This script is used by the unit testing framework.

In its current form, it deduces the locations of the various
stand-alone executables for the modules based on the structure of the
repository.  That is, this script assumes that is lives in the `bin/`
folder of the repository and that the executables will be found in
subfolders of the `src` folder, where `bin/` and `src/` are at the
same level of the repository folder structure.

	~/feff8l/> ls -R
	./bin:
	f85e*

	./src:
	ATOM/    EXCH/  FOVRG/    JSON/  PAR/   RDINP/  FortranCompilation.py
	COMMON/  FF2X/  GENFMT/   LIB/   PATH/  XSPH/   SConstruct
	DEBYE/   FMS/   HEADERS/  MATH/  POT/   attic/  README.md

	and so on...

Eventually, an installation scripts might be written, in which case a
successor to this script will need to know (or determine) the
installation locations of the various executables.

For now, I am assuming that primary purpose of the `f85e` script is to
aid unit testing during the development of Feff8L.

# feff85_light

This is a script written by Josh.  It makes calls to executables
called `opconsat`, `atomic`, and `aps2exc`, none of which I know how
to build from the source code provided.  The text of the script
suggests that Josh intended this to be an example of using
Feff8L, but its use is unclear.

