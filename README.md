feff85exafs
===========

Feff 8.50 for EXAFS: Open Source theoretical EXAFS calculations

Based on or developed using Distribution: FEFF8.5L
Copyright (c) [2013] University of Washington

For information, see [the feff85exafs wiki](https://github.com/xraypy/feff85exafs/wiki).

----

Compilation requires the [scons software construction tool](http://www.scons.org/).

To set compilation flags, edit the file `src/FortranCompilation.py`.

To compile, cd to the `src/` directory and type `scons`.

See [`src/README.md`](src/README.md) for details on compiling this
version of Feff, including compiling against
[json-fortran](https://github.com/jacobwilliams/json-fortran).

For the fortran entry point to the stand-alone F_eff calculation, the
C wrapper around it, or language bindings to C wrapper, see the
`wrapprs/` directory.
