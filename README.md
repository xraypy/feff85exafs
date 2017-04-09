Feff8L
======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20629.svg)](http://dx.doi.org/10.5281/zenodo.20629)


Feff8L: Open Source theoretical EXAFS calculations

Based on or developed using Distribution: FEFF8.5L
Copyright (c) [2013] University of Washington

* The
  [goals](https://github.com/xraypy/feff85exafs/wiki/Goals-of-the-feff85exafs-project)
  of this project are discussed on the
  [wiki](https://github.com/xraypy/feff85exafs/wiki).

Compilation uses fairly generic Makefiles.  Testing requires
[python's nose tool](https://nose.readthedocs.org/en/latest/).

To set compilation flags, edit the `Makefile`.

To compile, test, and build, do the following in the top directory

```
	  feff8l> make
	  feff8l> LD_LIBRARY_PATH='wrappers/python:src/GENFMT' && nosetests --verbosity=3
	  feff8l> sudo make install
```

See [`src/README.md`](src/README.md) for details on compiling this
version of Feff, including compiling against
[json-fortran](https://github.com/jacobwilliams/json-fortran).

For example programs using the fortran entry point to the stand-alone
F_eff calculations or programs using the C wrapper see the `wrappers/`
directory.  There you will also find language bindings to the C
wrapper, including python and perl.
