feff85exafs
===========

Feff 8.50 for EXAFS: Open Source theoretical EXAFS calculations

Based on or developed using Distribution: FEFF8.5L
Copyright (c) [2013] University of Washington

* The
  [goals](https://github.com/xraypy/feff85exafs/wiki/Goals-of-the-feff85exafs-project)
  of this project are discussed on the
  [wiki](https://github.com/xraypy/feff85exafs/wiki).

Compilation requires the [scons software construction tool](http://www.scons.org/).
Testing requires [python's nose tool](https://nose.readthedocs.org/en/latest/).

To set compilation flags, edit the file `src/FeffBuild.py`.

To compile, test, and build, do the following in the top directory

```
	  feff85exafs> scons
	  feff85exafs> LD_LIBRARY_PATH='wrappers/python:src/GENFMT' && nosetests --verbosity=3
	  feff85exafs> scons install
```

See [`src/README.md`](src/README.md) for details on compiling this
version of Feff, including compiling against
[json-fortran](https://github.com/jacobwilliams/json-fortran).

For the fortran entry point to the stand-alone F_eff calculation, the
C wrapper around it, or language bindings to C wrapper, see the
`wrappers/` directory.
