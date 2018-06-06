Feff8L
======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20629.svg)](https://doi.org/10.5281/zenodo.20629)

Feff8L: Open Source theoretical EXAFS calculations

Based on or developed using Distribution: FEFF8.5L Copyright (c) [2013] University of Washington

 * The [goals](https://github.com/xraypy/feff85exafs/wiki/Goals-of-the-feff85exafs-project) of this
   project are discussed on the [wiki](https://github.com/xraypy/feff85exafs/wiki).

 * This project is sometimes called Feff8L or Feff85exafs.  The "L" implies "Lite" (as in diluted
   beer, not illumination), "85" comes from the version of Feff (8.5) that this version is derived
   from, and "exafs" is because this version is used to calculate EXAFS only -- not XANES or other
   core-level spectroscopies.  Using either Feff8L or feff85exafs or some combination of these is
   fine with us.

* Build and Installation instructions are given in [`Install.md`](Install.md), with more details in
 [`src/README.md`](src/README.md) .

* There are extensive tests for Feff8L in the [`tests`](tests) directory.  To run these, you will
  need [Larch](https://github.com/xraypy/xraylarch)

* Example programs using the fortran entry point to the stand-alone F_eff calculations or programs
  using the C wrapper can be found in the `wrappers/` directory.  There you will also find language
  bindings to the C wrapper, including python and perl.
