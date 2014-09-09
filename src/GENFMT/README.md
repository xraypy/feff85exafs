
# Content of the GENFMT folder

This directory contains routines for generalized F matrix
calculations, as described in
[the 1990 Rehr-Albers paper](http://dx.doi.org/10.1103/PhysRevB.41.8139).

All routines in this directory are covered by the [LICENSE](../HEADERS/license.h)

# Simple static analysis

To make HTML files explaining data I/O for each fortran source file, do

	../src> ftnchek -mkhtml *.f

# Call graph

![call graph for the GENFMT folder](tree/genfmt.png)
