
# Content of the EXCH folder

This directory contains various routines to calculate self-energy
and exchange-correlation potentials.

All routines in this directory are covered by the [LICENSE](../HEADERS/license.h)

[Ankudinov and Rehr, Non-Local Self-Energy Models for XAS](https://doi.org/10.1051/jp4/1997099)
[Rehr, Kas, Prange, Sorini, Takimoto, Vila, Ab initio theory and calculations of X-ray spectra](http://dx.doi.org/10.1016/j.crhy.2008.08.004)

# Simple static analysis

To make HTML files explaining data I/O for each fortran source file, do

	../src> ftnchek -mkhtml *.f

# Call graph

![call graph for the EXCH folder](tree/EXCH.png)
