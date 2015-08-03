# Feff85EXAFS example: NiO

This is a data-bearing example, thus the tests will cover both the
computation and the results of the analysis.

NiO is a simple transition metal oxide in the rock-salt structure.
This is intended as a slightly more complex example than copper metal.
The XAS data are measured at the Ni K edge.

The NiO sample was kindly provided by Neil Hyatt and Martin Stennett,
Sheffield University and measured by Bruce at NSLS X23A2.
See [DOI: 10.1107/S1600577515013521](http://dx.doi.org/10.1107/S1600577515013521).

![Ball and stick figure of rock-salt NiO](NiO.png)

* [Atoms input file file for NiO](NiO_atoms.inp)
* [template for feff input file (generated from the atoms input file)](NiO.mustache) -- uses the [mustache](http://mustache.github.io/) templating system
* [XDI data file: chi(k)](NiO.chik)
* [Athena project file](NiO.prj) -- Data processed with a pre-release of Demeter 0.9.21 using the larch (0.9.24) backend and XDI 1.0

