# Feff85EXAFS example: Copper metal

This is a data-bearing example, thus the tests will cover both the
computation and the results of the analysis.

This is a copper foil measured in a displex near 10K.  Since 1995,
this has been a canonical example for testing EXAFS theory and analysis.

This copper data is the foil in reference 57 of the Zabinski et al
paper from 1995 which introduced Feff's approach to multiple
scattering to the world.  According to Google Scholar, this spectrum
has been cited at least 17 times!

For an amusing walk down memory lane, here's Matt and Bruce recalling
this measurement 22 years later:
https://github.com/XraySpectroscopy/XAS-Data-Interchange/issues/29


* [Atoms input file for copper](Copper_atoms.inp)
* [template for feff input file (generated from the atoms input file)](Copper.mustache) -- uses the [mustache](http://mustache.github.io/) templating system
* [XDI data file: chi(k)](Copper.chik)
* [The raw data file](cu.012)
* [Athena project file](Copper.prj) -- Data processed with a pre-release of Demeter 0.9.21 using the larch (0.9.23) backend and XDI 1.0

