# Feff85EXAFS example: 1-bromoadamantane

This is a data-bearing example, thus the tests will cover both the
computation and the results of the analysis.

Bromoadamantane is a single-bonded, 10-carbon frame with hydrogen
atoms accommodating the remaining charge on the carbon atoms.  One
hydrogen atom is replaced with a bromine atom.  The XAS data are
measured at the Br K edge.

This is an interesting analysis problem because the Br-H scattering
plays a significant and easily measurable role in the EXAFS analysis.

The bromoadamantane sample was kindly provided by Alessandra Leri,
Marymount Manhattan College and measured by Bruce at NSLS X23A2.

![Ball and stick figure of the 1-bromoadamantane molecule](bromoadamantane.png)


* [Chemical table file (cartesian coordinates)](bromoadamantane.sdf) -- [CID 79106](http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=79106)
* [template for feff input file (converted from SDF file)](bromoadamantane.mustache) -- uses the [mustache](http://mustache.github.io/) templating system
* [XDI data file: chi(k)](bromoadamantane.chik)
* [Athena project file](bromoadamantane.prj) -- Data processed with a pre-release of Demeter 0.9.21 using the larch (0.9.23) backend and XDI 1.0

