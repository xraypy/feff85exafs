# Feff85EXAFS example: uraninite, UO2

This is a data-bearing example, thus the tests will cover both the
computation and the results of the analysis.



Bromoadamantane is a single-bonded, 10-carbon frame with hydrogen
atoms accommodating the remaining charge on the carbon atoms.  One
hydrogen atom is replaced with a bromine atom.  The XAS data are
measured at the Br K edge.

This is an interesting analysis problem because the Br-H scattering
plays a significant and easily measurable role in the EXAFS analysis.

The uraninite data was measured at APS 10ID and kindly provided by
Shelly Kelly.


![Ball and stick figure of uraninite](uraninite.gif)

![Canada &#9829; uraninite!](stamp.jpg)


* [CIF file for uraninite](UO2.cif)
* [template for feff input file (generated from CIF file)](UO2.mustache) -- uses the [mustache](http://mustache.github.io/) templating system
* [XDI data file: chi(k)](UO2.chik)
* [Athena project file](UO2.prj) -- Data processed with a pre-release of Demeter 0.9.21 using the larch (0.9.23) backend and XDI 1.0

