# Feff85EXAFS example: Lanthanum cuprate, polarization parallel to the z-axis

This is not a data-bearing example, thus the tests will only cover the
computation.

This is La2CuO4, an orthorhombic crystal with the oxygen octahedron
around the copper atom elongated in the z direction.  In the plane,
there are 4 Cu-O scatterers at 1.90 angstrom.  In the z direction,
there are two Cu-O scatterers at 2.43 angstrom.

This tests making a polarizion dependent Feff calculation.  See also
the LCO-perp test which uses the same material to test polarization
and ellipticipty.

Structure from
[Radaelli at al.](http://dx.doi.org/10.1103/PhysRevB.49.4163).
EXAFS analysis from
[Haskel, et al.](http://dx.doi.org/10.1103/PhysRevB.56.R521)

![Ball and stick figure of La2CuO4 from Haskel, et al. ](LCO.png)


* [Atoms input file for LCO](LCO-para.inp)
* [template for feff input file (generated from the atoms input file)](LCO-para.mustache) -- uses the [mustache](http://mustache.github.io/) templating system

