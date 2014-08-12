List of materials for unit testing
==================================

## copper metal

Because copper

  * `cu10k.chi`: copper metal data at 10K
  * `cu_atoms.inp`: atoms input file
  * `cu.inp`: feff8 input file

## nickel oxide

This is a simple, cubic, metal oxide

  * `NiO_atoms.inp`: atoms input file
  * `NiO.inp`: feff8 input file
  
## uraninite

This is an f-electron system

  * `UO2.cif`: CIF file
  * `UO2.inp`: feff8 input file
  * `UO2.chik`: chi(k) data file
  
## zircon

This has a 4d metal and Si, good for testing at each edge

  * `ZrSiO4.inp`: atoms input file
  * `zircon.inp`: feff8 input file
  
## bromoadamantane

This is a small molecule which can be fit with a fairly simple model
of four paths, but for which the 6 nearby hydrogen scatterers seem to
play a big role in the fit.

  * `bromoadamantane.inp`: feff input file
  * `bromoadamantane.chik`: chi(k) data for bromoadamantane
  * `bromoadamantane.png`: a picture of the bromoadamantane molecule

## ferrocene macrocycle

This is a iron organometallic from Dinnebier et al, Organometallics
2001, 20, 5642-5647, doi:10.1021/om0105066

It has Fe in between two 5 member carbon rings, so it has 10 C
neighbors at a range of distances from 2.026 A to 2.099 A

  * `ferrocene-macrocycle_atoms.inp`: atoms input file
  * `ferrocene-macrocycle.inp`: feff8 input file
