# Dictionaries of Feff interchange files

Feff uses a slew of files to exchange data between the different parts of the program.

An early modification of feff8l changed the legacy structured
text files into JSON files.

old        | new         | purpose
---------- | ----------- | ------------------------------------------------------
mod1.inp   | pot.json    | input data for POT module taken from feff.inp
mod2.inp   | xsph.json   | input data for XSPH module taken from feff.inp
mod4.inp   | path.json   | input data for PATH module taken from feff.inp
mod5.inp   | genfmt.json | input data for GENFMT module taken from feff.inp
mod6.inp   | ff2x.json   | input data for FF2X module taken from feff.inp
global.dat | global.json | input data common to all modules, taken from feff.inp
atoms.dat  | atoms.dat   | geometric information about ATOMS list
geom.dat   | geom.dat    | geometric information about ATOMS list
xsect.bin  | xsect.json  | energy grid and cross-section data


The files in this folder are JSON files containing dictionaries of the
content of each of the interchange files.  Hopefully this will help
application developers who need to read from or write to these
interchange files.


## Other ascii interchange files that do not yet have json files

Some of these are throw-away diagnostic files

1. nstar.dat, see `GENFMT/genfmt.f`, but nstar.dat seems to be written but not read in feff8l
2. wscrn.dat, see `XSPH/xsph.f`
3. exc.dat + sigma.dat + ratio.dat + ratiop.dat + mpse.dat, see `XSPH/xsect.f` + `XSPH/phase.f`
4. various in `XSPH/rholat.f`
5. emesh.dat in `XSPH/phmesh.f`
6. axafs.dat in `XSPH/axafs.f`
7. various things written with PRINT cards, see comments lines 319 to 342 in `RDINP/rdinp_l.f`
8. rotmat.dat in `GENFMT/rot3i.f`
9. various in `FMS/xprep.f`, `FMS/inc.f`
10. various in `FF2X/ff2gen.f`
11. exc.dat in `EXCH/mpse.f`
12. various in `DEBYE/sigrem.f`
13. `spring.inp` in  `DEBYE/sigrem.f`
