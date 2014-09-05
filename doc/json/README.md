# Dictionaries of Feff interchange files

Feff uses a slew of files to exchange data between the different parts of the program.

An early modification of feff85exafs changed the legacy structured
text files into JSON files.

old        | new         | purpose
-----------+-------------+-------------------------------------------------------
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
