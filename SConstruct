## -*- python -*-
## feff85exafs build system based on scons
## see HEADERS/license.h for license information

import sys
sys.path.append( 'src' )
from FeffBuild import CompilationEnvironment, InstallEnvironment

import os

env  = CompilationEnvironment()

#env = Environment()
print '\n\tFortran compiler is ' + env['FORTRAN']
ienv = InstallEnvironment()
print "\tinstallation prefix is: " + ienv['i_prefix'] +"\n"


fortran = ['src/PAR/SConstruct',
            'src/COMMON/SConstruct',
            'src/json-fortran/SConstruct',
            'src/JSON/SConstruct',
            'src/MATH/SConstruct',
            'src/ATOM/SConstruct',
            'src/DEBYE/SConstruct',
            'src/EXCH/SConstruct',
            'src/FOVRG/SConstruct',
            'src/FMS/SConstruct',
            'src/RDINP/SConstruct',
            'src/OPCONSAT/SConstruct',
            'src/POT/SConstruct',
            'src/XSPH/SConstruct',
            'src/PATH/SConstruct',
            'src/GENFMT/SConstruct',
	    'src/FF2X/SConstruct',]
python   = [
    #'wrappers/python/SConstruct',
    'tests/SConstruct',
]

if os.name == 'nt':
    everything = fortran
else:
    everything = fortran + python

SConscript(everything)


#env = Environment(BUILDERS = {'MyBuild' : b})
#env.MyBuild(join(home, larch_dir, 'plugins', 'f85ut.py'), join('tests', 'f85ut.py'))


# RDINP  -> feff.inp reader
# POT    -> module 1
# XSPH   -> module 2
# FMS    -> module 3, executable is not a part of feff85exafs
# PATH   -> module 4
# GENFMT -> module 5
# FF2X   -> module 6
