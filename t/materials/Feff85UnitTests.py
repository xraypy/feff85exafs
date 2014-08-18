
from   os      import makedirs, chdir
from   os.path import realpath, isdir, join
from   shutil  import rmtree
import subprocess, pystache, json, re


import larch
larch.use_plugin_path('xafs')
from feffdat import feffpath
larch.use_plugin_path('math')
from mathutils import _deriv

import pylab
from termcolor import colored
from numpy import diff, gradient

def runfeff(folder, doscf=False, _larch=None):
    """
    Make a feff.inp from the mustache template, then run feff 8.5

    Take the name of a folder containing the test as the argument.
    """
    if not isdir(folder):
        print folder + " is not one of the available tests"
        exit

    here     = realpath(join('.'))
    testrun  = join(folder, 'testrun')
    repotop  = realpath(join('..','..'))

    if isdir(testrun): rmtree(testrun)
    makedirs(testrun)

    mat_json = json.load(open(join(folder, folder + '.json')))
    mat_json['doscf']='* '
    if doscf: mat_json['doscf']=''

    
    inp      = open(join(testrun,'feff.inp'), 'w')
    renderer = pystache.Renderer()
    inp.write(renderer.render_path( join(folder, folder + '.mustache'), # material/material.mustache
                                    mat_json ))                         # material/material.json
    inp.close

    chdir(testrun)
    f85escript = join(repotop, 'bin', 'f85e')
    subprocess.call(f85escript); # the f85e shell script emulates the behavior of the monolithic Feff application
    chdir(here)


def compare_nnnn(baseline, testrun, nnnn=1, part='mag_feff', doscf=False, _larch=None):

    #testrun  = join(realpath('.'), folder, 'testrun')
    #baseline = join(realpath('.'), folder, 'baseline')

    nnnndat = "feff%4.4d.dat" % nnnn

    blpath = feffpath(join(baseline, nnnndat))
    trpath = feffpath(join(testrun,  nnnndat))

    mbl = getattr(blpath._feffdat, 'mag_feff')
    mtr = getattr(trpath._feffdat, 'mag_feff')
    pbl = getattr(blpath._feffdat, 'pha_feff')
    ptr = getattr(trpath._feffdat, 'pha_feff')

    pylab.xlabel('wavenumber $\AA^{-1}$')
    pylab.ylabel('magnitude and phase')
    pylab.title(nnnndat)
    pylab.plot(blpath._feffdat.k, mbl,           label='magnitude of baseline')
    pylab.plot(trpath._feffdat.k, mtr,           label='magnitude of test run')
    pylab.plot(blpath._feffdat.k, gradient(pbl), label='grad(phase of baseline)')
    pylab.plot(trpath._feffdat.k, gradient(ptr), label='grad(phase of test run)')
    pylab.legend(loc='best')

    rfactor_mag = sum((mbl - mtr)**2) / sum(mbl**2)
    rfactor_pha = sum((pbl - ptr)**2) / sum(pbl**2)

    scf="without SCF"
    if doscf: scf="with SCF"

    print colored("\nComparing magnitude and phase of %s (%s)" % (nnnndat, scf), 'green', attrs=['bold'])
    print "magnitude R-factor = %.3f\nphase R-factor = %.3f" % (rfactor_mag, rfactor_mag)

    pylab.show(block=False)
    raw_input("done? ")


def clean_testrun(testrun):
    rmtree(testrun)

#def registerLarchPlugin(): # must have a function with this name!
#    return ('Feff85UnitTests', {'runfeff': runfeff, 'compare_nnnn': compare_nnnn})
