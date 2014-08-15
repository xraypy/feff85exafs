
from   os      import makedirs, chdir
from   os.path import realpath, isdir, join
from   shutil  import rmtree
import subprocess, pystache, json, re


import larch
larch.use_plugin_path('xafs')
from feffdat import feffpath

import pylab

def runfeff(folder, _larch=None):
    """
    Make a feff.inp from the mustache template, then run feff 8.5

    Take the name of a folder containing the test as the argument.
    """
    if not isdir(folder):
        print folder + " is not one of the available tests"
        exit

    here     = realpath(join('.'))
    testrun  = join(folder, 'testrun2')
    repotop  = realpath(join('..','..'))

    if isdir(testrun): 
        rmtree(testrun)
    makedirs(testrun)
    
    inp      = open(join(testrun,'feff.inp'), 'w')
    renderer = pystache.Renderer()
    inp.write(renderer.render_path( join(folder, folder + '.mustache'),                # material/material.mustache
                                    json.load(open(join(folder, folder + '.json'))) )) # material/material.json
    inp.close

    chdir(testrun)
    f85escript = join(repotop, 'bin', 'f85e')
    subprocess.call(f85escript); # the f85e shell script emulates the behavior of the monolithic Feff application
    chdir(here)


def compare_nnnn(folder, n=1, part='mag_feff', _larch=None):
    if not isdir(folder):
        print folder + " is not one of the available tests"
        exit

    here     = realpath(join('.'))
    testrun  = join(realpath('.'), folder, 'testrun')
    baseline = join(realpath('.'), folder, 'baseline')

    nnnn = "feff%4.4d.dat" % n

    blpath = feffpath(join(baseline, nnnn))
    trpath = feffpath(join(testrun,  nnnn))

    y1 = getattr(blpath._feffdat, part)
    y2 = getattr(trpath._feffdat, part)

    pylab.plot(blpath._feffdat.k, y1)
    pylab.plot(trpath._feffdat.k, y2)

    rfactor = sum((y1 - y2)**2) / sum(y1**2)

    word = re.compile("pha")
    which = 'magnitude'
    if word.match(part):
        which="phase"
    print which + " R-factor = %.3f" % rfactor

    pylab.show()
    #pause("= hit return to continue")
    
#def registerLarchPlugin(): # must have a function with this name!
#    return ('Feff85UnitTests', {'runfeff': runfeff, 'compare_nnnn': compare_nnnn})
