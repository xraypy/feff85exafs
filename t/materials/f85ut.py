
from   os      import makedirs, chdir, getcwd
from   os.path import realpath, isdir, join
from   shutil  import rmtree
import subprocess, pystache, json, re


from larch import (Group, Parameter, isParameter,
                   param_value, use_plugin_path, isNamedClass)
use_plugin_path('xafs')
from feffdat import feffpath
use_plugin_path('wx')
from plotter import (_newplot, _plot)
use_plugin_path('std')
from show import _show
#import pylab

from termcolor import colored
import numpy as np


class Feff85exafsUnitTestGroup(Group):

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='Feff85exafs unit test: %s' % folder)
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self.doplot     = True
        self.doscf      = False
        self.feffran    = False
        self.folder     = folder
        self.path       = realpath(folder)
        self.repotop    = realpath(join('..','..'))
        self.f85escript = join(self.repotop, 'bin', 'f85e')


    def run(self):
        """
        Make a feff.inp from the mustache template, then run feff 8.5
        
        Take the name of a folder containing the test as the argument.
        """

        if not isdir(self.folder):
            print self.folder + " is not one of the available tests"
            exit

        here     = getcwd()

        scf = 'noSCF'
        if self.doscf: scf = 'withSCF'
        self.testrun  = join(self.folder, 'testrun')


        if isdir(self.testrun): rmtree(self.testrun)
        makedirs(self.testrun)

        self.json = json.load(open(join(self.folder, self.folder + '.json')))
        self.json['doscf']='* '
        if self.doscf: self.json['doscf']=''

    
        inp      = open(join(self.testrun,'feff.inp'), 'w')
        renderer = pystache.Renderer()
        inp.write(renderer.render_path( join(self.path, self.folder + '.mustache'), # material/material.mustache
                                        self.json ))                                # material/material.json
        inp.close

        chdir(self.testrun)
        ok = subprocess.check_call(self.f85escript); # the f85e shell script emulates the behavior of the monolithic Feff application
        print colored("\nRan Feff85EXAFS on %s (%s)" % (self.folder, scf), 'yellow', attrs=['bold'])
        chdir(here)
        self.feffran = True


    def compare(self, nnnn=1):
        """
        compare a feffNNNN.dat file from the testrun with the same file from the baseline calculation
        """
        if not self.feffran:
            print colored("You have not yet run the test Feff calculation", 'magenta', attrs=['bold'])
            return

        if self.folder[-1] == '/': self.folder = self.folder[:-1]
        scf = 'noSCF'
        if self.doscf: scf = 'withSCF'
        self.baseline = join(self.folder, 'baseline', scf)
        self.testrun  = join(self.folder, 'testrun')

        nnnndat = "feff%4.4d.dat" % nnnn

        blpath = feffpath(join(self.baseline, nnnndat))
        trpath = feffpath(join(self.testrun,  nnnndat))

        mbl = getattr(blpath._feffdat, 'mag_feff') # + np.random.uniform(0,1,size=1)
        mtr = getattr(trpath._feffdat, 'mag_feff')
        pbl = getattr(blpath._feffdat, 'pha_feff') # + np.random.uniform(0,1,size=1)
        ptr = getattr(trpath._feffdat, 'pha_feff')

        self.rfactor_mag = sum((mbl - mtr)**2) / sum(mbl**2)
        self.rfactor_pha = sum((pbl - ptr)**2) / sum(pbl**2)

        scf="without SCF"
        if self.doscf: scf="with SCF"

        print colored("\nComparing magnitude and phase of %s (%s)" % (nnnndat, scf), 'green', attrs=['bold'])
        print "magnitude R-factor = %.3f\nphase R-factor = %.3f" % (self.rfactor_mag, self.rfactor_pha)

        if self.doplot:
            _newplot(blpath._feffdat.k, mbl) #,        label='magnitude of baseline', xlabel='wavenumber $\AA^{-1}$', ylabel='magnitude and phase', title=nnnndat)
            #_plot(trpath._feffdat.k, mtr) #,           label='magnitude of test run')
            #_plot(blpath._feffdat.k, np.gradient(pbl)) #, label='grad(phase of baseline)')
            #_plot(trpath._feffdat.k, np.gradient(ptr)) #, label='grad(phase of test run)')

        # if self.doplot:
        #     pylab.xlabel('wavenumber $\AA^{-1}$')
        #     pylab.ylabel('magnitude and phase')
        #     pylab.title(nnnndat)
        #     pylab.plot(blpath._feffdat.k, mbl, label='magnitude of baseline')
        #     pylab.plot(trpath._feffdat.k, mtr, label='magnitude of test run')
        #     pylab.plot(blpath._feffdat.k, gradient(pbl), label='grad(phase of baseline)')
        #     pylab.plot(trpath._feffdat.k, gradient(ptr), label='grad(phase of test run)')
        #     pylab.legend(loc='best')
        #     #pylab.ion()
        #     pylab.show()
        #raw_input("done? ")


    def clean(self):
        rmtree(self.testrun)
        self.feffran = False

def ut(folder=None, _larch=None, **kws):
    return Feff85exafsUnitTestGroup(folder=folder, _larch=_larch)

def registerLarchPlugin(): # must have a function with this name!
    return ('f85ut', { 'ut': ut })
    
