
## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

from   os        import makedirs, chdir, getcwd
from   os.path   import realpath, isdir, join
from   shutil    import rmtree
import sys, subprocess, glob, pystache, json, re
from   termcolor import colored
import numpy     as np
import importlib

from larch import (Group, Parameter, isParameter, param_value, use_plugin_path, isNamedClass, Interpreter)
use_plugin_path('xafs')
from feffdat import feffpath
use_plugin_path('wx')
from plotter import (_newplot, _plot)

class Feff85exafsUnitTestGroup(Group):
    """
    A group for performing unit tests on feff85exafs.

    Methods:
       run        : run the test feff calculation, no return value
       testpaths  : fill the paths attribute, no return value
       available  : returns True is a path index has a corresponding feffNNNN.dat file
       compare    : make a comparison of columns in feffNNNN.dat, returns True is no difference between test and baseline
       geometry   : write a description of scattering path to the screen
       radii      : fetch a list of muffin tin or norman radii for the unique potentials
       s02        : fetch the calculated value of s02 from the testrun or the baseline
       feffterms  : perform a test on various values in the header of feffNNNN.dat, returns True if no difference
       clean      : remove the testrun folder
    

    Attributes:
       doplot    :  boolean, True = make plots in compare method
       doscf     :  boolean, True = use self-consistency when running feff for unit tests
       verbose   :  boolean, True = write screen messages when running feff and performing tests
       feffran   :  boolean, True = feff has been run and testrun folder holds the output
       folder    :  string,  name of folder containing test materials
       testrun   :  string,  name of folder containing output of feff test run
       baseline  :  string,  name of folder containing the baseline feff calculation
       paths     :  list,    list of feffNNNN.dat files in testrun folder
       bpaths    :  string,  list of feffNNNN.dat files from baseline calculation
       path      :  string,  fully resolved path to folder
       repotop   :  string,  fully resolved path to top of feff85exafs repository
       json      :  json string used to configure the test feff run
       f85script :  string,  fully resolved path to the f85e script, which emulates monolithic feff
       rfactor   :  float,   R-factor computed from feffNNNN.dat columns in testrun compared to baseline
       rfactor_2 :  float,   second R-factor, used when compare called with part='feff'
       epsilon   :  float,   value for comparing columns from feffNNNN.dat with the baseline and other things
       count, datacount, feffcount : count number of tests
    """

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='Feff85exafs unit test: %s' % folder)
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self._larch     = Interpreter()
        self.doplot     = True  
        self.doscf      = False # True = use self-consistency
        self.verbose    = True  # True = print Feff's screen messages and other screenmessages
        self.feffran    = False # True = Feff calculation has been run
        self.count      = 0
        self.feffcount  = 0
        self.datacount  = 0
        self.failed     = list()
        self.folder     = folder
        if self.folder[-1] == '/': self.folder = self.folder[:-1]
        self.testrun    = realpath(join(self.folder, 'testrun'))
        self.testpaths()
        if not isdir(folder):
            print colored(folder + " is not one of the available tests", 'magenta', attrs=['bold'])
            return None
        self.path       = realpath(folder)
        self.repotop    = realpath(join('..','..'))
        # the f85e shell script emulates the behavior of the monolithic Feff application
        self.f85escript = join(self.repotop, 'bin', 'f85e')
        self.epsilon    = 0.00001
        self.epsfit     = 0.001


    def __repr__(self):
        if not isdir(self.folder):
            return '<Feff85exafs Unit Test Group (empty)>'
        if self.folder is not None:
            return '<Feff85exafs Unit Test Group: %s>' % self.folder
        return '<Feff85exafs Unit Test Group (empty)>'


    @property
    def baseline(self):
        scf = 'noSCF'
        if self.doscf: scf = 'withSCF'
        return realpath(join(self.folder, 'baseline', scf))

    @property
    def bpaths(self):
        """
        Gather a list of feffNNNN.dat files from the baseline calculation
        """
        here = getcwd()
        chdir(self.baseline)
        p = list()
        feffoutput = glob.glob("*")
        for f in sorted(feffoutput):
            tosave = re.compile("feff\d+\.dat")
            if tosave.match(f):
                p.append(f)
        chdir(here)
        return p

    def run(self):
        """
        Make a feff.inp from the mustache template, then run feff 8.5 in t he testrun folder
        """
        if not isdir(self.folder):
            print colored(self.folder + " is not one of the available tests", 'magenta', attrs=['bold'])
            return False

        if isdir(self.testrun): rmtree(self.testrun)
        makedirs(self.testrun)

        scf = 'without SCF'
        self.json = json.load(open(join(self.folder, self.folder + '.json')))
        self.json['doscf']='* '
        if self.doscf:
            scf = 'with SCF'
            self.json['doscf']=''
    
        renderer = pystache.Renderer()
        with open(join(self.testrun,'feff.inp'), 'w') as inp:
            inp.write(renderer.render_path( join(self.path, self.folder + '.mustache'), # material/material.mustache
                                            self.json ))                                # material/material.json

        here     = getcwd()
        chdir(self.testrun)
        if self.verbose:
            subprocess.check_call(self.f85escript);
        else :
            outout = subprocess.check_output(self.f85escript);

        if self.verbose: print colored("\nRan Feff85EXAFS on %s (%s)" % (self.folder, scf), 'yellow', attrs=['bold'])
        self.testpaths()

        chdir(here)
        self.feffran = True


    def testpaths(self):
        """
        Gather a list of feffNNNN.dat files from the testrun
        """
        self.paths = list()
        if isdir(self.testrun):
            here = getcwd()
            chdir(self.testrun)
            feffoutput = glob.glob("*")
            for f in sorted(feffoutput):
                tosave = re.compile("feff\d+\.dat")
                if tosave.match(f):
                    self.paths.append(f)
            chdir(here)
            self.feffran = True


    def available(self, which=0):
        """
        Test whether a path index was included in the test calculation.

           larch> print group.available(1)
           True

           larch> print group.available(112)
           False

        If the argument is 0 or negative, print a list of all available paths files and return None
        """
        if which>0:
            nnnn = "feff%4.4d.dat" % which
            if self.paths.count(nnnn) and self.bpaths.count(nnnn):
                return True
            else:
                return False
        else:
            print "The following paths are available:"
            for f in self.paths:
                print "\t"+f
            return None


    def compare(self, nnnn=1, part='feff', _larch=None):
        """
        compare a feffNNNN.dat file from the testrun with the same file from the baseline calculation

            group.compare(N, part)

        where N is the path index and part is one of 
            feff    :      test the magnitude AND phase of F_eff (default)
            amp     :      test the total amplitude of the path
            phase   :      test the total phase of the path
            lambda  :      test the mean free path
            caps    :      test the central atom phase shift
            redfact :      test the reduction factor
            rep     :      test the real part of the complex wavenumber

        """
        if self._larch is None:
            raise Warning("cannot do path comparison -- larch broken?")

        if not self.feffran:
            print colored("You have not yet run the test Feff calculation", 'magenta', attrs=['bold'])
            return

        if not self.available(nnnn):
            print colored("Path %d was not saved from the test Feff calculation" % nnnn, 'magenta', attrs=['bold'])
            return

        nnnndat = "feff%4.4d.dat" % nnnn

        blpath = feffpath(join(self.baseline, nnnndat))
        trpath = feffpath(join(self.testrun,  nnnndat))

        if part=='feff':
            baseline_1 = getattr(blpath._feffdat, 'mag_feff') # + np.random.uniform(0,1,size=1)
            testrun_1  = getattr(trpath._feffdat, 'mag_feff')
            baseline_2 = getattr(blpath._feffdat, 'pha_feff') # + np.random.uniform(0,1,size=1)
            testrun_2  = getattr(trpath._feffdat, 'pha_feff')
            ylabel     = 'magnitude and phase'
            label      = 'magnitude' 
        elif part=='amp':
            baseline_1 = getattr(blpath._feffdat, 'amp')
            testrun_1  = getattr(trpath._feffdat, 'amp')
            ylabel     = 'total amplitude'
            label      = 'amplitude' 
        elif part=='phase':
            baseline_1 = getattr(blpath._feffdat, 'pha')
            testrun_1  = getattr(trpath._feffdat, 'pha')
            ylabel     = 'total phase shift'
            label      = 'phase' 
        elif part=='lambda':
            baseline_1 = getattr(blpath._feffdat, 'lam')
            testrun_1  = getattr(trpath._feffdat, 'lam')
            ylabel     = 'mean free path'
            label      = 'MFP' 
        elif part=='caps':
            baseline_1 = getattr(blpath._feffdat, 'real_phc')
            testrun_1  = getattr(trpath._feffdat, 'real_phc')
            ylabel     = 'central atom phase shift'
            label      = 'CAPS' 
        elif part=='redfact':
            baseline_1 = getattr(blpath._feffdat, 'red_fact')
            testrun_1  = getattr(trpath._feffdat, 'red_fact')
            ylabel     = 'reduction factor'
            label      = 'reduction factor'
        elif part=='rep':
            baseline_1 = getattr(blpath._feffdat, 'rep')
            testrun_1  = getattr(trpath._feffdat, 'rep')
            ylabel     = 'real part of p(k)'
            label      = 'Re[p(k)]'
        else:
            if self.verbose: print colored("Unknown choice of parts \"%s\"\nmust be one of (feff|amp|phase|lambda|caps|redfact|rep)\nusing feff" % part,
                                           'magenta', attrs=['bold'])
            part       = 'feff'
            baseline_1 = getattr(blpath._feffdat, 'mag_feff')
            testrun_1  = getattr(trpath._feffdat, 'mag_feff')
            baseline_2 = getattr(blpath._feffdat, 'pha_feff')
            testrun_2  = getattr(trpath._feffdat, 'pha_feff')
            ylabel     = 'magnitude and phase'
            label      = 'magnitude' 
 

        scf="without SCF"
        if self.doscf: scf="with SCF"

        self.rfactor_2 = 0
        self.rfactor = sum((baseline_1 - testrun_1)**2) / sum(baseline_1**2)
        if self.verbose: 
            print colored("\nComparing %s of %s (%s)" % (label, nnnndat, scf), 'yellow', attrs=['bold'])
            self.geometry(blpath)

        if self.verbose: 
            color = 'green'  if self.rfactor < self.epsilon else 'magenta'
            print label + " R-factor = " + colored("%.3f" % self.rfactor, color, attrs=['bold'])
        if part=='feff':
            self.rfactor_2 = sum((testrun_2 - testrun_2)**2) / sum(baseline_2**2)
            if self.verbose: 
                color = 'green'  if self.rfactor_2 < self.epsilon else 'magenta'
                print "phase R-factor = " + colored("%.3f" % self.rfactor_2, color, attrs=['bold'])
        if self.verbose: print ""

        if self.doplot:
            _newplot(blpath._feffdat.k, baseline_1, _larch=self._larch, label=label+' of test run',
                     xlabel='wavenumber $\AA^{-1}$', ylabel=ylabel, title=nnnndat, show_legend=True, legend_loc='best')
            _plot   (trpath._feffdat.k, testrun_1,  _larch=self._larch, label=label+' of test run')
            if part=='feff':
                _plot(blpath._feffdat.k, np.gradient(baseline_2), _larch=self._larch, label='grad(phase of baseline)')
                _plot(trpath._feffdat.k, np.gradient(testrun_2),  _larch=self._larch, label='grad(phase of test run)')

        if part=='feff':
            return self.rfactor < self.epsilon and self.rfactor_2 < self.epsilon
        else:
            return self.rfactor < self.epsilon


    def geometry(self, path):
        """
        Print out a table showing the scattering geometry of a path
        """
        print "path (reff=%.4f nlegs=%d) geometry:" % (path.reff, path.nleg)
        for atom in path.geom:
            print "\t%-2s  %7.4f  %7.4f  %7.4f  %d" % (atom[0], atom[3], atom[4], atom[5], atom[2])
        print "\t%-2s  %7.4f  %7.4f  %7.4f  %d\n" % (path.geom[0][0], path.geom[0][3], path.geom[0][4],
                                                     path.geom[0][5], path.geom[0][2])



    def s02(self, folder='testrun'):
        """
        Fetch the value of S02 calculated by the testrun or the baseline.
           testrun:    larch> group.s02()
           baseline:   larch> group.s02('baseline')
        """
        if folder == 'testrun':
            chidat = join(self.testrun, 'chi.dat')
        else:
            chidat = join(self.baseline, 'chi.dat')

        with open(chidat, 'r') as searchfile:
            for line in searchfile:
                m = re.search('S02=(\d\.\d+)', line)
                if m:
                    return float(m.group(1))
        return -1


    def radii(self, folder='testrun', radius='muffintin'):
        """
        Return a list of the selected radii for each unique potential
           larch> group.s02(folder, radius)
        where folder = testrun|baseline and radius = norman|muffintin

        """
        if folder == 'testrun':
            chidat = join(self.testrun, 'chi.dat')
        else:
            chidat = join(self.baseline, 'chi.dat')

        values = list()
        with open(chidat, 'r') as searchfile:
            for line in searchfile:
                if radius == 'norman':
                    m = re.search('Rnm=\s*(\d\.\d+)', line)
                else:
                    m = re.search('Rmt=\s*(\d\.\d+)', line)
                if m:
                    values.append(float(m.group(1)))
        return values


    def feffterms(self, nnnn=1):
        """
        Fetch the various numbers in the feffNNNN.dat header for both the
        baseline and testrun.  A table of these values is printed.
        Return True is all are equal.
        """
        nnnndat = "feff%4.4d.dat" % nnnn
        bl      = feffpath(join(self.baseline, nnnndat))
        tr      = feffpath(join(self.testrun,  nnnndat))

        termdict = {'edge':   'energy threshold relative to atomic value',
                    'gam_ch': 'core level energy width',
                    'kf':     'k value at Fermi level',
                    'mu':     'Fermi level in eV',
                    'rs_int': 'interstitial radius',
                    'vint':   'interstitial potential'}
        if self.verbose:
            print "%-8s %-42s   %10s  %10s  %s" % ('key', 'component', 'baseline', 'testrun', 'same?')
            print "-----------------------------------------------------------------------------------"
        ok = True
        for key in termdict:
            same = getattr(bl._feffdat, key) == getattr(tr._feffdat, key)
            if self.verbose: print "%-6s   %-42s : %10s  %10s  %s" % (key, termdict[key], getattr(bl._feffdat, key),
                                                                      getattr(tr._feffdat, key), same)
            ok = ok and same
        return ok



    def fit(self):
        """
        Perform a canned fit using the baseline and testrun Feff calculations
        """
        sys.path.append(self.folder)
        module = importlib.import_module(self.folder, package=None)

        print "\tfitting %s using %s" % (self.folder, 'baseline')
        self.blfit = module.do_fit(self, 'baseline')
        print "\tfitting %s using %s" % (self.folder, 'testrun')
        self.trfit = module.do_fit(self, 'testrun')


    def clean(self):
        """
        Remove testrun folder
        """
        rmtree(self.testrun)
        self.feffran = False


######################################################################

def ut(folder=None, _larch=None, **kws):
    """
    Make a Feff85exafsUnitTestGroup group given a folder containing a baseline calculation
    """
    return Feff85exafsUnitTestGroup(folder=folder, _larch=_larch)

def ir(folder=None, _larch=None, **kws):
    """
    _i_mport and _r_un
    """
    utobj=ut(folder)
    utobj.run()
    return utobj

def irc(folder=None, _larch=None, **kws):
    """
    _i_mport, _r_un, and _c_ompare
    """
    utobj=ut(folder)
    utobj.run()
    utobj.compare(1)
    return utobj

def irf(folder=None, _larch=None, **kws):
    """
    _i_mport, _r_un, and _f_it
    """
    utobj=ut(folder)
    utobj.run()
    utobj.fit()
    return utobj

def registerLarchPlugin(): # must have a function with this name!
    return ('f85ut', { 'ut': ut, 'ir': ir, 'irc': irc, 'irf': irf })
    
