
from   os        import makedirs, chdir, getcwd
from   os.path   import realpath, isdir, join
from   shutil    import rmtree
import subprocess, glob, pystache, json, re
from   termcolor import colored
import numpy     as np

from larch import (Group, Parameter, isParameter, param_value, use_plugin_path, isNamedClass)
use_plugin_path('xafs')
from feffdat import feffpath
use_plugin_path('wx')
from plotter import (_newplot, _plot)


## :TODO:
##   1. see which paths come in the baseline, available method chould check baseline and testrun
##   2. cull s02 calculation from f85e.log or from chi.dat, test its value
##   3. tests for muffin and norman radii, mu, kf, vint, rs_int
##   4. compare() should include a report of path geometry, ie feffpath._feffdat.geom
##   5. begin writing data/fitting tests
##   6. need a way of validating the R-factors, that is, changes may introduce negligible numerical differences -- define negligible
##   7. capture and interpret feff's screen messages to notice when a feff run fails
##   8.    "     "      "       "       "       "    to use number of SCF iterations as a unit test
##   9. someting like perl's Term::Twiddle would be nice to use when self.quiet is True


class Feff85exafsUnitTestGroup(Group):

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='Feff85exafs unit test: %s' % folder)
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self._larch     = _larch
        self.doplot     = True  # True = plot comparisons
        self.doscf      = False # True = use self-consistency
        self.quiet      = False # True = suppress Feff's screen messages
        self.feffran    = False # True = Feff calculation has been run
        self.folder     = folder
        self.testrun    = join(self.folder, 'testrun')
        if not isdir(folder):
            print colored(folder + " is not one of the available tests", 'magenta', attrs=['bold'])
            return None
        self.path       = realpath(folder)
        self.repotop    = realpath(join('..','..'))
        # the f85e shell script emulates the behavior of the monolithic Feff application
        self.f85escript = join(self.repotop, 'bin', 'f85e')

    def __repr__(self):
        if not isdir(self.folder):
            return '<Feff85exafs Unit Test Group (empty)>'
        if self.folder is not None:
            return '<Feff85exafs Unit Test Group: %s>' % self.folder
        return '<Feff85exafs Unit Test Group (empty)>'


    def run(self, _larch=None):
        """
        Make a feff.inp from the mustache template, then run feff 8.5
        
        Take the name of a folder containing the test as the argument.
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

    
        inp      = open(join(self.testrun,'feff.inp'), 'w')
        renderer = pystache.Renderer()
        inp.write(renderer.render_path( join(self.path, self.folder + '.mustache'), # material/material.mustache
                                        self.json ))                                # material/material.json
        inp.close

        here     = getcwd()
        chdir(self.testrun)
        if self.quiet:
            outout = subprocess.check_output(self.f85escript);
        else :
            subprocess.check_call(self.f85escript);

        print colored("\nRan Feff85EXAFS on %s (%s)" % (self.folder, scf), 'yellow', attrs=['bold'])
        #print "ok = %s" % ok

        self.paths = list()
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
            if self.paths.count(nnnn):
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

        if self.folder[-1] == '/': self.folder = self.folder[:-1]
        scf = 'noSCF'
        if self.doscf: scf = 'withSCF'
        self.baseline = join(self.folder, 'baseline', scf)
        self.testrun  = join(self.folder, 'testrun')

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
            print colored("Unknown choice of parts \"%s\"\nmust be one of (feff|amp|phase|lambda|caps|redfact|rep)\nusing feff" % part,
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
        print colored("\nComparing %s of %s (%s)" % (label, nnnndat, scf), 'green', attrs=['bold'])
        print label + " R-factor = %.3f" % self.rfactor
        if part=='feff':
            self.rfactor_2 = sum((testrun_2 - testrun_2)**2) / sum(baseline_2**2)
            print "phase R-factor = %.3f" % self.rfactor_2
        print "\n"

        if self.doplot:
            _newplot(blpath._feffdat.k, baseline_1, _larch=self._larch, label=label+' of test run',
                     xlabel='wavenumber $\AA^{-1}$', ylabel=ylabel, title=nnnndat, show_legend=True, legend_loc='best')
            _plot   (trpath._feffdat.k, testrun_1,  _larch=self._larch, label=label+' of test run')
            if part=='feff':
                _plot(blpath._feffdat.k, np.gradient(baseline_2), _larch=self._larch, label='grad(phase of baseline)')
                _plot(trpath._feffdat.k, np.gradient(testrun_2),  _larch=self._larch, label='grad(phase of test run)')



    def clean(self):
        rmtree(self.testrun)
        self.feffran = False

######################################################################

def ut(folder=None, _larch=None, **kws):
    return Feff85exafsUnitTestGroup(folder=folder, _larch=_larch)

def registerLarchPlugin(): # must have a function with this name!
    return ('f85ut', { 'ut': ut })
    
