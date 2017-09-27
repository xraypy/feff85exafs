
## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

from   os        import makedirs, chdir, getcwd, unlink, listdir
from   os.path   import realpath, isdir, isfile, join
from   shutil    import rmtree, copy
import sys, subprocess, glob, pystache, json, re
from   termcolor import colored
import numpy     as np
import importlib
from   distutils.spawn import find_executable

from larch import (Group, Parameter, isParameter, param_value,
                   isNamedClass, Interpreter)
from larch_plugins.xafs.feffdat import feffpath
from larch_plugins.xafs.feffrunner import feffrunner
from larch_plugins.wx import (_newplot, _plot)

WRAPPER_AVAILABLE = False
try:
    from larch_plugins.xafs.feff8lpath import Feff8L_XAFSPath
    WRAPPER_AVAILABLE = True
except:
    pass


WARN_COLOR = 'yellow'
WARN_COLOR = 'blue'

def print_error(text):
    print(colored(text, 'magenta', attrs=['bold']))

def print_warn(text):
    print(colored(text, WARN_COLOR, attrs=['bold']))

def print_good(text):
    print(colored(text, 'green', attrs=['bold']))

def test_text(text, cond):
    color = 'green' if cond else 'magenta'
    return colored(text, color, attrs=['bold'])

class Feff85exafsUnitTestGroup(Group):
    """
    A group for performing unit tests on feff85exafs.

    Activity methods:
       run            : run the test feff calculation, no return value
       fit            : run fits using the baselline and new Feff caclautions
       clean          : remove the testrun folder

    Testing methods:
       available      : returns True is a path index has a corresponding feffNNNN.dat file
       compare        : make a comparison of columns in feffNNNN.dat, returns True if no difference
       feffterms      : perform a test on values in the header of feffNNNN.dat, ret. True if no difference
       radii          : fetch a list of muffin tin or norman radii for the unique potentials
       s02            : fetch the calculated value of s02 from the testrun or the baseline
       print_geometry : write a description of scattering path to the screen

    Attributes:
       doplot    :  boolean, True = make plots in compare method
       doscf     :  boolean, True = use self-consistency when running feff for unit tests
       verbose   :  boolean, True = write screen messages when running feff and performing tests
       feffran   :  boolean, True = feff has been run and testrun folder holds the output
       folder    :  string,  name of folder containing test materials
       testrun   :  string,  name of folder containing output of feff test run
       fefflog   :  string,  name of log file from feff run
       baseline  :  string,  name of folder containing the baseline feff calculation
       paths     :  list,    list of feffNNNN.dat files in testrun folder
       bpaths    :  list,    list of feffNNNN.dat files from baseline calculation
       path      :  string,  fully resolved path to folder
       repotop   :  string,  fully resolved path to top of feff85exafs repository
       json      :  json string, used to configure the test feff run
       rfactor   :  float,   R-factor computed from feffNNNN.dat column comparison, set by compare()
       rfactor_2 :  float,   second R-factor, used when compare called with part='feff', set by compare()
       epsilon   :  float,   value for comparing columns from feffNNNN.dat with the baseline and other things
       count     \\
       datacount  > integers, used to count number of tests run
       feffcount /
    """

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='Feff85exafs unit test: %s' % folder)
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self._larch     = Interpreter()
        self.doplot     = True
        self.doscf      = False # True = use self-consistency
        self.verbose    = True  # True = print Feff's screen messages and other screen messages
        self.feffran    = False # True = Feff calculation has been run
        self.count      = 0
        self.feffcount  = 0
        self.datacount  = 0
        self.failed     = list()
        if folder[-1] == '/': folder = folder[:-1] # strip trailing /
        self.folder     = folder

        if not isdir(folder):
            folder = join('tests', folder)
        if not isdir(folder):
            print_error(folder + " isn't one of the available tests")
            return None

        self.path       = realpath(folder)
        self.testrun    = realpath(join(self.path, 'testrun'))
        self.fefflog    = realpath(join(self.path, 'testrun', 'feff8l.log'))
        self.__testpaths()
        self.repotop    = getcwd()
        if not self.repotop.endswith('feff85exafs'):  self.repotop = realpath(join('..'))
        # the f85e shell script emulates the behavior of the monolithic Feff application
        self.eps5       = 0.00001
        self.eps4       = 0.0001
        self.eps3       = 0.001
        self.epsilon    = self.eps4
        self.epsfit     = self.eps3
        self.epserr     = 5.0 * self.epsfit
        self.firstshell = False
        self.fittest    = None
        if WRAPPER_AVAILABLE:
            self.sp = Feff8L_XAFSPath(_larch=self._larch)

    def __repr__(self):
        if not isdir(self.folder):
            return '<Feff85exafs Unit Test Group (empty)>'
        if self.folder is not None:
            return '<Feff85exafs Unit Test Group: %s>' % self.folder
        return '<Feff85exafs Unit Test Group (empty)>'


    @property
    def baseline(self):
        scf = 'withSCF' if self.doscf else 'noSCF'
        return realpath(join(self.path, 'baseline', scf))

    @property
    def bpaths(self):
        """
        Gather a list of feffNNNN.dat files from the baseline calculation

        """
        p   = list()
        owd = getcwd()
        try:
            chdir(self.baseline)
            feffoutput = glob.glob("*")
            for f in sorted(feffoutput):
                tosave = re.compile("feff\d+\.dat")
                if tosave.match(f):
                    p.append(f)
        finally:
            chdir(owd)
        return p

    def run(self):
        """
        Make a feff.inp from the mustache template, then run feff 8.5 in
        the testrun folder

        """
        if not isdir(self.path):
            print_error(self.folder + " is not one of the available tests")
            return False

        ## clean up an earlier run
        if isdir(self.testrun): rmtree(self.testrun)
        makedirs(self.testrun)

        ## mustache to feff.inp
        scf = 'without SCF'
        self.json = json.load(open(join(self.path, self.folder + '.json')))
        self.json['doscf']='* '
        if self.doscf:
            scf = 'with SCF'
            self.json['doscf']=''
        renderer = pystache.Renderer()
        with open(join(self.testrun,'feff.inp'), 'w') as inp:
            inp.write(renderer.render_path( join(self.path, self.folder + '.mustache'), # material/material.mustache
                                            self.json ))                                # material/material.json

        ## run feff with feffrunner
        owd = getcwd()
        try:
            chdir(self.testrun)
            if isfile(self.fefflog):
                unlink(self.fefflog)

            self.feffrunner=feffrunner(feffinp=join(self.testrun,'feff.inp'), verbose=self.verbose, repo=self.repotop, _larch=self._larch)
            self.feffrunner.run()
            for f in glob.glob("*"):
                if f.startswith('log'):
                    unlink(f)
            if self.verbose:
                print_warn("\nRan Feff85EXAFS on %s (%s)" % (self.folder, scf))
            self.__testpaths()
        finally:
            chdir(owd)
        self.feffran = True

    def run_opconsat(self):
        if not self.feffran:
            print_error("You need to run the rest of the feff calculation first.")
            return False
        ## run feff with feffrunner
        owd = getcwd()
        try:
            chdir(self.testrun)
            self.feffrunner=feffrunner(feffinp=join(self.testrun,'feff.inp'), verbose=self.verbose, repo=self.repotop, _larch=self._larch)
            copy(join(self.testrun, "..", "opconsat", "baseline", "exc.inp"), self.testrun)
            self.feffrunner.run('opconsat')
        finally:
            chdir(owd)
        if not isfile(join(self.testrun, "exc.dat")):
            print_error("Failed to run opconsat.")
            return False
        return True

    def __testpaths(self):
        """
        Gather a list of feffNNNN.dat files from the testrun

        """
        self.paths = list()
        if isdir(self.testrun):
            owd = getcwd()
            try:
                chdir(self.testrun)
                feffoutput = glob.glob("*")
                for f in sorted(feffoutput):
                    tosave = re.compile("feff\d+\.dat")
                    if tosave.match(f):
                        self.paths.append(f)
                self.feffran = True
            finally:
                chdir(owd)


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


    def __snarf_geometry(self, index=1):
        pathsdat = self.testrun+'/paths.dat'
        try:
            lines = open(pathsdat, 'r').readlines()
        except:
            print( 'Error reading file %s ' % pathsdat)
            return

        degmatch = re.compile("degeneracy")
        mode = 'searching'
        count = 1
        for line in lines:
            line = line.replace(',', ' ').strip()
            words  = line.split()
            if len(words) < 5:
                continue
            elif mode == 'searching':
                #print '> %s %s' % (line, str(degmatch.match(line)))
                if words[5] == 'degeneracy' and line.startswith(str(index)):
                    self.sp.index = index
                    nlegs = int(words[1])
                    self.sp.deg = float(words[2])
                    mode = 'parsing'
                    continue
            elif mode == 'parsing':
                #print '| %s %d %d' % (line, count, nlegs)
                if words[5] == 'degeneracy':
                    break
                elif count == nlegs:
                    break
                elif line.startswith('x'):
                    continue
                else:
                    words = line.split()
                    self.sp.add_scatterer(x=float(words[0]),
                                          y=float(words[1]),
                                          z=float(words[2]),
                                          ipot=int(words[3]))
                    count+=1

    def compare(self, nnnn=1, part='feff', use_wrapper=False, _larch=None):
        """
        Compare a feffNNNN.dat file from the testrun with the same file
        from the baseline calculation

            group.compare(N, part, use_wrapper)

        where N is the path index, part is one of
            feff    :      test the magnitude AND phase of F_eff (default)
            amp     :      test the total amplitude of the path
            phase   :      test the total phase of the path
            lam     :      test the mean free path
            caps    :      test the central atom phase shift
            red_fact:      test the reduction factor
            rep     :      test the real part of the complex wavenumber

        and use_wrapper is True if the test is to use the python interface to
        the feffpath library or False if the test is to use the feffrunner and
        the Feff executables.  (Currently, only feffrunner is working...)

        Sets self.rfactor (and self.rfactor_2).

        """
        if self._larch is None:
            raise Warning("cannot do path comparison -- larch broken?")

        if not self.feffran:
            print_error("You have not yet run the test Feff calculation")
            return

        if not self.available(nnnn):
            print_error("Path %d was not saved from the test Feff calculation" % nnnn)
            return

        how     = 'wrapper' if use_wrapper else 'executables'
        nnnndat = "feff%4.4d.dat" % nnnn

        blpath = feffpath(join(self.baseline, nnnndat), _larch=self._larch)
        if use_wrapper: # make the feffNNNN.dat file on the fly
            self.sp.phase_file = join(self.testrun, 'phase.pad')
            self.sp.nnnn=True
            self.sp.index=nnnn
            self.sp.verbose=self.verbose
            self.__snarf_geometry(nnnn)
            owd = getcwd()
            try:
                chdir(self.testrun)
                self.sp.calculate_xafs()
            finally:
                chdir(owd)
            n3nndat = "feff%4.4d.dat" % nnnn
            trpath = feffpath(join(self.testrun,  n3nndat), _larch=self._larch)
        else:                   # the feffNNNN.dat file was made by the monolithic feff run
            trpath = feffpath(join(self.testrun,  nnnndat), _larch=self._larch)

        if part=='feff':
            baseline_1 = getattr(blpath._feffdat, 'mag_feff') # + np.random.uniform(0,1,size=1)
            testrun_1  = getattr(trpath._feffdat, 'mag_feff')
            baseline_2 = np.sin(getattr(blpath._feffdat, 'pha_feff')) # + np.random.uniform(0,1,size=1)
            testrun_2  = np.sin(getattr(trpath._feffdat, 'pha_feff'))
            ylabel     = 'magnitude and phase'
            label      = 'magnitude'
        elif part=='amp':
            baseline_1 = getattr(blpath._feffdat, 'amp')
            testrun_1  = getattr(trpath._feffdat, 'amp')
            ylabel     = 'total amplitude'
            label      = 'amplitude'
        elif part=='phase':
            baseline_1 = np.sin(getattr(blpath._feffdat, 'pha'))
            testrun_1  = np.sin(getattr(trpath._feffdat, 'pha'))
            ylabel     = 'total phase shift'
            label      = 'phase'
        elif part=='lam':
            baseline_1 = getattr(blpath._feffdat, 'lam')
            testrun_1  = getattr(trpath._feffdat, 'lam')
            ylabel     = 'mean free path'
            label      = 'lambda'
        elif part=='phc':
            baseline_1 = getattr(blpath._feffdat, 'real_phc')
            testrun_1  = getattr(trpath._feffdat, 'real_phc')
            ylabel     = 'central atom phase shift'
            label      = 'phc'
        elif part=='red_fact':
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
            if self.verbose:
                print_error("Unknown choice of parts \"%s\"\nmust be one of (feff|amp|phase|lam|phc|red_fact|rep)\nusing feff" % part)
            part       = 'feff'
            baseline_1 = getattr(blpath._feffdat, 'mag_feff')
            testrun_1  = getattr(trpath._feffdat, 'mag_feff')
            baseline_2 = np.sin(getattr(blpath._feffdat, 'pha_feff'))
            testrun_2  = np.sin(getattr(trpath._feffdat, 'pha_feff'))
            ylabel     = 'magnitude and phase'
            label      = 'magnitude'


        self.rfactor_2 = 0
        self.rfactor = sum((baseline_1 - testrun_1)**2) / sum(baseline_1**2)
        if self.verbose:
            print_warn("\nComparing %s of %s (%s) (using %s)" % (label, nnnndat, "with SCF" if self.doscf else "without SCF", how))
            self.print_geometry(blpath)
            print_warn("%s:  %f v %f " % (blpath, sum(baseline_1**2) ,  sum(testrun_1**2)))

        if self.verbose:
            print label + " R-factor = " + test_text("%.9g" % self.rfactor, (self.rfactor < self.epsilon))
        if part=='feff':
            self.rfactor_2 = sum((baseline_2 - testrun_2)**2) / sum(baseline_2**2)
            if self.verbose:
                print(" -- ",self.epsilon, sum((baseline_2 - testrun_2)**2),  sum(baseline_2**2))
                print "phase R-factor = " + test_text("%.9g" % self.rfactor_2, self.rfactor_2 < self.epsilon)

        if self.verbose: print ""

        if self.doplot:
            _newplot(blpath._feffdat.k, baseline_1, _larch=self._larch, label=label+' of baseline',
                     xlabel='wavenumber $\AA^{-1}$', ylabel=ylabel, title=nnnndat, show_legend=True, legend_loc='best')
            _plot   (trpath._feffdat.k, testrun_1,  _larch=self._larch, label=label+' of test run')
            if part=='feff':
                _plot(blpath._feffdat.k, np.gradient(baseline_2), _larch=self._larch, label='grad(phase of baseline)')
                _plot(trpath._feffdat.k, np.gradient(testrun_2),  _larch=self._larch, label='grad(phase of test run)')

        if use_wrapper and WRAPPER_AVAILABLE:
            self.sp.reset()

        if part=='feff':
            return self.rfactor < self.epsilon and self.rfactor_2 < self.epsilon
        else:
            return self.rfactor < self.epsilon


    def print_geometry(self, path):
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
        Return the value of S02 calculated by the testrun or the baseline.
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
        Fetch the various numbers in the feffNNNN.dat header for both
        the baseline and testrun.  A table of these values is printed
        if the verbose flag is True.  Return True if all are equal.

        """
        if (not self.feffran):
            #if self.verbose: print colored("You have not yet made the test run of Feff", 'magenta', attrs=['bold'])
            raise Exception("You have not yet made the test run of Feff")
        nnnndat = "feff%4.4d.dat" % nnnn
        bl      = feffpath(join(self.baseline, nnnndat), _larch=self._larch)
        tr      = feffpath(join(self.testrun,  nnnndat), _larch=self._larch)

        termdict = {'edge':   'energy threshold relative to atomic value',
                    'gam_ch': 'core level energy width',
                    'kf':     'k value at Fermi level',
                    'mu':     'Fermi level in eV',
                    'rs_int': 'interstitial radius',
                    'vint':   'interstitial potential'}
        if self.verbose:
            print "%-8s %-42s   %10s  %10s  %s" % ('key', 'description', 'baseline', 'testrun', 'same?')
            print "-----------------------------------------------------------------------------------"
        ok = True
        for key in termdict:
            same = abs(getattr(bl._feffdat, key) - getattr(tr._feffdat, key))/getattr(bl._feffdat, key) < self.epsilon
            if self.verbose: print "%-6s   %-42s : %10s  %10s  %s" % (key, termdict[key], getattr(bl._feffdat, key),
                                                                      getattr(tr._feffdat, key), same)
            ok = ok and same
        return ok



    def fit(self):
        """
        Perform a canned fit using the baseline and testrun Feff
        calculations.

        Sets self.blfit and self.trfit, containing the feffit groups
        for the fits suing the baseline and test run Feff
        calculations, respectively.

        """
        sys.path.append(self.path)
        module = importlib.import_module(self.folder, package=None)

        print "\tfitting %s using %s" % (self.folder, 'baseline')
        self.blfit = module.do_fit(self, 'baseline')
        print "\tfitting %s using %s" % (self.folder, 'testrun')
        self.trfit = module.do_fit(self, 'testrun')


    def fitcompare(self):
        """
        Perform a canned fit using a sequence of Feff calculations

        This is not used for unit testing feff85exafs.  It is a
        generalization of the fit() method intended for use with a
        generic testing framework, e.g. https://github.com/bruceravel/SCFtests

        Uses self.fittest, the name of a folder containing a sequence of
        Feff calculations.

        Sets self.models, a list of subfolders containing the identifying
        strings of the sequence of calculations.

        Returns a bunch of self.<id>, where <id> are feffit groups
        containing the fits in the sequence.
        """
        sys.path.append(self.folder)
        module = importlib.import_module(self.folder, package=None)
        save = self.doscf

        if self.fittest is None:
            raise Exception("You must select a test folder.")

        self.models = []
        for d in sorted(listdir(join(self.folder, self.fittest))):
            if isdir(join(self.folder, self.fittest, d)):
                print '>>>>>>>>> fitting with model: %s' % d
                self.models.append(d)
                this = module.do_fit(self, d, firstshell=self.firstshell, fittest=self.fittest)
                setattr(self, d, this)


    def clean(self):
        """
        Remove testrun folder
        """
        if isdir(self.testrun): rmtree(self.testrun)
        self.feffran = False


######################################################################

def ut(folder=None, _larch=None, **kws):
    """
    Make a Feff85exafsUnitTestGroup group given a folder containing a
    baseline calculation

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
    """_i_mport, _r_un, and _f_it

    """
    utobj=ut(folder)
    utobj.run()
    utobj.fit()
    return utobj

def registerLarchPlugin(): # must have a function with this name!
    return ('f85ut', { 'ut': ut, 'ir': ir, 'irc': irc, 'irf': irf })
