## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

import os
import sys
from os.path import isfile, isdir, join, exists
import re

from math import isnan

import nose
import pytest

import larch
from larch_plugins.xafs.feffdat import feffpath
from larch_plugins.io import read_ascii

from f85ut import ut

ALL_FOLDERS = ('Copper', 'NiO', 'UO2', 'Zircon', 'ferrocene',
               'bromoadamantane', 'LCO-para', 'LCO-perp')

class Feff8Test:
    def __init__(self, folder, verbose=True, doplot=False, doscf=None):
        if folder not in ALL_FOLDERS:
            raise ValueError("Must choose one of the available tests:\n  %s\n" % repr(ALL_FOLDERS))

        if doscf is None:
            envar = os.getenv('FEFF_TEST_SCF', 'False')
            doscf = envar.lower() in ("yes", "true", "on", "y", "t", "1")

        self.folder = folder
        self.test = ut(folder)
        self.test.verbose = verbose
        self.test.doplot = doplot
        self.test.doscf = doscf

    def test_feffrun(self):
        "run feff"
        self.test.run()
        scf = 'without SCF'
        if self.test.doscf:
            scf = 'with SCF'
        # feff8l.log shouldn't be more than a couple thousand lines long
        # (ferrocene w/SCF is 1096)
        with open(join(self.test.testrun, 'feff8l.log'), 'r') as log:
            lines = log.readlines()
            endlines = '_'.join(lines[-5:])
            m = 'Done with module 6:' in endlines

        assert m and self.test.available(1), "feff run on %s (%s) not successful" % (self.folder, scf)


    def test_feff(self):
        """check F_eff for each path"""
        assert self.test.feffran, "no test results found for %s" % self.folder
        for path in self.test.paths:
            index = int(path[4:8])
            for part in ('feff', 'amp', 'phase',
                         'lam', 'phc', 'red_fact', 'rep'):
                result = self.test.compare(index, part=part)
                assert result, "comparison of %s for path %d in %s" % (part, index, self.folder)

    def test_terms(self):
        "check various quantities computed by feff"
        assert self.test.feffran, "no test results found for %s" % self.folder
        for npath in range(1, 10):
            basefile = join(self.test.baseline, 'feff%4.4i.dat' % npath)
            if not exists(basefile):
                continue
            bl = feffpath(join(self.test.baseline, 'feff%4.4i.dat' % npath), _larch=self.test._larch)
            tr = feffpath(join(self.test.testrun,  'feff%4.4i.dat' % npath), _larch=self.test._larch)
            for term in ('edge', 'gam_ch', 'kf', 'mu', 'rs_int', 'vint'):
                blval = getattr(bl._feffdat, term)
                trval = getattr(tr._feffdat, term)
                tdiff = blval - trval
                assert abs(tdiff) < 2.0e-4, "feff term %s not close enough for %s" % (term, self.folder)

        for radius in ('muffintin', 'norman'):
            bl = self.test.radii('baseline', radius)
            tr = self.test.radii('testrun',  radius)
            assert bl == tr, "list of %s radii are different for %s" % (radius, self.folder)
        assert self.test.s02() == self.test.s02('baseline'), "s02 calculated incorrectly for %s" % self.folder

    def test_clean(self):
        """remove the test run"""
        self.test.clean()
        assert not isdir(self.test.testrun), "clean up %s calculation not successful" % self.folder

    def test_columns_wrapper(self):
        assert self.test.feffran, "no test results found for %s" % self.folder
        msg = "comparison of %s for path %i in %s using wrapper"

        for path in self.test.paths:
            index = int(path[4:8])
            for part in ('feff', 'amp', 'phase', 'lam', 'phc', 'red_fact', 'rep'):
                result = self.test.compare(index, part=part, use_wrapper=True)
                assert result, msg % (part, index, self.folder)


    def test_fit(self):
        """check various fitting and statistical parameters from a canned fit"""
        assert self.test.feffran, "no test results found for %s" % self.folder

        test_script = join(self.test.path, self.test.folder+'.py')
        stat_msg = "statistic %s evaluated inconsistently for %s (%.5f %.5f)"
        param_msg = "%s of fitting parameter %s evaluated inconsistently for %s (%.5f %.5f)"

        if isfile(test_script):
            self.test.fit()
            eps = self.test.epsfit
            for stat in ('chi_reduced', 'chi_square', 'rfactor'):
                bl = getattr(self.test.blfit, stat)
                tr = getattr(self.test.trfit, stat)
                close = abs((bl-tr)/bl) < eps
                assert close, stat_msg % (stat, self.folder, bl, tr)

            for par in self.test.blfit.var_names:
                for part, eps in (('value',  self.test.epsfit),
                                  ('stderr', self.test.epserr)):
                    bl = getattr(self.test.blfit.params[par], part)
                    tr = getattr(self.test.trfit.params[par], part)
                    close = abs((bl-tr)/bl) < eps
                    assert close, param_msg % (part, par, self.folder, bl, tr)


    @pytest.mark.skip(reason="not testing opconsat")
    def test_feffrun_opconsat(self):
        oca_ok = self.test.run_opconsat()
        assert oca_ok, "opconsat run on %s not successful" % folder

    @pytest.mark.skip(reason="not testing opconsat")
    def test_opconsat(self):
        check_skip("not checking opconstat")
        self.test.run_opconsat()
        for f in self.tests:
            pass # check_opconsat(f)


################################################################################


def check_opconsat(folder):
    orig = read_ascii(join(tests[folder].testrun, "..", "opconsat", "baseline", "exc.dat"), labels='a b c d', _larch=tests[folder]._larch)
    if not isfile(join(tests[folder].testrun, "exc.dat")):
        assert False, 'Failed to run opconsat'
    new  = read_ascii(join(tests[folder].testrun, "exc.dat"), labels='a b c d', _larch=tests[folder]._larch)
    #show(orig, _larch=self._larch)
    #show(new, _larch=self._larch)
    rf1  = sum((orig.a - new.a)**2) / sum(orig.a**2)
    rf2  = sum((orig.b - new.b)**2) / sum(orig.b**2)
    rf3  = sum((orig.c - new.c)**2) / sum(orig.c**2)
    rf4  = sum((orig.d - new.d)**2) / sum(orig.d**2)
    which = []
    failed = False
    if isnan(rf1):
        which.append('1')
        failed = True
    if isnan(rf2):
        which.append('2')
        failed = True
    if isnan(rf3):
        which.append('3')
        failed = True
    if isnan(rf4):
        which.append('4')
        failed = True
    if failed:
        assert False, "columns of exc.dat are not a number: %s" % ', '.join(which)

    which = []
    failed = False
    if rf1 < tests[folder].epsfit:
        which.append('1')
        failed = True
    if rf2 < tests[folder].epsfit:
        which.append('2')
        failed = True
    if rf3 < tests[folder].epsfit:
        which.append('3')
        failed = True
    if rf4 < tests[folder].epsfit:
        which.append('4')
        failed = True
    assert failed, "columns of exc.dat are not equal to baseline: %s" % ', '.join(which)

if __name__ == '__main__':
    TEST_FOLDERS = ('Copper', 'NiO', 'Zircon', 'ferrocene', 'LCO-para', 'LCO-perp')
    # TEST_FOLDERS = ALL_FOLDERS
    for folder in TEST_FOLDERS:
        t =  Feff8Test(folder)
        t.test_feffrun()
        t.test_feff()
        t.test_terms()
        t.test_columns_wrapper()
        t.test_fit()
        t.test_clean()
