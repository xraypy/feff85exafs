## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

import os
import sys
from os.path import isfile, isdir, join, exists

from F85_Tester import Feff85exafsUnitTestGroup

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
        self.test = Feff85exafsUnitTestGroup(folder, verbose=verbose,
                                             doplot=doplot, doscf=doscf)


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

    def test_feffresults(self):
        """check F_eff for each path"""
        self.test.compare_path_results()

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



################################################################################

if __name__ == '__main__':
    TEST_FOLDERS = ('Copper', 'NiO', 'Zircon', 'ferrocene', 'LCO-para', 'LCO-perp')
    scfopts = (False, True)
    for folder in TEST_FOLDERS:
        for doscf in scfopts:
            t =  Feff8Test(folder, doscf=doscf)
            t.test_feffrun()
            t.test_feffresults()
            t.test_columns_wrapper()
            t.test_fit()
            t.test_clean()
