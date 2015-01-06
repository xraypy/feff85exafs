## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

import larch
from f85ut import ut
larch.use_plugin_path('xafs')
from feffdat import feffpath

import re
from os import getenv
from os.path import isfile, isdir, join

#folders = ('Copper',)
#folders = ('Copper', 'NiO', 'UO2', 'Zircon', 'ferrocene', 'bromoadamantane')
folders = ('Copper', 'NiO', 'UO2', 'Zircon', 'ferrocene', 'bromoadamantane', 'LCO-para', 'LCO-perp')
tests   = dict()
doscf   = getenv('FEFF_TEST_SCF', 'False')

for f in folders:
    tests[f]         = ut(f)
    tests[f].verbose = False
    tests[f].doplot  = False
    tests[f].doscf   = doscf.lower() in ("yes", "true", "on", "y", "t", "1")


## run feff
def test_feffrun():
    for f in folders:
        yield check_feffrun, f

## check the feffNNNN.dat columns that are the same for all paths
def test_columns():
    for f in folders:
        for part in ('lambda', 'caps', 'redfact', 'rep'):
            yield check_columns, f, part

def test_columns_wrapper():
    if not tests[folders[0]].wrapper_available:
        yield check_false, "wrapper unavailable"
        return;
    for f in folders:
        for part in ('lambda', 'caps', 'redfact', 'rep'):
            yield check_columns_wrapper, f, part

## check F_eff for each path
def test_feff():
    for f in folders:
        for path in tests[f].paths:
            index = int(path[4:8])
            for part in ('feff', 'amp', 'phase'):
                yield check_feff, f, index, part

def test_feff_wrapper():
    if not tests[folders[0]].wrapper_available:
        yield check_false, "wrapper unavailable"
        return;
    for f in folders:
        for path in tests[f].paths:
            index = int(path[4:8])
            for part in ('feff', 'amp', 'phase'):
                if f == ('LCO-perp'):
                    tests[f].sp.evec  = (1,1,0)
                    tests[f].sp.xivec = (0,0,1)
                    tests[f].sp.elpty = 1
                elif f == ('LCO-para'):
                    tests[f].sp.evec  = (0,0,1)
                yield check_feff_wrapper, f, index, part

## check norman and muffin tin radii of the ipots from feff
def test_radii():
    for f in folders:
        yield check_radii, f, 'muffintin'
        yield check_radii, f, 'norman'


## check various quantities computed by feff
def test_terms():
    for f in folders:
        for term in ('edge', 'gam_ch', 'kf', 'mu', 'rs_int', 'vint'):
            yield check_feffterms, f, term

## check that S02 was computed correctly
def test_s02():
    for f in folders:
        yield check_s02, f

## check various fitting and statistical parameters from a canned fit
def test_fit():
    for f in folders:
        #if not tests[f].feffran: tests[f].run()
        if isfile(join(tests[f].path, tests[f].folder+'.skip')):
            yield check_true, "skipping data test for %s" % f
        elif isfile(join(tests[f].path, tests[f].folder+'.py')):
            tests[f].fit()
            for p in ('chi_reduced', 'chi_square', 'rfactor'):
                yield check_stat, f, p
            if (hasattr(tests[f].blfit.params, 'covar_vars')):
                for p in tests[f].blfit.params.covar_vars:
                    yield check_param, f, p, 'value'
                    yield check_param, f, p, 'stderr'
            else:
                yield check_false, "fit could not evaluate uncertainties for %s" % f
        else:
            yield check_true, "no data tests for %s" % f

## remove the test run
def test_clean():
    for f in folders:
        tests[f].clean()
        yield check_clean, f


################################################################################


def check_feffrun(folder):
    tests[folder].run()
    scf = 'without SCF'
    if tests[folder].doscf:
        scf = 'with SCF'
    with open(join(tests[folder].testrun,'f85e.log'), 'r') as log:
        lines = log.readlines() # f85e.log shouldn't be more than a couple thousand lines long (ferrocene w/SCF is 1096)
        m = re.search('Done with module 6:', lines[-1])
    assert m and tests[folder].available(1), "feff run on %s (%s) not successful" % (folder, scf)

def check_columns(folder, part):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    this = tests[folder].compare(1, part=part)
    assert this, "comparison of %s for path 1 in %s" % (part, folder)
def check_columns_wrapper(folder, part):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    this = tests[folder].compare(1, part=part, use_wrapper=True)
    assert this, "comparison of %s for path 1 in %s using wrapper" % (part, folder)

def check_feff(folder, index, part):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    this = tests[folder].compare(index, part=part)
    assert this, "comparison of %s for path %d in %s" % (part, index, folder)
def check_feff_wrapper(folder, index, part):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    this = tests[folder].compare(index, part=part, use_wrapper=True)
    assert this, "comparison of %s for path %d in %s using wrapper" % (part, index, folder)

def check_radii(folder, radius):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    bl = tests[folder].radii('baseline', radius)
    tr = tests[folder].radii('testrun',  radius)
    assert bl == tr, "list of %s radii are different for %s" % (radius, folder)

def check_feffterms(folder, term):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    bl = feffpath(join(tests[folder].baseline, 'feff0001.dat'))
    tr = feffpath(join(tests[folder].testrun,  'feff0001.dat'))
    assert abs(getattr(bl._feffdat, term) - getattr(tr._feffdat, term)) < tests[folder].eps4, "feff term %s calculated incorrectly for %s" % (term, folder)

def check_s02(folder):
    if not tests[folder].feffran: assert False, "failed to find results of feff calculation for %s" % folder
    assert tests[folder].s02() == tests[folder].s02('baseline'), "s02 calculated incorrectly for %s" % folder
    
def check_true(msg):
    assert True, msg

def check_false(msg):
    assert False, msg

def check_param(folder, param, part):
    bl=getattr(getattr(tests[folder].blfit.params, param), part)
    tr=getattr(getattr(tests[folder].trfit.params, param), part)
    assert abs(bl-tr) < tests[folder].epsfit, "%s of fitting parameter %s evaluated inconsistently for %s" % (part, param, folder)

## this tests for < 0.1% deviation in statistic value, scaled to the size of the statistic
def check_stat(folder, param):
    bl=getattr(tests[folder].blfit.params, param)
    tr=getattr(tests[folder].trfit.params, param)
    assert (abs(bl-tr)/bl) < tests[folder].epsfit, "statistic %s evaluated inconsistently for %s (%.5f %.5f)" % (param, folder, bl, tr)

def check_clean(folder):
    assert not isdir(tests[folder].testrun), "clean up %s calculation not successful" % folder
