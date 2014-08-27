## feff85exafs unit testing system using larch
## see HEADERS/license.h for feff's license information

import larch
from f85ut import ut

from   os.path   import isfile, isdir, join

fut           = ut('NiO')
fut.verbose   = False
fut.doplot    = False
fut.doscf     = False

skip_data_test = False

def test_columns():
    if not fut.feffran: fut.run()
    for part in ['lambda', 'caps', 'redfact', 'rep']:
        yield check_columns, part

def check_columns(part):
    this = fut.compare(1, part=part)
    assert this, "comparison of %s for path 1 in NiO" % part

################################################################################

def test_feff():
    if not fut.feffran: fut.run()
    for path in fut.paths:
        index = int(path[4:8])
        for part in ['feff', 'amp', 'phase']:
            yield check_feff, index, part

def check_feff(index, part):
    this = fut.compare(index, part=part)
    assert this, "comparison of %s for path %d in NiO" % (part, index)


################################################################################

def test_terms():
    assert fut.feffterms(), "some feff terms calculated incorrectly for NiO"

def test_s02():
    assert fut.s02() == fut.s02('baseline'), "s02 calculated incorrectly for NiO"

################################################################################

def test_fit():
    if not fut.feffran: fut.run()
    if skip_data_test:
        yield check_null, "skipping data test for NiO"
    elif isfile(join(fut.path, fut.folder+'.py')):
        fut.fit()
        for p in fut.blfit.params.covar_vars:
            yield check_param, p, 'value'
            yield check_param, p, 'stderr'
        for p in ['chi_reduced', 'chi_square', 'rfactor']:
            yield check_stat, p
    else:
        yield check_null, "no data tests for NiO"

def check_null(msg):
    assert True, msg


def check_param(param, part):
    bl=getattr(getattr(fut.blfit.params, param), part)
    tr=getattr(getattr(fut.trfit.params, param), part)
    assert abs(bl-tr) < fut.epsilon, "%s of fitting parameter %s evaluated inconsistently for NiO" % (part, param)

def check_stat(param):
    bl=getattr(fut.blfit.params, param)
    tr=getattr(fut.trfit.params, param)
    assert abs(bl-tr) < fut.epsilon, "statistic %s evaluated inconsistently for NiO" % param

################################################################################

def test_clean():
    fut.clean()
    assert not isdir(fut.testrun), "successfully cleaned up NiO calculation"
