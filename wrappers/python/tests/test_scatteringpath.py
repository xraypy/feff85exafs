from scatteringpath import scatpath

a=scatpath()
epsilon=1e-4

def test_phpad():
    yield check_phpad

def test_degen():
    yield check_degen, 0
    yield check_degen, 48

def test_index():
    yield check_index, 0
    yield check_index, 1

def test_iorder():
    yield check_iorder, 0
    yield check_iorder, 1


def test_nnnn():
    yield check_nnnn, True
    yield check_nnnn, False

def test_verbose():
    yield check_verbose, True
    yield check_verbose, False


def test_evec():
    yield check_evec

def test_xivec():
    yield check_xivec

def test_atom():
    a.clear()
    yield check_index, 0
    a.degen = 48
    a.index = 4
    a.phpad = "../fortran/phase.pad"
    a.nnnn  = False
    a.atom(0, 0, -3.61, 1)
    yield check_atom, 2
    a.atom(-1.805, 0, -1.805, 1)
    yield check_atom, 3
    a.make()
    yield check_ri, 'ri'
    yield check_beta, 'beta'
    yield check_param, 'edge'
    yield check_param, 'gam_ch'
    yield check_param, 'kf'
    yield check_param, 'mu'
    yield check_param, 'rnorman'
    yield check_param, 'rs_int'
    yield check_param, 'vint'
    yield check_param, 'exch'

def test_atom_errors():
    a.clear();
    a.atom(0, 0, -3.61, -1)     # ipot negative
    yield check_error, 1

    a.clear();
    a.atom(0, 0, -3.61, 9)      # ipot too big
    yield check_error, 2

    a.clear();
    a.atom(0, 0, -3.61, 1)
    a.atom(0, 0, -3.71, 1)
    yield check_error, 4        # atoms too close together

    a.clear();
    a.atom(0, 0, 0, 1)
    a.atom(0, 0, -3.61, 1)
    a.make()
    yield check_error, 1        # first atom absorber

    a.clear();
    a.atom(0, 0, -3.61, 1)
    a.atom(0, 0, 0, 1)
    a.make()
    yield check_error, 2        # last atom absorber

    a.clear();
    try:
        a.degen=-12
    except ValueError:
        yield check_pythonerror, 'degen'        # degeneracy negative
    else:
        yield check_wrong_pythonerror, 'degen'

    try:
        a.index=40000
    except ValueError:
        yield check_pythonerror, 'index'        # bad index
    else:
        yield check_wrong_pythonerror, 'index'

    try:
        a.elpty=-0.5
    except ValueError:
        yield check_pythonerror, 'elpty'        # bad elpty
    else:
        yield check_wrong_pythonerror, 'elpty'

    try:
        a.iorder=-1
    except ValueError:
        yield check_pythonerror, 'iorder'        # bad elpty
    else:
        yield check_wrong_pythonerror, 'iorder'

    try:
        a.phpad='foo.bar'
    except ValueError:
        yield check_pythonerror, 'phpad'        # bad phpad
    else:
        yield check_wrong_pythonerror, 'phpad'


################################################################################


def check_phpad():
    pb = "../fortran/phase.pad";
    a.phpad = pb
    assert a.wrapper.phpad == pb

def check_degen(degen):
    if degen:
        a.degen = degen
        assert a.wrapper.degen == degen
    else:
        assert a.degen == 1

def check_index(ind):
    if ind:
        a.index = ind
        assert a.wrapper.index == ind
    else:
        assert a.index == 9999

def check_iorder(iorder):
    if iorder:
        a.iorder = iorder
        assert a.wrapper.iorder == iorder
    else:
        assert a.iorder == 2

def check_nnnn(nnnn):
    a.nnnn = nnnn
    assert a.wrapper.nnnn == nnnn

def check_verbose(verbose):
    a.verbose = verbose
    assert a.wrapper.verbose == verbose

def check_evec():
    a.evec = [0,0,1]
    assert a.evec[0] == 0 and a.evec[1] == 0 and a.evec[2] == 1

def check_xivec():
    a.evec = [1,1,0]
    assert a.evec[0] == 1 and a.evec[1] == 1 and a.evec[2] == 0

def check_atom(nleg):
    assert a.nleg == nleg

def check_ri(which):
    assert abs(a.ri[0] - 3.61) < epsilon and abs(a.ri[1] - 2.55266) < epsilon and abs(a.ri[2] - 2.55266) < epsilon

def check_beta(which):
    assert abs(a.beta[0] - 135) < epsilon and abs(a.beta[1] - 90) < epsilon and abs(a.beta[2] - 135) < epsilon

def check_error(err):
    assert a.errorcode == err

def check_pythonerror(which):
    assert True
def check_wrong_pythonerror(which):
    assert False

def check_param(which):
    if (which == 'edge'):
        assert abs(a.edge + 3.82035) < epsilon
    elif (which == 'gam_ch'):
        assert abs(a.gam_ch - 1.72919) < epsilon
    elif (which == 'kf'):
        assert abs(a.kf - 1.824) < epsilon
    elif (which == 'mu'):
        assert abs(a.mu + 3.82035) < epsilon
    elif (which == 'rnorman'):
        assert abs(a.rnorman - 2.63173) < epsilon
    elif (which == 'rs_int'):
        assert abs(a.rs_int - 1.98947) < epsilon
    elif (which == 'vint'):
        assert abs(a.vint + 16.48140) < epsilon
    elif (which == 'exch'):
        assert a.exch == 'H-L exch'
    else:
        assert False;
