import numpy as np
import ctypes
from ctypes import (POINTER, pointer, c_int, c_long, c_char, c_char_p, c_double,
                    c_bool, c_byte, c_void_p, Structure)

import six

from matplotlib import pylab
FLIB = ctypes.cdll.LoadLibrary('../../local_install/lib/libfeff8lpath.dylib')
print("FLIB: ", FLIB)

FEFF_maxpts = 150  # nex
FEFF_maxpot = 11   # nphx
FEFF_maxleg = 9    # legtot
BOHR = 0.5291772490

string_attrs = ('exch', 'version')

def Py2tostr(val):
    return str(val)

def Py2tostrlist(address, nitems):
    return [str(i) for i in (nitems*c_char_p).from_address(address)]

def Py3tostr(val):
    if isinstance(val, str):
        return val
    if isinstance(val, bytes):
        return str(val, 'utf-8')
    return str(val)

def Py3tostrlist(address, nitems):
    return [str(i, 'latin_1') for i in (nitems*c_char_p).from_address(address)]

tostr  = Py2tostr
tostrlist = Py2tostrlist
if six.PY3:
    tostr = Py3tostr
    tostrlist = Py3tostrlist


def with_phase_file(fcn):
    """decorator to ensure that the wrapped function either
    has a non-None 'phase_file' argument or that that
    self.phase_file is not None
    """
    errmsg = "function '%s' needs a non-None phase_file"
    def wrapper(*args, **keywords):
        "needs phase_file"
        phase_file = keywords.get('phase_file', None)
        if phase_file is None:
            phase_file = getattr(args[0], 'phase_file', None)
            if phase_file is None:
                raise AttributeError(errmsg % fcn.__name__)
        else:
            setattr(args[0], 'phase_file', phase_file)
        # raise Warning(errmsg % fcn.__name__)
        return fcn(*args, **keywords)
    wrapper.__doc__ = fcn.__doc__
    wrapper.__name__ = fcn.__name__
    wrapper.__filename__ = fcn.__code__.co_filename
    wrapper.__dict__.update(fcn.__dict__)
    return wrapper

class Feff8PathStruct(Structure):
    _fields_  = [
        # INPUTS:
        ('phpad', c_char_p),     # path to phase.pad file
        ('index', c_int),         # path index
        ('nleg',  c_int),         # number of legs in path
        ('degen', c_double),      # path degeneracy
        ('rat',   c_void_p),      # path geometry (array)
        ('ipot',  c_void_p),      # path ipotentials (array)
        ('iorder', c_int),        # order of approx in GENFMT (default=2)
        ('nnnn',   c_byte),       # flag to write feffNNNN.dat
        ('xdi',    c_byte),       # flag to write feffNNNN.xdi
        ('verbose',  c_byte),     # flag to write screen messages
        ('ipol',  c_byte),        # flag for polarized calc
        ('evec',  c_void_p),      # polarization vector
        ('elpty',  c_double),     # polarization ellipticity
        ('xivec',  c_void_p),   # X-ray propagation direction

        # OUTPUTS:
        ('edge', c_double),       # energy threshold relative to atomic value
        ('gam_ch', c_double),     # core level energy width
        ('kf', c_double),         # k value at Fermi level
        ('mu', c_double),         # Fermi level in eV
        ('rnorman', c_double),    # Norman radiu
        ('rs_int', c_double),     # interstitial radius
        ('vint',  c_double),      # interstitial potential
        ('exch',  c_char_p),      # description of exchange model
        ('version', c_char_p),    # Feff version
        ('iz', c_void_p),         # geom: atomic number for atoms in path (array)
        ('ri', c_void_p),         # geom: leg lengths (array)
        ('beta', c_void_p),     # geom: beta angles (array)
        ('eta', c_void_p),      # geom: eta angles  (array)
        ('reff', c_double),       # geom: half path length
        ('ne',  c_int),           # number of energy points
        ('k',  c_void_p),       # k grid (array)
        ('real_phc', c_void_p), # central atom phase shifts (array of k)
        ('mag_feff', c_void_p), # magnitude of Feff (array of k)
        ('pha_feff', c_void_p), # phase of Feff (array of k)
        ('red_fact', c_void_p), # reduction factor (array of k)
        ('lam', c_void_p),      # mean free path (array of k)
        ('rep', c_void_p),      # real part of complex wavenumber (array of k)
        ('errorcode', c_int),     # error code
        ('errormessage', c_char_p), # error string
        ]


class ScatteringPath(object):
    """A Scatering Path for calculating a XAFS signal with Feff

    A calculation requires a Potentials and Phase Shift calculation
    in PAD format from Feff85, and a list of scattering paths

    Usage:
    ------
       # create path
       path  = ScatteringPath(phase_file='phase.pad')

       # list 'ipot' and labels for absorber, scattererers
       path.list_scatterers()

       # set coords for absorbing atom
       path.set_absorber(x=0., y=0., z=0.)

       # add scattering atom
       path.add_scatterer(x=1.5, y=1.5, z=1.5, ipot=1)

       # calculate basic (unaltered) XAFS contributions
       path.calcuate_xafs()

    """
    def __init__(self, phase_file=None):
        self.phase_file = phase_file
        self.clear()

    def clear(self):
        """reset all path data"""

        self.index   = 9999
        self.degen   = 1.
        self.nnnn_out = False
        self.xdi_out = False
        self.json_out = False
        self.verbose  = False
        self.ipol   = 0
        self.ellip  = 0.
        self.nepts  = 0
        self.genfmt_order = 2
        self.genfmt_vers  = ""
        self.exch_label   = ""
        self.rs     = 0.
        self.vint   = 0.
        self.xmu    = 0.
        self.edge   = 0.
        self.kf     = 0.
        self.rnorman = 0.
        self.gamach = 0.
        self.nepts  = FEFF_maxpts

        dargs = dict(dtype=np.float64, order='F')
        largs = dict(dtype=np.int32, order='F')

        self.evec   = np.zeros(3, **dargs)
        self.xivec  = np.zeros(3, **dargs)
        self.ipot   = np.zeros(1+FEFF_maxleg, **largs)
        self.beta   = np.zeros(1+FEFF_maxleg, **dargs)
        self.eta    = np.zeros(2+FEFF_maxleg, **dargs)
        self.ri     = np.zeros(FEFF_maxleg, **dargs)
        self.rat    = np.zeros((3, 2+FEFF_maxleg), **dargs)
        self.iz     = np.zeros(1+FEFF_maxpot, **largs)
        self.kfeff  = np.zeros(FEFF_maxpts, **dargs)
        self.real_phc = np.zeros(FEFF_maxpts, **dargs)
        self.mag_feff = np.zeros(FEFF_maxpts, **dargs)
        self.pha_feff = np.zeros(FEFF_maxpts, **dargs)
        self.red_fact = np.zeros(FEFF_maxpts, **dargs)
        self.lam      = np.zeros(FEFF_maxpts, **dargs)
        self.rep      = np.zeros(FEFF_maxpts, **dargs)
        self.nleg = 1

    @with_phase_file
    def list_scatterers(self, phase_file=None):
        """list Feff Potentials atoms ('ipots') fo phase file"""
        atoms = []
        with open(self.phase_file,'r') as fh:
            line1_words = fh.readline().strip().split()
            text = fh.readlines()
        nphases = int(line1_words[4])
        for line in text[4:]:
            if line.startswith('$'): continue
            words = line.split()
            atoms.append((int(words[1]), words[2]))
            if len(atoms) > nphases:
                break
        out = ["# Potential   Z   Symbol"]
        for ipot, atom in enumerate(atoms):
            out.append("    %2i      %3i     %s" % (ipot, atom[0], atom[1]))
        return "\n".join(out)

    @with_phase_file
    def set_absorber(self, x=0., y=0., z=0., phase_file=None):
        """set coordinates for absorbing atom ('ipot'=0)"""
        self.rat[0, 0] = x
        self.rat[1, 0] = y
        self.rat[2, 0] = z

    @with_phase_file
    def add_scatterer(self, x=0., y=0., z=0., ipot=1, phase_file=None):
        self.rat[0, self.nleg] = x
        self.rat[1, self.nleg] = y
        self.rat[2, self.nleg] = z
        self.ipot[self.nleg] = ipot
        self.nleg += 1
        # set final atom coords to same as absorber
        self.rat[0, self.nleg] = self.rat[0, 0]
        self.rat[1, self.nleg] = self.rat[1, 0]
        self.rat[2, self.nleg] = self.rat[2, 0]
        self.ipot[self.nleg]   = self.ipot[0]

    @with_phase_file
    def calculate_xafs(self, phase_file=None):

        class args: pass

        # strings / char*.  Note fixed length to match Fortran
        args.phase_file     = (self.phase_file + ' '*256)[:256]
        args.exch_label     = ' '*8
        args.genfmt_version = ' '*30

        # integers, including booleans
        for attr in ('index', 'nleg', 'genfmt_order', 'ipol', 'nnnn_out',
                     'json_out', 'verbose', 'nepts'):
            setattr(args, attr, pointer(c_long(int(getattr(self, attr)))))

        # doubles
        for attr in ('degen', 'rs', 'vint', 'xmu', 'edge', 'kf', 'rnorman',
                     'gamach', 'ellip'):
            setattr(args, attr, pointer(c_double(getattr(self, attr))))

        # integer arrays
        for attr in ('ipot', 'iz'):
            arr = getattr(self, attr)
            cdata = arr.ctypes.data_as(POINTER(arr.size*c_int))
            setattr(args, attr, cdata)

        # double arrays
        self.rat = self.rat/BOHR
        for attr in ('evec', 'xivec', 'rat', 'ri', 'beta', 'eta',
                     'kfeff', 'real_phc', 'mag_feff', 'pha_feff',
                     'red_fact', 'lam', 'rep'):
            arr = getattr(self, attr)
            cdata = arr.ctypes.data_as(POINTER(arr.size*c_double))
            setattr(args, attr, cdata)

        print(" --> MAKE ONE PATH")
        x = FLIB.calc_onepath(args.phase_file, args.index, args.nleg,
                              args.degen, args.genfmt_order, args.exch_label,
                              args.rs, args.vint, args.xmu, args.edge, args.kf,
                              args.rnorman, args.gamach, args.genfmt_version,
                              args.ipot, args.rat, args.iz, args.ipol,
                              args.evec, args.ellip, args.xivec, args.nnnn_out,
                              args.json_out, args.verbose, args.ri, args.beta,
                              args.eta, args.nepts, args.kfeff, args.real_phc,
                              args.mag_feff, args.pha_feff, args.red_fact,
                              args.lam, args.rep)

        self.exch_label   = args.exch_label.strip()
        self.genfmt_version = args.genfmt_version.strip()

        for attr in ('index', 'nleg', 'genfmt_order', 'degen', 'rs',
                     'vint', 'xmu', 'edge', 'kf', 'rnorman', 'gamach',
                     'ipol', 'ellip', 'nnnn_out', 'json_out', 'verbose',
                     'nepts'):
            setattr(self, attr, getattr(args, attr).contents.value)

        for attr in ('ipot', 'evec', 'xivec', 'beta', 'eta', 'ri', 'rat',
                    'iz', 'kfeff', 'mag_feff', 'pha_feff', 'red_fact',
                    'lam', 'rep'):
            cdata = getattr(args, attr).contents[:]
            setattr(self, attr, np.array(cdata))

        # some data needs recasting, reformatting
        self.nnnn_out = bool(self.nnnn_out)
        self.json_out = bool(self.json_out)
        self.verbose  = bool(self.verbose)
        self.rat = self.rat.reshape((2+FEFF_maxleg, 3)).transpose()*BOHR


if __name__ == '__main__':
    path = ScatteringPath(phase_file='phase.pad')
    path.set_absorber(x=0.01, y=0.1, z=0.01)
    path.add_scatterer(x=1.806, y=0.1, z=1.806, ipot=1)
    path.degen = 12
    path.calculate_xafs()
    print 'calculate xafs ', path.phase_file, path.nleg
    print 'Atom  IPOT   X, Y, Z'
    for i in range(path.nleg):
        print i, path.ipot[i], path.rat[:,i]

    print path.index, path.degen, path.xmu, path.kf, path.verbose
    print path.ipot
    print path.rat
    npts = 1 + max(np.where(path.kfeff > 0)[0])
    print len(path.kfeff), npts
    print path.kfeff[:5], path.mag_feff[:5], path.rep[:5]

    print "-----------------------"

    # pylab.plot(path.kfeff[:path.nepts], path.mag_feff[:path.nepts])
    # pylab.show()
