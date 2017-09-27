#!/usr/bin/python
"""
Feff85L EXAFS Scatternig Path for Python
"""
from __future__ import print_function
import os
import sys
import ctypes
from ctypes import POINTER, pointer, c_int, c_long, c_char_p, c_double
import numpy as np

def load_feff8lpath():
    dllform = 'lib{:s}.so'
    pathsep = ':'
    loadlib = ctypes.cdll
    if sys.platform.lower().startswith('darwin'):
        dllform = 'lib{:s}.dylib'

    if os.name == 'nt':
        dllform = '{:s}.dll'
        pathsep = ';'
        loadlib = ctypes.windll

    dllname = dllform.format('feff8lpath')
    spath = ['.']
    spath.extend(os.environ.get('LD_LIBRARY_PATH','').split(pathsep))
    spath.extend(os.environ.get('PATH', '').split(pathsep))
    spath.extend(['../../local_install/lib/', '../../src/GENFMT/lib'])
    for dname in spath:
        fullname = os.path.join(dname, dllname)
        if os.path.exists(fullname):
            return loadlib.LoadLibrary(fullname)
    return None


FLIB = load_feff8lpath()


FEFF_maxpts = 150  # nex
FEFF_maxpot = 11   # nphx
FEFF_maxleg = 9    # legtot
BOHR = 0.5291772490


# bytes/str conversion
str2bytes = bytes2str = str
if sys.version_info[0] == 3:
    def bytes2str(val):
        if isinstance(val, str):
            return val
        if isinstance(val, bytes):
            return str(val, 'utf-8')
        return str(val)

    def str2bytes(val):
        if isinstance(val, bytes):
            return val
        return bytes(val, 'utf-8')

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

        self.index   = 1
        self.degen   = 1.
        self.nnnn_out = False
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
        args.exch_label     = str2bytes(' '*8)
        args.genfmt_version = str2bytes(' '*30)

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

        x = FLIB.calc_onepath(args.phase_file, args.index, args.nleg,
                              args.degen, args.genfmt_order,
                              args.exch_label, args.rs, args.vint,
                              args.xmu, args.edge, args.kf, args.rnorman,
                              args.gamach, args.genfmt_version, args.ipot,
                              args.rat, args.iz, args.ipol, args.evec,
                              args.ellip, args.xivec, args.nnnn_out,
                              args.json_out, args.verbose, args.ri,
                              args.beta, args.eta, args.nepts, args.kfeff,
                              args.real_phc, args.mag_feff, args.pha_feff,
                              args.red_fact, args.lam, args.rep)

        self.exch_label = bytes2str(args.exch_label).strip()
        self.genfmt_version = bytes2str(args.genfmt_version).strip()

        for attr in ('index', 'nleg', 'genfmt_order', 'degen', 'rs',
                     'vint', 'xmu', 'edge', 'kf', 'rnorman', 'gamach',
                     'ipol', 'ellip', 'nnnn_out', 'json_out', 'verbose',
                     'nepts'):
            setattr(self, attr, getattr(args, attr).contents.value)

        for attr in ('ipot', 'evec', 'xivec', 'beta', 'eta', 'ri', 'rat',
                    'iz', 'kfeff', 'real_phc', 'mag_feff', 'pha_feff',
                    'red_fact', 'lam', 'rep'):
            setattr(self, attr, np.array(getattr(args, attr).contents[:]))

        # some data needs recasting, reformatting
        self.nnnn_out = bool(self.nnnn_out)
        self.json_out = bool(self.json_out)
        self.verbose  = bool(self.verbose)
        self.rat = self.rat.reshape((2+FEFF_maxleg, 3)).transpose()*BOHR


if __name__ == '__main__':
    path = ScatteringPath(phase_file='phase.pad')
    path.set_absorber( x=0.01,   y=0.1,   z=0.01)
    path.add_scatterer(x=1.8058, y=0.005, z=1.8063, ipot=1)
    path.degen = 12
    path.calculate_xafs()
    print('# Calculate EXAFS with PhaseFile: {:s}'.format(path.phase_file))
    print('# Path Geometry: \n#  IPOT  IZ     X        Y        Z')
    for i in range(path.nleg):
        ipot = path.ipot[i]
        iz   = path.iz[ipot]
        rat  = path.rat[:,i]
        print("#   %2i   %2i  %8.4f %8.4f %8.4f" % (ipot,iz, rat[0], rat[1], rat[2]))

    print("# Polarization: {:d}, ellipticity={:4f}".format(path.ipol, path.ellip))
    print("# Polarization E Vector = {:s}".format(", ".join(["%.4f" % a for a in path.evec])))
    print("# Polarization X Vector = {:s}".format(", ".join(["%.4f" % a for a in path.xivec])))
    print("# Path Settings")
    for attr in ('rs', 'vint', 'xmu', 'edge', 'kf', 'rnorman', 'gamach'):
          print("#   {:8s} = {:+4f} ".format(attr, getattr(path, attr)))
    for attr in ('exch_label', 'genfmt_version'):
          print("#   {:8s} = {:s} ".format(attr, getattr(path, attr)))

    print("Path settings:  degen=%10.5f,  xmu=%10.5f, kf=%10.5f" % (path.degen, path.xmu, path.kf))
    npts = 1 + max(np.where(path.kfeff > 0)[0])
    print("# k         rep          real_phc     phase_feff   mag_feff     red_factor   lambda ")
    fmt = " %6.3f  %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  %11.7f"
    for i in range(int(npts/3.0)):
        print(fmt % (path.kfeff[i], path.rep[i], path.real_phc[i], path.pha_feff[i],
              path.mag_feff[i], path.red_fact[i], path.lam[i]))
        # print(fmt.format(path.kfeff[i], path.rep[i], path.real_phc[i],
        #                  path.mag_feff[i], path.red_fact[i], path.lam[i]))
