#!/usr/bin/env python
"""A larchified, simplified, feature-full wrapper around the thin
(but weird) SWIG wrapper for the feffpath library.

Use of this plugin presumes that the feffpath C library is installed
and that the SWIG wrapper has been built, compiled, and installed.

  add_plugin("scatteringpath")
  a = scatteringpath()
  try:
     a.index = 1
  except ValueError:
     print "Bad index value, index still %" : a.index
  #endtry
  try:
     a.degen = 48
  except ValueError:
     print "Bad degeneracy value, degen still %" : a.degen
  #endtry
  a.atom(1.805,0,1.805,1)
  if a.errorcode:
     print a.errormessage
     exit()
  #endif
  a.make()
  if a.errorcode:
     print a.errormessage
     exit()
  #endif
  a.clear()

Methods:

  atom:   add a scattering atom to a path and increment nleg
  make:   compute F_eff for a scattering path
  clear:  reinitialize  the scattering path

Attributes (input):

  phpad        : string   path to `phase.pad`                       `phase.pad`       
  index        : integer  path index                                9999              
  degen        : float    path degeneracy                           required input    
  nleg         : integer  number of legs in path                    set using atom method (read-only)
  rat          : list of tuples  cartesian positions of atoms in path      set using atom method
  ipot         : list of int     unique potentials of atoms in path        set using atom method
  iorder       : integer  order of approximation in genfmt          2
         
  nnnn         : boolean  flag to write `feffNNNN.dat` file         False             
  json         : boolean  flag to write `feffNNNN.json` file        False             
  verbose      : boolean  flag to write screen messages             False             

  ipol         : boolean  flag to do polarization calculation       False (set True when evec|xivec|elpty set)
  evec         : array    polarization vector                       (0,0,0)           
  elpty        : float    ellipticity                               0                 
  xivec        : array    direction of X-ray propagation            (0,0,0)           

Attributes (output):

  ri           : array    leg lengths                                                
  beta         : array    beta angles                                                
  eta          : array    eta angles                                                 
  reff         : float    half path length (computed from ri)

  ne           : integer  number of energy points actually used by Feff (read-only)
  k            : array    k grid for feff path calculation, column 1 in `feffNNNN.dat`
  real_phc     : array    central atom phase shifts. column 2 in `feffNNNN.dat`
  mag_feff     : array    magnitude of F_eff, column 3 in `feffNNNN.dat`
  pha_feff     : array    phase of F_eff, column 4 in `feffNNNN.dat`
  red_fact     : array    reduction factor, column 5 in `feffNNNN.dat`
  lam          : array    mean free path, column 6 in `feffNNNN.dat`
  rep          : array    real part of complex momentum, column 7 in `feffNNNN.dat`

Output attributes related to Larch's feffdat Group
  geom         : list     path geometry: list of (Symbol, Z, ipot, x, y, z)
  amp          : array    mag_feff * red_fact
  pha          : array    pha_feff + real_phc

Error handling attributes
  errorcode    : integer  error code from `atom` or `make`
  errormessage : string   error message from `atom` or `make`

The following attributes are related to Feff's potential model
  edge         : float    energy threshold relative to atomic value (a poor estimate)
  gam_ch       : float    core level energy width
  kf           : float    k value at Fermi level
  mu           : float    Fermi level, eV
  rnorman      : float    Norman radius
  version      : string   Feff version
  exch         : string   brief description of the electronic exchange model
  rs_int       : float    interstitial radius
  vint         : float    interstitial potential
  version      : string   the version of feff and the feffpath revision

This attribute of Larch's _feffdat group is not captured
  potentials   : list     path potentials: list of (ipot, z, r_MuffinTin, r_Norman)

"""

# LICENSE AND COPYRIGHT
#
# To the extent possible, the authors have waived all rights granted by
# copyright law and related laws for the code and documentation that
# make up the Python Interface to the feffpath library.  While information
# about Authorship may be retained in some files for historical reasons,
# this work is hereby placed in the Public Domain.  This work is
# published from: United States.
#
# Note that the onepath and feffpath libraries themselves are NOT
# public domain, nor is the Fortran source code for Feff that it
# relies upon.
#
# Author: Bruce Ravel (bravel AT bnl DOT gov).
# Last update: 12 December, 2014

import larch
from larch import (Group, Parameter, isParameter, ValidateLarchPlugin, param_value,
                   use_plugin_path, isNamedClass, Interpreter)
use_plugin_path('xafs')
use_plugin_path('xray')
from   xraydb_plugin import atomic_symbol

import sys
from   os            import name
from   os.path       import isfile
from   numpy         import array

#installdir = larch.site_config.larchdir
installdir = larch.larchlib.larchdir
dlldir     = larch.larchlib.get_dlldir()


dllfile=installdir+'/dlls/'+dlldir
if not dllfile in sys.path:
    sys.path.append(dllfile)
moddir=installdir+'/modules'
if not moddir in sys.path:
    sys.path.append(moddir)

import feffpathwrapper



class FeffPathBoolean(object):
    """A descriptor for boolean-valued FeffPath attributes"""
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        return getattr(instance.wrapper, self._name)

    def __set__(self, instance, val):
        setattr(instance.wrapper, self._name, val)

class FeffPathUnsettableInteger(object):
    """A descriptor for integer-valued FeffPath attributes which are not
    intended to be set by the user, instead are set by calls to the
    atom or make methods.

    """
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        return int(getattr(instance.wrapper, self._name))

    def __set__(self, instance, val):
        pass

class FeffPathUnsettableString(object):
    """A descriptor for string-valued FeffPath attributes which are not
    intended to be set by the user, instead are set by calls to the
    make method.

    """
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        return str(getattr(instance.wrapper, self._name))

    def __set__(self, instance, val):
        pass

class FeffPathUnsettableFloat(object):
    """A descriptor for float-valued FeffPath attributes which are not
    intended to be set by the user, instead are set by calls to the
    make method.

    """
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        return float(getattr(instance.wrapper, self._name))

    def __set__(self, instance, val):
        pass

class FeffPathNlegList(object):
    """A descriptor for list-valued FeffPath attributes which have nleg
    members and describe the geometry of the path, i.e. ri, beta, and
    eta

    """
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        lst=[]
        for i in range(instance.wrapper.nleg):
            x = eval('feffpathwrapper.get_'+self._name+'(instance.wrapper,i)')
            lst.append(float(x))
            #lst.append(feffpathwrapper.get_ri(self.wrapper,i))
        return lst

    def __set__(self, instance, val):
        pass

class FeffPathColumn(object):
    """A descriptor for the array-valued FeffPath attributes which are the
    columns of feffNNNN.dat"""
    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        arr=[]
        for i in range(instance.wrapper.ne):
            x = eval('feffpathwrapper.get_'+self._name+'(instance.wrapper,i)')
            arr.append(float(x))
            #arr.append(feffpathwrapper.get_k(instance.wrapper,i))
        return array(arr)

    def __set__(self, instance, val):
        pass



class FeffPath(Group):
    """
    A larchified, simplified, feature-full wrapper around the thin
    (but weird) SWIG wrapper for the feffpath library.
    """

    nnnn      = FeffPathBoolean('nnnn')
    json      = FeffPathBoolean('json')
    verbose   = FeffPathBoolean('verbose')
    ipol      = FeffPathBoolean('ipol')

    nleg      = FeffPathUnsettableInteger('nleg')
    ne        = FeffPathUnsettableInteger('ne')
    errorcode = FeffPathUnsettableInteger('errorcode')

    exch      = FeffPathUnsettableString('exch')
    version   = FeffPathUnsettableString('version')
    errormessage = FeffPathUnsettableString('errormessage')

    edge      = FeffPathUnsettableFloat('edge')
    gam_ch    = FeffPathUnsettableFloat('gam_ch')
    kf        = FeffPathUnsettableFloat('kf')
    mu        = FeffPathUnsettableFloat('mu')
    rnorman   = FeffPathUnsettableFloat('rnorman')
    rs_int    = FeffPathUnsettableFloat('rs_int')
    vint      = FeffPathUnsettableFloat('vint')

    ri        = FeffPathNlegList('ri')
    beta      = FeffPathNlegList('beta')
    eta       = FeffPathNlegList('eta')

    k         = FeffPathColumn('k')
    real_phc  = FeffPathColumn('real_phc')
    mag_feff  = FeffPathColumn('mag_feff')
    pha_feff  = FeffPathColumn('pha_feff')
    red_fact  = FeffPathColumn('red_fact')
    lam       = FeffPathColumn('lam')
    rep       = FeffPathColumn('rep')

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='FeffPath wrapper')
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self._larch     = Interpreter()
        self.wrapper    = feffpathwrapper.FEFFPATH()
        feffpathwrapper.create_path(self.wrapper)
        self.wrapper.phpad = ''
        self.ipot = []
        self.rat  = []
        self.geom = []


    ## ---- scalar valued attributs, these get their own properties so
    ##      tailored ValueErrors can be issued

    @property
    def phpad(self):
        return self.wrapper.phpad
    @phpad.setter
    def phpad(self,value):
        if not isfile(value):
            raise ValueError("%s is not a readable file" % value)
        self.wrapper.phpad = value

    @property
    def index(self):
        return int(self.wrapper.index)
    @index.setter
    def index(self,value):
        value = int(value)
        if value < 1 or value > 9999:
            raise ValueError("Index values must be between 1 and 9999: %s" % value)
        self.wrapper.index = long(value)

    @property
    def degen(self):
        return self.wrapper.degen
    @degen.setter
    def degen(self,value):
        if value < 1:
            raise ValueError("Negative degeneracies not allowed: %s" % value)
        self.wrapper.degen = value

    @property
    def iorder(self):
        return int(self.wrapper.iorder)
    @iorder.setter
    def iorder(self,value):
        if value < 0 or value > 10:
            raise ValueError("Iorder value must be 10 or less: %s" % value)
        self.wrapper.iorder = long(value)

    @property
    def elpty(self):
        return self.wrapper.elpty
    @elpty.setter
    def elpty(self,value):
        if value < 0:
            raise ValueError("elpty must be between 0 and 1: %s" % value)
        elif value > 1:
            raise ValueError("elpty must be between 0 and 1: %s" % value)
        self.wrapper.elpty = value
        self.wrapper.ipol  = True

    ## ---- error handling (this is the only string-valued attribute, it is returned from atom and make)

    # @property
    # def errormessage(self):
    #     return self.wrapper.errormessage
    # @errormessage.setter
    # def errormessage(self,value):
    #     pass


    ## ---- 3-vec valued attributes (with tailored ValueErrors)

    @property
    def evec(self):
        return [feffpathwrapper.get_evec(self.wrapper,0),
                feffpathwrapper.get_evec(self.wrapper,1),
                feffpathwrapper.get_evec(self.wrapper,2)]
    @evec.setter
    def evec(self, vec):
        if type(vec).__name__ not in ('list', 'tuple'):
            raise ValueError("evec must be 3 element list (or tuple) with (x,y,z) of the electric vector (not a list/tuple)")
        elif len(vec) != 3:
            raise ValueError("evec must be 3 element list (or tuple) with (x,y,z) of the electric vector (not of length 3)")
        feffpathwrapper.set_evec(self.wrapper,0,vec[0])
        feffpathwrapper.set_evec(self.wrapper,1,vec[1])
        feffpathwrapper.set_evec(self.wrapper,2,vec[2])
        self.wrapper.ipol = True

    @property
    def xivec(self):
        return [feffpathwrapper.get_xivec(self.wrapper,0),
                feffpathwrapper.get_xivec(self.wrapper,1),
                feffpathwrapper.get_xivec(self.wrapper,2)]
    @xivec.setter
    def xivec(self, vec):
        if type(vec).__name__ not in ('list', 'tuple'):
            raise ValueError("xivec must be 3 element list/tuple with (x,y,z) of the Poynting vector (not a list/tuple)")
        elif len(vec) != 3:
            raise ValueError("xivec must be 3 element list/tuple with (x,y,z) of the Poynting vector (not of length 3)")
        feffpathwrapper.set_xivec(self.wrapper,0,vec[0])
        feffpathwrapper.set_xivec(self.wrapper,1,vec[1])
        feffpathwrapper.set_xivec(self.wrapper,2,vec[2])
        self.wrapper.ipol = True


    @property
    def iz(self):
        lst = [];
        for i in range(feffpathwrapper.nphx):
            lst.append(int(feffpathwrapper.get_iz(self.wrapper,i)))
        return lst
    @iz.setter
    def iz(self, vec):
        pass

    @property
    def amp(self):
        return self.mag_feff * self.red_fact
    @amp.setter
    def amp(self, value):
        pass

    @property
    def pha(self):
        return self.pha_feff + self.real_phc
    @pha.setter
    def pha(self, value):
        pass


    ## ---- methods

    def atom(self, x, y, z, ip):
        """
        Add a scatterer to a scattering path and increment nleg.

          a = scatteringpath()
          a.atom(0, 0, -3.61, 1)

        Arguments are x,y,z coordinates of the scatterer and it ipot value

        The Cartesian coordinates are relative to the ABSORBER AT 0,0,0
        """
        self.ipot.append(ip)
        self.rat.append((x, y, z))
        feffpathwrapper.add_scatterer(self.wrapper, x, y, z, ip)

    def make(self):
        """
        Compute F_eff for a path, make columns of feffNNNN.dat available as
        properties of the object.  Write feffNNNN.dat file if verbose attribute
        set to True.

          a = scatteringpath()
          a.atom(0, 0, -3.61, 1)
          a.make()

        """
        feffpathwrapper.make_path(self.wrapper)
        if self.errorcode == 0:
            self.geom = []
            for i in range(self.nleg-1):
                ip = self.ipot[i]
                this = (str(atomic_symbol(self.iz[ip], _larch=self._larch)),
                        self.iz[ip], ip, self.rat[i][0], self.rat[i][1], self.rat[i][2])
                self.geom.append(this)

    def clear(self):
        """
        Reinitialize a path.

          a.clear()

        """
        feffpathwrapper.clear_path(self.wrapper)
        self.ipot = []
        self.rat  = []
        self.geom = []



######################################################################

def scatpath(_larch=None, **kws):
    """
    Make and return an empty FeffPath group
    """
    return FeffPath(_larch=_larch)

def registerLarchPlugin(): # must have a function with this name!
    return ('_xafs', { 'scatteringpath': scatpath })
