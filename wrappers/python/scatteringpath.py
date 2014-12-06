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
     a.deg = 48
  except ValueError:
     print "Bad degeneracy value, deg still %" : a.deg
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

  phbin        : string   path to `phase.bin`                       `phase.bin`       
  index        : integer  path index                                9999              
  deg          : float    path degeneracy                           required input    
  nleg         : integer  number of legs in path                    set using atom method (read-only)
  rat          : array    cartesian positions of atoms in path      set using atom method (not yet available)
  ipot         : integer  unique potentials of atoms in path        set using atom method (not yet available)
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

  errorcode    : integer  error code from `atom` or `make`
  errormessage : string   error message from `atom` or `make`

"""

# Missing features:
#
#  * Getter method(s) for rat and ipot (this is a missing feature in
#    the feffpath library


# LICENSE AND COPYRIGHT
#
# To the extent possible, the authors have waived all rights granted by
# copyright law and related laws for the code and documentation that
# make up the Python Interface to the feffpath library.  While information
# about Authorship may be retained in some files for historical reasons,
# this work is hereby placed in the Public Domain.  This work is
# published from: United States.
#
# Note that the onepath library itself is NOT public domain, nor is the
# Fortran source code for Feff that it relies upon.
#
# Author: Bruce Ravel (bravel AT bnl DOT gov).
# Last update: 4 December, 2014

from larch import (Group, Parameter, isParameter, ValidateLarchPlugin, param_value,
                   use_plugin_path, isNamedClass, Interpreter)
use_plugin_path('xafs')
import feffpathwrapper
from   os.path   import isfile

class FeffPath(Group):
    """
    A larchified, simplified, feature-full wrapper around the thin
    (but weird) SWIG wrapper for the feffpath library.
    """

    def __init__(self, folder=None, _larch=None, **kws):
        kwargs = dict(name='FeffPath wrapper')
        kwargs.update(kws)
        Group.__init__(self,  **kwargs)
        self._larch     = Interpreter()
        self.wrapper    = feffpathwrapper.FEFFPATH()
        feffpathwrapper.create_path(self.wrapper)
        self.wrapper.phbin = ''


    ## ---- scalar valued attributs

    @property
    def phbin(self):
        return self.wrapper.phbin
    @phbin.setter
    def phbin(self,value):
        if not isfile(value):
            raise ValueError("%s is not a readable file" % value)
        self.wrapper.phbin = value

    @property
    def index(self):
        return int(self.wrapper.index)
    @index.setter
    def index(self,value):
        value = int(value)
        if value < 1:
            raise ValueError("Non-positive value of index not allowed: %s" % value)
        elif value > 9999:
            raise ValueError("Index value must be four digits or less: %s" % value)
        self.wrapper.index = long(value)

    @property
    def nleg(self):
        return int(self.wrapper.nleg)
    @nleg.setter
    def nleg(self,value):
        pass

    @property
    def deg(self):
        return self.wrapper.deg
    @deg.setter
    def deg(self,value):
        if value < 1:
            raise ValueError("Negative degeneracies not allowed: %s" % value)
        self.wrapper.deg = value

    @property
    def iorder(self):
        return int(self.wrapper.iorder)
    @iorder.setter
    def iorder(self,value):
        if value < 0:
            raise ValueError("Negative value of iorder not allowed: %s" % value)
        elif value > 10:
            raise ValueError("Iorder value must be 10 or less: %s" % value)
        self.wrapper.iorder = long(value)

    @property
    def nnnn(self):
        return self.wrapper.nnnn
    @nnnn.setter
    def nnnn(self,value):
        self.wrapper.nnnn = self.bool(value)

    @property
    def json(self):
        return self.wrapper.json
    @json.setter
    def json(self,value):
        self.wrapper.json = self.bool(value)

    @property
    def verbose(self):
        return self.wrapper.verbose
    @verbose.setter
    def verbose(self,value):
        self.wrapper.verbose = self.bool(value)

    @property
    def ipol(self):
        return self.wrapper.ipol
    @ipol.setter
    def ipol(self,value):
        self.wrapper.ipol = self.bool(value)

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

    @property
    def ne(self):
        return int(self.wrapper.ne)
    @ne.setter
    def ne(self,value):
        pass


    ## ---- error handling

    @property
    def errorcode(self):
        return int(self.wrapper.errorcode)
    @errorcode.setter
    def ne(self,value):
        pass

    @property
    def errormessage(self):
        return self.wrapper.errormessage
    @errormessage.setter
    def ne(self,value):
        pass


    ## ---- array valued attributes

    @property
    def evec(self):
        return [feffpathwrapper.get_evec(self.wrapper,0), feffpathwrapper.get_evec(self.wrapper,1), feffpathwrapper.get_evec(self.wrapper,2)]
    @evec.setter
    def evec(self, vec):
        if type(value) not in ('list', 'tuple'):
            raise ValueError("evec must be 3 element list (or tuple) with (x,y,z) of the electric vector")
        elif len(value) != 3:
            raise ValueError("evec must be 3 element list (or tuple) with (x,y,z) of the electric vector")
        feffpathwrapper.set_evec(self.wrapper,0,vec[0])
        feffpathwrapper.set_evec(self.wrapper,1,vec[1])
        feffpathwrapper.set_evec(self.wrapper,2,vec[2])
        self.wrapper.ipol = True

    @property
    def xivec(self):
        return [feffpathwrapper.get_xivec(self.wrapper,0), feffpathwrapper.get_xivec(self.wrapper,1), feffpathwrapper.get_xivec(self.wrapper,2)]
    @xivec.setter
    def xivec(self, vec):
        if type(value) not in ('list', 'tuple'):
            raise ValueError("xivec must be 3 element list (or tuple) with (x,y,z) of the Poynting vector")
        elif len(value) != 3:
            raise ValueError("xivec must be 3 element list (or tuple) with (x,y,z) of the Poynting vector")
        feffpathwrapper.set_xivec(self.wrapper,0,vec[0])
        feffpathwrapper.set_xivec(self.wrapper,1,vec[1])
        feffpathwrapper.set_xivec(self.wrapper,2,vec[2])
        self.wrapper.ipol = True

    @property
    def ri(self):
        arr=[]
        for i in range(self.wrapper.nleg):
            arr.append(feffpathwrapper.get_ri(self.wrapper,i))
        return arr
    @ri.setter
    def ri(self, x, y, z):  pass

    @property
    def beta(self):
        arr=[]
        for i in range(self.wrapper.nleg):
            arr.append(feffpathwrapper.get_beta(self.wrapper,i))
        return arr
    @beta.setter
    def beta(self, x, y, z):  pass

    @property
    def eta(self):
        arr=[]
        for i in range(self.wrapper.nleg):
            arr.append(feffpathwrapper.get_eta(self.wrapper,i))
        return arr
    @eta.setter
    def eta(self, x, y, z):  pass


    ## ---- columns of feffNNNN.dat

    @property
    def k(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_k(self.wrapper,i))
        return arr
    @k.setter
    def k(self, value): pass

    @property
    def real_phc(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_real_phc(self.wrapper,i))
        return arr
    @real_phc.setter
    def real_phc(self, value): pass

    @property
    def mag_feff(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_mag_feff(self.wrapper,i))
        return arr
    @mag_feff.setter
    def mag_feff(self, value): pass

    @property
    def pha_feff(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_pha_feff(self.wrapper,i))
        return arr
    @pha_feff.setter
    def pha_feff(self, value): pass

    @property
    def red_fact(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_red_fact(self.wrapper,i))
        return arr
    @red_fact.setter
    def red_fact(self, value): pass

    @property
    def lam(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_lam(self.wrapper,i))
        return arr
    @lam.setter
    def lam(self, value): pass

    @property
    def realp(self):
        arr=[]
        for i in range(self.wrapper.ne):
            arr.append(feffpathwrapper.get_realp(self.wrapper,i))
        return arr
    @realp.setter
    def realp(self, value): pass



    ## ---- methods

    def atom(self, x, y, z, ip):
        """
        Add a scatterer to a scattering path and increment nleg.

          a = scatteringpath()
          a.atom(0, 0, -3.61, 1)

        Arguments are x,y,z coordinates of the scatterer and it ipot value

        The Cartesian coordinates are relative to the ABSORBER AT 0,0,0
        """
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

    def clear(self):
        """
        Reinitialize a path.

          a.clear()

        """
        feffpathwrapper.clear_path(self.wrapper)



######################################################################

def scatpath(_larch=None, **kws):
    """
    Make and return an empty FeffPath group
    """
    return FeffPath(_larch=_larch)

def registerLarchPlugin(): # must have a function with this name!
    return ('_xafs', { 'scatteringpath': scatpath })
