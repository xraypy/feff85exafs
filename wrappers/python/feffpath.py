#!/usr/bin/env python
"""
Compute Feff scattering paths
"""
import os
import ctypes
import ctypes.util

class FEFFPATHStruct(ctypes.Structure):
    "emulate the FEFFPATH struct"
    _fields_ = [('index',    ctypes.c_long),
                ('nleg',     ctypes.c_long),
                ('deg',      ctypes.c_double),
                ('iorder',   ctypes.c_long),
                ('nnnn',     ctypes.c_bool),
                ('json',     ctypes.c_bool),
                ('verbose',  ctypes.c_bool),
                ('ipol',     ctypes.c_bool),
                ('evec',     ctypes.c_void_p),
                ('elpty',    ctypes.c_double),
                ('xivec',    ctypes.c_void_p),
                ('ri',       ctypes.c_void_p),
                ('beta',     ctypes.c_void_p),
                ('eta',      ctypes.c_void_p),
                ('reff',     ctypes.c_double),
                ('ne',       ctypes.c_long),
                ('k',        ctypes.c_void_p),
                ('real_phc', ctypes.c_void_p),
                ('mag_feff', ctypes.c_void_p),
                ('pha_feff', ctypes.c_void_p),
                ('red_fact', ctypes.c_void_p),
                ('lam',      ctypes.c_void_p),
                ('rep',      ctypes.c_void_p)]


def get_dll(libname):
    """find and load a shared library, swiped from Larch;s lib/larchlib.py"""
    _paths = {'PATH': '', 'LD_LIBRARY_PATH': '', 'DYLD_LIBRARY_PATH':''}
    _dylib_formats = {'win32': '%s.dll', 'linux2': 'lib%s.so',
                      'darwin': 'lib%s.dylib'}
    thisdir = os.path.abspath(os.path.join(sys_larchdir, 'dlls',
                                           get_dlldir()))
    dirs = [thisdir]

    loaddll = ctypes.cdll.LoadLibrary

    if sys.platform == 'win32':
        loaddll = ctypes.windll.LoadLibrary
        dirs.append(sys_larchdir)

    if hasattr(sys, 'frozen'): # frozen with py2exe!!
        dirs.append(os.path.dirname(sys.executable))

    for key in _paths:
        for d in dirs:
            _paths[key] = add2path(key, d)

    # normally, we expect the dll to be here in the larch dlls tree
    # if we find it there, use that one
    fname = _dylib_formats[sys.platform] % libname
    dllpath = os.path.join(thisdir, fname)
    if os.path.exists(dllpath):
        return loaddll(dllpath)

    # if not found in the larch dlls tree, try your best!
    return loaddll(ctypes.util.find_library(libname))


FEFFPATHLIB = None
def get_feffpathlib():
    """make initial connection to FEFFPATH dll"""
    global FEFFPATHLIB
    if FEFFPATHLIB is None:
        FEFFPATHLIB = get_dll('feffpath')
    return FEFFPATHLIB



class FeffPath(object):
    """ Feff path:

    See https://github.com/xraypy/feff85exafs

    for further details

    >>> path = FeffPath()

    Principle methods:
      this(): 
      that(): 
    """
    _invalid_msg = "invalid data for '%s':  was expecting %s, got '%s'"

    def __init__(self, filename=None):
        self.feffpath = get_feffpathlib()
        self.status = None


    def create(self):
        """map the FEFFPATH struct onto a python object
        """
        pfeffpath = ctypes.pointer(FEFFPATHStruct())
        self.status = out = self.xdilib.XDI_readfile(filename, pxdi)
        if out < 0:
            msg =  self.xdilib.XDI_errorstring(out)
            msg = 'Error reading XDIFile %s\n%s' % (filename, msg)
            raise XDIFileException(msg        



