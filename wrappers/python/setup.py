#!/usr/bin/env python

"""
setup.py file for feffpathwrapper (via SWIG)
"""

from distutils.core import setup, Extension
import sys, os
from os.path import join

import larch
#installdir = larch.site_config.larchdir
installdir = larch.larchlib.larchdir
dlldir =  larch.larchlib.get_dlldir()

if os.name == 'nt':
    feffpath_module = Extension('_feffpathwrapper',
                                sources=['feffpath_wrap.c'],
                                include_dirs = ['C:\\Program Files\\mingw-w64\\x86_64-4.9.2-win32-seh-rt_v3-rev1\\mingw64\\lib\\gcc\\x86_64-w64-mingw32\\4.9.2\\include',],
                                library_dirs = [join('..','..','src','GENFMT')],
                                libraries=['feffpath'],
                            )
else:
    feffpath_module = Extension('_feffpathwrapper',
                                sources=['feffpath_wrap.c'],
                                library_dirs = [join('..','..','src','GENFMT')],
                                libraries=['feffpath'],
                            )

setup (name         = 'scatteringpath',
       version      = '0.1',
       author       = "Bruce Ravel",
       author_email = 'bravel@bnl.gov',
       url          = 'https://github.com/xraypy/feff85exafs',
       download_url = 'https://github.com/xraypy/feff85exafs',
       description  = """A Larch wrapper around the feffpath library""",
       ext_modules  = [feffpath_module],
       data_files   = [(join(installdir,'plugins','xafs'), ["scatteringpath.py",]),
                       (join(installdir,'modules'),        ["feffpathwrapper.py",]),
                       (join(installdir,'dlls',dlldir),    ["_feffpathwrapper.so",])],
       )
