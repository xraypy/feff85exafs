#!/usr/bin/env python

"""
setup.py file for feffpathwrapper (via SWIG)
"""

from distutils.core import setup, Extension
import sys, os
from platform import uname, architecture
from os.path import join

import larch
if os.name == 'nt':
    installdir = larch.site_configdata.win_installdir
else:
    installdir = larch.site_configdata.unix_installdir


## swiped from larch's dylibs/configure.py
system = uname()[0]
arch   = architecture()[0]
dlldir = None
if os.name == 'nt':
    dlldir = 'win32'
    if arch.startswith('64'):
        dlldir = 'win64'
else:
    if system.lower().startswith('linu'):
        dlldir = 'linux32'
        if arch.startswith('64'):    dlldir = 'linux64'
    elif system.lower().startswith('darw'):
        dlldir = 'darwin'

#print installdir
#print dlldir


feffpath_module = Extension('_feffpathwrapper',
                            sources=['feffpath_wrap.c'],
                            include_dirs = ['C:\\Program Files\\mingw-w64\\x86_64-4.9.2-win32-seh-rt_v3-rev1\\mingw64\\lib\\gcc\\x86_64-w64-mingw32\\4.9.2\\include',],
                            library_dirs = [join('..','..','src','GENFMT')],
                            libraries=['feffpath'],
                           )
# feffpath_module = Extension('_feffpathwrapper',
#                             sources=['feffpath_wrap.c'],
#                             library_dirs = [join('..','..','src','GENFMT')],
#                             libraries=['feffpath'],
#                            )

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
