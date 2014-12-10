#!/usr/bin/env python

"""
setup.py file for feffpathwrapper (via SWIG)
"""

from distutils.core import setup, Extension
import sys

feffpath_module = Extension('_feffpathwrapper',
                            sources=['feffpath_wrap.c'],
                            libraries=['feffpath'],
                           )

setup (name = 'scatteringpath',
       version = '0.1',
       author      = "Bruce Ravel",
       author_email = 'bravel@bnl.gov',
       url          = 'https://github.com/xraypy/feff85exafs',
       download_url = 'https://github.com/xraypy/feff85exafs',
       description = """A Larch wrapper around the feffpath library""",
       ext_modules = [feffpath_module],
       data_files = [('/usr/local/share/larch/plugins/xafs/', ["feffpathwrapper.py","scatteringpath.py","_feffpathwrapper.so"])],
       )
