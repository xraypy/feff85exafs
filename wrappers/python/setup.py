#!/usr/bin/env python

"""
setup.py file for feffpathwrapper (via SWIG)
"""

from distutils.core import setup, Extension


feffpath_module = Extension('_feffpathwrapper',
                            sources=['feffpath_wrap.c'],
                            libraries=['feffpath'],
                           )

setup (name = 'feffpathwrapper',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Thin Python wrapper around the feffpath library""",
       ext_modules = [feffpath_module],
       py_modules = ["feffpathwrapper"],
       )
