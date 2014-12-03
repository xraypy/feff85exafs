#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


feffpath_module = Extension('_feffpath',
                            sources=['feffpath_wrap.c'],
                            libraries=['feffpath'],
                           )

setup (name = 'feffpath',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Feffpath""",
       ext_modules = [feffpath_module],
       py_modules = ["feffpath"],
       )
