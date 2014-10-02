# Wrapper source

This folder is for the SWIG-generated wrapper around the `feffpath`
library.  Once the files in the `src/GENFMT` folder of the feff85exafs
distribution are built and installed, this folder will be populated
with these files:

 * `feffpath.h`: header file
 * `feffpath_wrap.c`: SWIG wrapper

These will be compiled as part of the Perl build procedure into the
share object library used by the Perl wrapper.
