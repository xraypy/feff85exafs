Python interface using ctypes: A work in progress 

Basically, this makes a direct connection to the
Fortran routine "onepath()" to calculate the XAFS
contributions for a scattering path, given path
geometry and a PAD-format file of the Potentials
and Phase-Shifts calculated by Feff.

Status:
==========

   Basic Functionality: Connecting to and running
   "onepath" from the libfeffgenfmt or libonepath
   shared library works.

   A simple calculation for 1 path appears to be
   approximately correct, but needs testing.

    
To Do:
=========

  1. Add a Larch Group interface.  This should
     better match Bruce's scatteringpath module,

  2. Get the unit tests to run!
  
  3. More testing, including Multiple Scattering

  4. Add "calculate_chi()" method to convert results of
     "calculate_xafs()" to chi(k).

  5. Mix with Larch's FeffPath, so that a Path
     can be specified as either FeffDatFile or as
     PotentialsFile + Geometry, and have
     calculation of chi(k) being automatic (or
     whenever geometry changes).

  6. Change the Fortran onepath() to be broken into
     genfmt_prep() and genfmt_calc(), so that the
     Potentials File can be loaded to memory only once.

