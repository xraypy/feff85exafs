c     Periodic table of the elements
c     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

c     atwts(iz)  single precision fn, returns atomic weight
c     atwtd(iz)  double precision fn, returns atomic weight
c     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd (iz)
      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real function atwts (iz)
      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      character*2 sym
      save /atsyco/
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
c     PeRiodic TaBle Block Data

c     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)

      character*2 sym
      save /atsyco/
      common /atsyco/ sym(103)

      data weight /
     1   1.0079, 4.0026, 6.941,  9.0122, 10.81,   12.01,
     2   14.007, 15.999, 18.998, 20.18,  22.9898, 24.305,
     3   26.982, 28.086, 30.974, 32.064, 35.453,  39.948,
     4   39.09,  40.08,  44.956, 47.90,  50.942,  52.00,
     5   54.938, 55.85,  58.93,  58.71,  63.55,   65.38,
     6   69.72,  72.59,  74.922, 78.96,  79.91,   83.80,
     7   85.47,  87.62,  88.91,  91.22,  92.91,   95.94,
     8   98.91,  101.07, 102.90, 106.40, 107.87,  112.40,
     9   114.82, 118.69, 121.75, 127.60, 126.90,  131.30,
     x   132.91, 137.34, 138.91, 140.12, 140.91,  144.24,
     1   145,    150.35, 151.96, 157.25, 158.92,  162.50,
     2   164.93, 167.26, 168.93, 173.04, 174.97,  178.49,
     3   180.95, 183.85, 186.2,  190.20, 192.22,  195.09,
     4   196.97, 200.59, 204.37, 207.19, 208.98,  210,
     5   210,    222,    223,    226,    227,     232.04,
     6   231,    238.03, 237.05, 244,    243,     247,
     7   247,    251,    254,    257,    256,     254,
     8   257/

      data sym /  'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F', 'Ne',
     1            'Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K', 'Ca',
     2            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr',
     4            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5            'Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr','Nd',
     6            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7            'Lu','Hf','Ta','W', 'Te','Os','Ir','Pt','Au','Hg',
     8            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9            'Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     x            'Md','No','Lw'/

      end
