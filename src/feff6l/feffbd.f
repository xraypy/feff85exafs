      block data feffbd

      implicit double precision (a-h, o-z)

      character*10 shole(0:9)
      character*8 sout(0:6)
      common /labels/ shole, sout

      include 'vers.h'
c     character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
c     common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

      data shole /'no hole',    'K shell',
     1            'LI shell',   'LII shell',
     2            'LIII shell', 'MI shell',
     3            'MII shell',  'MIII shell',
     4            'MIV shell',  'MV shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 
     1           'DH - HL ', 'DH + HL ', 'HLnoimag', 'Gs HL   '/

c                   123456789012
      data vfeff  /'  Feff 6.02l'/
      data vpotph /'  potph 4.12'/
      data vpaths /'  paths 3.05'/
      data vgenfm /' genfmt 1.44'/
      data vff2ch /' ff2chi 2.01'/

c     6.01l EXAFS only lite version 10/02 jjr
c     5.05a is current working version
c     5.05j is jjr's version 6/93
c     6.00 Alexey's polarization and XANES 
c     6.01 Release version of FEFF6 including bug fixes ala and jjr
c     4.04 Major code reorganization.  Muffin tin finder modified -- now
c     uses average of all possible muffin tin radii instead of minimum.
c     26 March, 1991   Steven Zabinsky
c     4.05 Yet another improvement to muffin tin finder, now averages
c     based on volume of lense-shaped overlapping region, April, 1991
c     4.06 Bug fix in sumax, april 1991
c     4.07 Several minor changes involving non-standard F77 6/6/91, siz
c     4.08 ION card added 7/24/91, siz
c     4.08a, bug in header for ION card fixed 9/10/91, siz
c     4.09, quinn correction added to imhl, interstitial calculation
c           corrected, rmt modified to handle too few neighbors and
c           error msg in phase about hard test in fovrg modified,
c           folp card added
c     POTPH 4.1  Same as feff4.09, but version hacked to work with
c     module potph of feff5, Mar 1992, siz
c
c     new version common added, siz, Mar 1992
c     feff 5.03, first 'real' release, lots of little changes.
c                4 criteria added is the big change.  siz, April 1992
c     feffx 5.04, intermediate intermittent version of code with
c                 background, xsect, xmu, timereversal, lots
c                 of input cards, xanes, etc.  July 1992, siz
c     e REQUIRE card removed, Oct 92, siz
c     f, and paths 3.04, new crits, 9 points. Oct 92
c     g: major bug in xsect -  ixc not passed to xcpot, beginning with
c        5.04g, it's fixed.
c     h use gs for xsect (hard coded)
c     i fixed init and final state mixup in xsect
c     Feff 5.05, release version with all of the above in it.  XANES
c        is turned off in RDINP for the release -- turn it back on
c        there for development.
c     Feff 6 includes polarization (Alexey) and XANES (Steve Z.)
c     Feff 6.01 is the first release version of FEFF6.
c     Feff 6.01l EXAFS only lite version 10/02 jjr

      end
