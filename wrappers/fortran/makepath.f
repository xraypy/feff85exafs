      program makepath
      implicit double precision (a-h, o-z)

c     taken from feff's HEADERS/dim.h
      integer nex
      parameter (nex = 150)
      integer npatx
      parameter (npatx = 8)
      integer legtot
      parameter (legtot=npatx+1)
      integer nphx
      parameter (nphx = 11)
c      parameter (pi = 3.14159 26535 89793 23846 26433d0)


c+---------------------------------------------------------------------
c     block of parameter declarations for onepath
      character*256 phpad
      character*8 cxc
      character*30 versn
      integer index, iorder, innnn, ixdi, ivrbse, ixc
      double precision evec(3), xivec(3)
      double precision elpty, rs, vint, xmu, edge, xkf, rnrmav, gamach

      double precision rat(3,0:legtot+1)
      integer ipot(0:legtot), iz(0:nphx)

      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
c+---------------------------------------------------------------------

c     initialize everything
      call inipath(index, nleg, deg, iorder,
     &       cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &       versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &       innnn, ixdi, ivrbse, ri, beta, eta,
     &       ne,col1,col2,col3,col4,col5,col6,col7)
      innnn  = 1
      ixdi   = 1
      ivrbse = 1
c      phpad  = 'phase_orig.pad'
      phpad  = 'phase.pad'

c+---------------------------------------------------------------------
c  compute a single path, generating the F matrix then returning the 
c  information contained in a feffNNNN.dat file
c
c  INPUT:
c    index:    path index                            integer
c    nleg:     number of legs in path                integer
c    deg:      path degeneracy                       double
c    iorder:   order of approximation in genfmt      integer
c    ipot:     array of unique potentials            integer(legtot)
c    rat:      cartesian coordinates of scatterers   double(3,0:legtot+1)
c    ipol:     flag to do polarization               integer
c    evec:     polarization vector                   double(3)
c    elpty:    ellipticity                           double
c    xivec:    direction of travel                   double(3)
c    innnn:    flag to write feffNNNN.dat file       integer
c    ixdi:     flag to write feffNNNN.xdi file       integer
c    ivrbse:   flag to write screen messages         integer
c
c    also requires a phase.pad file from an earlier run of xsph
c
c  OUTPUT
c
c    ri:       leg lengths                           double(legtot)
c    beta:     beta angles                           double(legtot+1)
c    eta:      eta angles                            double(legtot+2)
c    ne:       number of k-grid points               integer
c    col1:     k-grid                                double(nex)
c    col2:     central atom phase shifts             double(nex)
c    col3:     magnitude of F_eff                    double(nex)
c    col4:     phase of F_eff                        double(nex)
c    col5:     reduction factor                      double(nex)
c    col6:     mean free path                        double(nex)
c    col7:     real partof complex momentum          double(nex)
c+---------------------------------------------------------------------

c     compute first shell of Copper (SS, deg=12)
      index = 1
      nleg  = 2
      deg   = 12
      call addatom(1,  1.805, 0.,  1.805, 1, ipot, rat)

      call onepath(phpad, index, nleg, deg, iorder,
     &     cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &     versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &     innnn, ixdi, ivrbse, ri, beta, eta,
     &     ne,col1,col2,col3,col4,col5,col6,col7)

c     this bit writes the data table from feff0001.dat to the screen
c         do 100 i=1,ne
c            write(*,400) col1(i), col2(i), col3(i), col4(i), col5(i),
c        &          col6(i), col7(i)
c    100  continue
c    400  format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
c        1       2(1pe11.4,1x))


c     compute fourth shell of Copper (DS, deg=48)
      call inipath(index, nleg, deg, iorder,
     &       cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &       versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &       innnn, ixdi, ivrbse, ri, beta, eta,
     &       ne,col1,col2,col3,col4,col5,col6,col7)
      innnn  = 1
      ixdi   = 1
      ivrbse = 1

      index = 4
      nleg  = 3
      deg   = 48
      call addatom(1,  0.,     0., -3.61,  1, ipot, rat)
      call addatom(2, -1.805,  0., -1.805, 1, ipot, rat)

      call onepath(phpad, index, nleg, deg, iorder,
     &       cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &       versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &       innnn, ixdi, ivrbse, ri, beta, eta,
     &       ne,col1,col2,col3,col4,col5,col6,col7)

      end


      subroutine addatom(leg, x, y, z, ip, ipot, rat)
      implicit double precision (a-h, o-z)
      integer npatx, legtot
      parameter (npatx = 8, legtot=npatx+1)

      integer leg, ip
      real x, y, z
      double precision rat(3,0:legtot+1)
      integer ipot(0:legtot)
      double precision bohr
      parameter (bohr = 0.529 177 249d0)

      ipot(leg)  = ip
      rat(1,leg) = dble(x) / bohr
      rat(2,leg) = dble(y) / bohr
      rat(3,leg) = dble(z) / bohr
      return
      end


      subroutine inipath(index, nleg, deg, iorder,
     &       cxc, rs, vint, xmu, edge, xkf, rnrmav, gamach,
     &       versn, ipot, rat, iz, ipol, evec, elpty, xivec,
     &       innnn, ixdi, ivrbse, ri, beta, eta,
     &       ne,col1,col2,col3,col4,col5,col6,col7)
      implicit double precision (a-h, o-z)

c     taken from feff's HEADERS/dim.h
      integer nex, npatx, legtot, ixc
      parameter (nex = 150, npatx = 8, legtot=npatx+1, nphx=11)

      integer index, iorder, innnn, ivrbse
      double precision evec(3), xivec(3)
      double precision elpty,rs, vint, xmu, edge, xkf, rnrmav, gamach

      double precision rat(3,0:legtot+1)
      integer ipot(0:legtot), iz(0:nphx)

      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      character*256 phpad
      character*8 cxc
      character*30 versn

      index  = 9999
      nleg   = 0
      deg    = 1.
      iorder = 2
      innnn  = 0
      ixdi   = 0 
      ivrbse = 0
      ipol   = 0 
      elpty  = 0.
      ne     = 0

      versn  = ""
      cxc    = ""
      rs     = 0.
      vint   = 0.
      xmu    = 0.
      edge   = 0.
      xkf    = 0.
      rnrmav = 0.
      gamach = 0.
      
      do 5  i=1,3
         evec(i)  = 0
         xivec(i) = 0
 5    continue
      do 10 i=0, legtot
         rat(1,i) = 0
         rat(2,i) = 0
         rat(3,i) = 0
         ipot(i)  = 0
         eta(i)   = 0
         if (i>0) then
            ri(i)   = 0
            beta(i) = 0
            eta(i)  = 0
         endif
 10   continue
      beta(legtot+1) = 0
      eta(legtot+1)  = 0

      do 15 i=0,nphx
         iz(i) = 0
 15   continue

      do 20 i=1, nex
         col1(i) = 0
         col2(i) = 0
         col3(i) = 0
         col4(i) = 0
         col5(i) = 0
         col6(i) = 0
         col7(i) = 0
 20   continue

      return
      end
