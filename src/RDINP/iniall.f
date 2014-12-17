      subroutine iniall
c     initializes all input variables contained in the
c     common blocks of the header file allinp.h 
c     written by Alexei Ankudinov , march 2001.
c     following the suggested by Bruce Ravel subroutine iorini
      implicit double precision (a-h, o-z)

      real szero, sone
      double precision dzero, done
      parameter(szero = 0.e0, dzero = 0.d0)
      parameter(sone = 1.e0,  done = 1.d0)
      include '../HEADERS/dim.h'
      include '../RDINP/allinp.h'

c     initialize character string arrays
      do 10 i=1,nheadx
        title(i) = ' '
 10   continue

c  initialize integer scalars
      iGrid = 0 ! Josh Kas
      ntitle = 0
      nat = 0
      natt = 0
      nph = 0

      iafolp = 0
      idwopt = -1
      ihole = 1
      inters = 0
      iorder = 2
      ipr1 = 0
      ipr2 = 0
      ipr3 = 0
      ipr4 = 0
      ipr5 = 0
      ipr6 = 0
      ipse = 0
      ipsk = 0
      ispec = 0
      ixc = 0
      ixc0 = -1
      jumprm = 0
      lfms1 = 0
      lfms2 = 0
      minv = 0
      lreal = 0
      mbconv = 0
      mchi = 1
      mfeff = 1
      mfms = 1
      mpath = 1
      mphase = 1
      mldos = 0
      mpot = 1
      ms = 0
      iPlsmn = 0 ! Josh Kas
      mso2conv = 0 ! Josh Kas
      nlegxx = 10
      nmix = 1
      nohole = -1
      nscmt = 0
      icoul = 0
      iunf = 0
      izstd = 0
      ifxc = 0
      ipmbse = 0
      itdlda = 0
      nonlocal = 0
      ibasis = 0

cc initialize reals
      critpw = 2.5*sone
      pcritk = szero
      pcrith = szero
      rmax = -1 * sone
      rfms1 = -1 * sone
      rfms2 = -1 * sone
      rdirec = -1 * sone
      toler1 = 1.e-3
      toler2 = 1.e-3

cc initialize double precision scalars
      alphat = dzero
      thetae = dzero
      ca1 = dzero
      critcw = 4*done
      eimag = -1*done
      ecv = -40*done 
      emax = dzero
      emin = 1000*done
      rclabs = dzero
      rgrd = 0.05 * done
      s02 = done
      sig2g = dzero
      thetad = dzero
      tk = dzero
      totvol = dzero
      vr0 = dzero
      vi0 = dzero
      vicorr = dzero
      vrcorr = dzero
      xkmax = 20*done
      xkstep = 0.07*done
      vixan = dzero
      wsigk = dzero ! Josh Kas
      cen = dzero ! Josh Kas
      
cc initialize logicals
      wnstar = .false.


c  initialize loops of number of potentials
      do 110 i=0,nphx
        xnatph(i) = dzero
        spinph(i) = dzero
        iz(i) = 0
        xion(i) = dzero
        folp(i) = done
        novr(i) = 0
        lmaxsc(i) = 0
        lmaxph(i) = 0
        potlbl(i) = ' '
 110  continue

      do 114 i=0,nphx
         do 112 j=1,novrx
            iphovr(j,i)=0
            nnovr(j,i)=0
            rovr(j,i)=dzero
 112     continue
 114  continue

c  initialize polarization data
      ipol = 0
      ispin = 0
      le2 = 0
      l2lp = 0
      elpty = dzero
      angks = dzero
      do 130 i=1,3
        evec(i) = dzero
        xivec(i) = dzero
        spvec(i) = dzero
 130  continue
      do 150 i=-1,1
        do 140 j=-1,1
          ptz(j,i) = dcmplx(dzero,dzero)
 140    continue
 150  continue

c  initialize atom list data
      do 170 i=1,nattx
 170  iphatx(i) = -1

c  initialize character strings - Josh Kas
      cfname = 'NULL'
      
c  initialize EELS variables !KJ 1-06
      ebeam=dzero
        aconv=dzero
        acoll=dzero
        nqr=0
        nqf=0
        magic=0
        emagic=dzero
        eels=0
        relat=1
        cross=1
        aver=0
      thetax=dzero
        thetay=dzero
        ipmin=1
        ipmax=9
        ipstep=1
        iinput=1  !5/6
       spcol=4
c KJ 

c for ABSOLUTE card  !KJ 3-06
        absolu=0  !KJ 3-06

      return
      end
c  end subroutine iniall
