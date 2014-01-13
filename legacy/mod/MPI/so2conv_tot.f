c///////////////////////////////////////////////////////////////////////
c Distribution:  COMMON 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: COMMON 1.0
c      COMMON 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c     sub-program ffmod9
      program  ffmod9
c     subroutine ffmod9

c     Calculation of S_0^2 
c     written by Luke Campbell 2002
c     modified by Luke Campbell 2005 for new I/O structure

c     INPUT: s02.inp mod6.inp and any set of spectroscopy output files
c            (xmu.dat, chi.dat, chipNNNN.dat, feffNNNN.dat)
c     OUTPUT: specfunct.dat and the input spectroscopy files (overwritten)

      implicit double precision (a-h, o-z)
      character*12 cfname

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}


      call par_begin
      if (worker) go to 400
c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='logso2.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logso2.dat', 'feff')

c     read  s02.inp
      call res02( mso2conv, ispec, ipr6, ipse, ipsk, wsigk, cen, cfname)

      if (mso2conv.eq.1) then 
        call wlog(' Calculating s_0^2 ...')

c       calculate S_0^2
        call so2conv (ispec, ipr6, ipse, ipsk, wsigk, cen, cfname)

        call wlog(' Done with s_0^2.')
      endif

  400 call par_barrier
      call par_end

c     sub-program ffmod9
      stop
c     return

      end
      subroutine sfconv(ekp,mu,gammach,npts2,wpts2,xchi,npts1,wpts1,
     2                spectf,weights,cchi,phase,iasym,icut,intout,omp)
* Convolutes the array xchi (signal) with the array spectf (spectral
* function or asymmetric broadening), where spectf
* has a delta function contribution of magnitude |weights(1)+i weights(2)|
* + weights(3).  Complex phase of weights(1) + i weights(2) is given 
* to overall phase of the otherwise real valued array spectf.
* Input: ekp - Photoelectron energy neglecting collective excitations.
*        mu - Chemical potential, position of edge.
*        gammach - Core hole lifetime.
*        npts2 - Dimension of signal array xchi.
*        wpts2 - Energy grid of signal.
*        xchi - Signal array.
*        npts1 - Dimension of spectral function array spectf.
*        wpts1 - Energy grid of spectral function.
*        weights - Array of weights of components of the spectral
*            function.  Only the delta function weights, weights(1)
*            and weights(3), and the imaginary part of the extrinsic
*            delta function, weights(2) are needed.
*        cchi - The magnitude of the convoluted signal.
*        phase - The phase of the convoluted signal.
*        iasym - set to 1 to include quasiparticle phase as an 
*            asymmetric 1/omega term to the etrinsic satellite
*            rather than as a complex spectral weight.  This is
*            necessary when convoluting with a real valued
*            function or one whose imaginary part is not known.
*        icut - set to 0 to avoid truncating the spectral function
*            at energies where there is insufficient energy to create
*            excitations.
*        intout - Flag to write diagnostic file, the running integration
*            of the convolution.
*        omp - Plasma frequency.
      implicit none
      integer npts1,npts2,i,j,intout,iasym,icut
      double precision wpts2(npts2),xchi(npts2),wpts1(npts1),
     2                 spectf(npts1),weights(8),cchi,phase,ekp,mu,
     3                 gammach,spectf2(npts1),wtmp,qpr,omp
*       spectf2 - The spectral function with a cutoff to prevent 
*           excitation energies higher than the total available 
*           energy in the system (ekp-mu).
*       wtmp - The new quasiparticle weight in going 
*           from spectf to spectf2.
*       qpr - Quasiparticle reduction, from renormalizing 
*           spectral function.
      double precision xrcchi,xicchi,xnorm,dw,w,www,xfact,
     2                 xrchi,am,phasez,efrac,pi,eV,store,
     3                 amp,del,lam
*       xrcchi - Real part of convoluted signal.
*       xicchi - Imaginary part of convoluted signal.
*       xnorm - Total weight of spectral function with cutoff,
*           used to normalize the spectral function.
*       dw - Energy interval.
*       w - Exitation energy (omega)
*       www - Quasiparticle energy, available energy minus excitation energy.
*       xfact - Used to store intermediate calculations.
*       xrchi - Convoluted signal before quasiparticle phase shift.
*       am - Magnitude of quasiparticle weight.
*       phasez - Quasiparticle phase.
*       efrac - Fractional distance from one energy grid point to
*           the next, used in interpolation.
*       pi - Ratio of circumference to diameter of a circle in
*           euclidian geometry.
*       eV - Unit conversion to electron volts.
*       store - Dummy variable used to store a value that will be 
*           overwritten but needs to be used again.
*       amp - "Amplitude", signal value at first data point,
*           used for extrapolating signal to lower energies to avoid
*           artifacts.
*       lam - "Lambda", decay constant for extrapolating signal below the
*           first data point.
*       del - "Delta", energy difference used in extrapolation of signal
*           below the first data point.
      parameter (eV=1.d0/27.21160d0)

      pi=dacos(-1.d0)
      xrcchi=0.d0
      xicchi=0.d0
*      iasym=0
* Extrinsic delta function (quasiparticle) weight.
      if (iasym.eq.1) then
        am=weights(1)
      else
        am=dsqrt(weights(1)**2+weights(2)**2)
      endif
* Quasiparticle phase.
      if (weights(1).ne.0.d0.and.iasym.ne.1) then
        phasez=datan(weights(2)/weights(1))
      else
        phasez=0.d0
      endif
* Reduce quasiparticle weight in accord with the same energy cutoff
* constraints used for the rest of the spectral function.
      if (icut.eq.0) then
        qpr=1.d0
      elseif (ekp-mu.ne.0.d0) then
        qpr=datan2(gammach,mu-ekp)/pi
      else
        qpr=0.5d0
      endif
      wtmp=qpr*(am+weights(3))
      xnorm=wtmp
* Cut off that portion of spectral function with more energy than is
* available in the system.
      do i=1,npts1
        if (i.eq.1) then
          dw=wpts1(2)-wpts1(1)
        elseif (i.eq.npts1) then
          dw=wpts1(npts1)-wpts1(npts1-1)
        else
          dw=(wpts1(i+1)-wpts1(i-1))/2.d0
        endif
        w=wpts1(i)
        www=ekp-w
* Construct weight function, cutoff at chemical potential, width 
* equal to core hole lifetime.
        if (icut.eq.0) then
          xfact=1.d0
        elseif (www-mu.ne.0.d0) then
          xfact=datan2(gammach,mu-www)/pi
        else
          xfact=0.5d0
        endif
* Multiply spectral function by cutoff weight function.
        if (icut.eq.0) then
          spectf2(i)=spectf(i)
        elseif (w.ge.0.d0) then
          spectf2(i)=spectf(i)*xfact
        else
          spectf2(i)=max(0.d0,spectf(i)*xfact)
        endif
        if (iasym.eq.1) then
          spectf2(i)=spectf2(i)-qpr*(weights(2)/(pi*am*dw))
     2      *log(((w+dw/2.d0)**2+(3.d0*dw)**2)
     3         /((w-dw/2.d0)**2+(3.d0*dw)**2))
     4         *exp(-((w)/(2*omp))**2)/2.d0
        endif
* Integration to find total spectral weight.
        xnorm=xnorm+spectf2(i)*dw
      enddo
* Main convolution loop.
      do i=1,npts1
        if (i.eq.1) then
          dw=wpts1(2)-wpts1(1)
        elseif (i.eq.npts1) then
          dw=wpts1(npts1)-wpts1(npts1-1)
        else
          dw=(wpts1(i+1)-wpts1(i-1))/2.d0
        endif
        w=wpts1(i)
        www=ekp-w
        xrchi=0
        if (www.gt.wpts2(npts2)) then
* Extrapolate signal to avoid artifacts.
          xrchi=xchi(npts2)
        elseif (www.le.wpts2(1)) then
* Extrapolate signal to avoid artifacts.
          amp=xchi(1)
          del=mu-wpts2(1)
          lam=del**2/(pi*dabs(amp)*(del**2+gammach**2))
          xrchi=amp*exp(lam*(www-wpts2(1)))
        else
* Interpolate signal onto spectral function energy.
          do j=1,npts2-1
            if (www.gt.wpts2(j).and.www.le.wpts2(j+1)) then
              efrac=(www-wpts2(j))/(wpts2(j+1)-wpts2(j))
              xrchi=xchi(j)+(xchi(j+1)-xchi(j))*efrac
              goto 50
            endif
          enddo
        endif
50      continue
        if ((w+wpts1(i-1))/2.d0.lt.0.d0.and.
     2      (w+wpts1(i+1))/2.d0.ge.0.d0) then
* Add delta function contribution.
          xrcchi=xrcchi+wtmp*xrchi
        endif
* Convolution integration.
        xrcchi=xrcchi+spectf2(i)*dw*xrchi
* Print diagnostic files of running integration if requested.
        if (intout.ne.0) then
          if ((w+wpts1(i-1))/2.d0.lt.0.d0.and.
     2        (w+wpts1(i+1))/2.d0.ge.0.d0) then
            write(28,500) www,xrcchi/xnorm,xrchi,am+weights(3),
     2           spectf2(i)/xnorm
          else
            write(28,500) www,xrcchi/xnorm,xrchi,spectf(i),
     2            spectf2(i)/xnorm
          endif
        endif
      enddo
* Add overall phase equal to quasiparticle phase.
      store=xrcchi
      xrcchi=store*dcos(phasez)-xicchi*dsin(phasez)
      xicchi=xicchi*dcos(phasez)+store*dsin(phasez)
      xrcchi=xrcchi/xnorm
      xicchi=xicchi/xnorm
      cchi=dsqrt(xrcchi**2+xicchi**2)
      phase=datan2(xicchi,xrcchi)

 500  format(1x,5(e12.5,1x))
 501  format(1x,7(e10.3,1x))
      return
      end
      subroutine croots (a,b,c,d,x1,x2,x3,nrroots)
* Solves the cubic polynomial ax^3+bx^2+cx+d=0.
* Returns the number of real roots nrroots and the three roots
* of the cubic polynomial.
      implicit none
      integer nrroots
      double precision a,b,c,d,p,q,r,ar,br,disc,disc2
      complex*16 AA,BB,yr,yi,y1,y2,y3,x1,x2,x3
      x1=0.d0
      x2=0.d0
      x3=0.d0
      if (a.eq.0.d0) then
        if (b.eq.0.d0) then
          if (c.eq.0.d0) then
            nrroots=0
            return
          endif
          nrroots=1
          x1=-d/c
          return
        endif
        disc=c**2-4.d0*d*b
        if (disc.ge.0.d0) then
          nrroots=2
        else
          nrroots=0
        endif
        x1=(-c+sqrt(dcmplx(disc,0.d0)))/(2.d0*b)
        x2=(-c-sqrt(dcmplx(disc,0.d0)))/(2.d0*b)
        return
      endif
      p=b/a
      q=c/a
      r=d/a
      ar=q-p**2/3.d0
      br=(2.d0*p**3-9.d0*p*q)/27.d0+r
      disc=br**2/4.d0+ar**3/27.d0
      if (disc.gt.0.d0) then
        nrroots=1
        disc2=-br/2.d0+sqrt(disc)
        if (disc2.ge.0.d0) then
          AA=dcmplx(disc2**(1.d0/3.d0),0.d0)
        else
          AA=dcmplx(-((-disc2)**(1.d0/3.d0)),0.d0)
        endif
        disc2=-br/2.d0-sqrt(disc)
        if (disc2.ge.0.d0) then
          BB=dcmplx(disc2**(1.d0/3.d0),0.d0)
        else
          BB=dcmplx(-((-disc2)**(1.d0/3.d0)),0.d0)
        endif
      else
        nrroots=3
        if (br.lt.0.d0) then
          AA=dcmplx(-br/2.d0,sqrt(-disc))**(1.d0/3.d0)
          BB=dcmplx(dble(AA),-dimag(AA))
        else
          AA=-dcmplx(br/2.d0,-sqrt(-disc))**(1.d0/3.d0)
          BB=dcmplx(dble(AA),-dimag(AA))
        endif
      endif
      yr=-(AA+BB)/2.d0
      yi=(AA-BB)*sqrt((-3.d0,0.d0))/2.d0
      y1=yr+yi
      y2=yr-yi
      y3=AA+BB
      x1=y1-p/3.d0
      x2=y2-p/3.d0
      x3=y3-p/3.d0
      return
      end
c**********************************************************************
c   This is Steve White's rewrite of Mike Teter's integration routine.
c   Modified by J. Rehr for complex integration.
c   The following is a listing of the arguments in the initial function
c   statement:
c      fn    -- routine requires external function statement in MAIN
c      xmin  -- lower limit
c      xmax  -- upper limit
c      abr   -- absolute tolerable error
c      rlr   -- relative tolerable error
c      nsing -- number of singularities or regions requiring
c                   special attention
c      xsing -- array of locations of singularities or endpoints
c                   of special regions
c      error -- output for routine error messages
c      numcal-- the number of times fn was called
c      maxns -- the maximum number of regions being considered simultaneously
c       function grater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
c       fn declared double precision
c       double precision function grater(fn,xmin,xmax,abr,rlr,
c       fn declared complex*16
c      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1

       double precision function grater(fn,xmin,xmax,abr,rlr,
     1 nsing,xsing,error,numcal,maxns)

       implicit double precision (a-h,o-z)
       parameter (mx=1500)
       dimension xleft(mx),fval(3,mx),dx(3),wt(3)
       dimension wt9(9), xsing(20)
       external fn
        logical atsing
        save dx,wt,wt9
        data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
        data wt/0.277777777777777778  ,0.4444444444444444444  ,
     1                               0.2777777777777777778  /
        data wt9/0.0616938806304841571  ,0.108384229110206161  ,
     1           0.0398463603260281088  ,0.175209035316976464  ,
     2           0.229732989232610220  ,0.175209035316976464  ,
     3           0.0398463603260281088  ,0.108384229110206161  ,
     4           0.0616938806304841571  /
c nstack is the number of different intervals into which the
c integration region is currently divided. The number of regions can
c grow if more accuracy is needed by dividing the right-most region
c into three regions. The number of regions shrinks when the integral
c over the right-most region is accurate enough, in which case that
c integral is added to the total (stored in grater) and the region
c is removed from consideration (and a new region is the right-most).
        nstack=nsing+1
        maxns = nstack
        error=0.
        grater=0.
c The array xleft stores the boundary points of the regions.
c The singular points just govern the initial placement of the regions.
        xleft(1)=xmin
        xleft(nsing+2)=xmax
        if(nsing.gt.0) then
          do 9 j=1,nsing
9           xleft(j+1)=xsing(j)
        endif
c For each region, calculate the function and store at three selected points.
        do 1 jj=1,nstack
          del=xleft(jj+1)-xleft(jj)
c         print*, 'fn call j= ,'
          do 1 j=1,3
c         print*, 'fn call in grater j= ',j
1           fval(j,jj)=fn(xleft(jj)+del*dx(j))
c         print*, 'output of fn call, fval(j,jj)',fval(j,jj)
        numcal = nstack * 3
6       continue
          if(nstack+3.ge.mx) then
            write(*,*) 'TOO MANY REGIONS'
            stop 0006
          endif
c Divide the rightmost region into three subregions.
          del=xleft(nstack+1)-xleft(nstack)
          xleft(nstack+3)=xleft(nstack+1)
          xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
          xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
c The three data points already found for the region become the
c middle data points (number 2 in first index of fval) for each region.
          fval(2,nstack+2)=fval(3,nstack)
          fval(2,nstack+1)=fval(2,nstack)
          fval(2,nstack)=fval(1,nstack)
c Now do the integral over the right-most region in two different ways-
c a three point integral (valu) over each of the three subregions
c and a more accurate nine-point integral (value) over whole region.
c valu is used only for the error estimate.
          icount=0
          value=0.
          valu=0.
          do 3 j=nstack,nstack+2
            del1=xleft(j+1)-xleft(j)
c         print*, 'fn call 2'
            fval(1,j)=fn(xleft(j)+dx(1)*del1)
            fval(3,j)=fn(xleft(j)+dx(3)*del1)
c         print*, 'fn call 2'
            numcal = numcal + 2
            do 5 k=1,3
              icount=icount+1
              value=value+wt9(icount)*fval(k,j)*del
5             valu=valu+fval(k,j)*wt(k)*del1
3         continue
          dif=abs(value-valu)
c If the following condition is true, add in this integral to the total,
c and reduce the number of regions under consideration.
          frac = del / (xmax - xmin)
          atsing = .false.
          if(frac .le. 1.0e-8) atsing = .true.
          if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.
     1       (atsing .and.
     2     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
c The following commented out line is Teeter's old error criterion.
c          if(dif.le.abr.or.dif.le.rlr*abs(value))then
            grater=grater+value
            error=error+abs(dif)
            nstack=nstack-1
c If no more regions, we are done.
            if(nstack.le.0) return
          else
c If the integration is insufficiently accurate, make each of the
c three subregions of the right-most region into regions.
c On next pass the right-most of these is the new current region.
            nstack=nstack+2
            maxns = max(maxns,nstack)
          endif
        go to 6
        end
      subroutine interpsf(npts,epts,wpts,spectf,cspec)
* Interpolates the spectral function calculated on a minimal grid
* to a uniform grid that can be handled by the convolution subroutine.
* input: npts - number of grid points the spectral function will
*               be interpolated on to.
*        epts - energy values of the minimal grid
*        wpts - energy values of the uniform grid
*        spectf - the spectral function on the minimal grid
*        cspec - the spectral function on the uniform grid
      implicit none
      integer npts,nsfpts,i,j
*      parameter (nsfpts=80)
      parameter (nsfpts=110)
      double precision spectf(8,nsfpts),cspec(npts),epts(nsfpts),
     2                 wpts(npts),wmin,wmax,dw,sfhi,sflo,delw
      double precision pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      common /convsf/ pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      wmin=epts(1)
      wmax=epts(nsfpts)
      dw=(wmax-wmin)/(npts-1)
      wpts(1)=wmin
      do i=2,npts
        wpts(i)=wmin+dw*(i-1)
      enddo
      cspec(1)=spectf(2,1)+spectf(5,1)-2.d0*spectf(4,1)
      cspec(npts)=spectf(2,nsfpts)+spectf(5,nsfpts)
     2            -2.d0*spectf(4,nsfpts)
      do i=2,npts-1
        do j=2,nsfpts
          if (wpts(i).ge.epts(j-1).and.wpts(i).lt.epts(j)) then
            delw=wpts(i)-epts(j-1)
            sfhi=spectf(2,j)+spectf(5,j)-2.d0*spectf(4,j)
            sflo=spectf(2,j-1)+spectf(5,j-1)-2.d0*spectf(4,j-1)
            cspec(i)=sflo+(sfhi-sflo)*delw/(epts(j)-epts(j-1))
            goto 10
          endif
        enddo
 10     continue
      enddo
      return
      end
      subroutine mkrmu(xmu,xmu0,rmu,wpts,npts)
* This subroutine does a Cramers-Kronig transform on the array xmu
* and returns an array rmu which is the real part of the analytic
* function whose imaginary part is xmu.  This is needed to get the proper
* phase shift for a convolution with a real function.
      implicit none
      integer npts,i,j
      double precision xmu(npts),xmu0(npts),rmu(npts),wpts(npts)
      double precision dw,pi
      pi=dacos(-1.d0)
      do j=1,npts
        rmu(j)=0.d0
        do i=1,npts
          if (i.eq.1) then
            dw=wpts(2)-wpts(1)
          elseif (i.eq.npts) then
            dw=wpts(npts)-wpts(npts-1)
          else
            dw=(wpts(i+1)-wpts(i-1))/2.d0
          endif
          if (i.ne.j) then
            rmu(j)=rmu(j)+dw*(xmu(i)-xmu0(i))/(wpts(i)-wpts(j))
          endif
        enddo
        rmu(j)=rmu(j)/pi
      enddo
      rmu(20)=(rmu(20)+rmu(21))/2.d0
      rmu(21)=rmu(20)
 500  format(1x,5(e12.5,1x))
      return
      end
      double precision function xmkesat(w)
* Find the extrinsic satellite spectral function.  This is the extrinsic 
* spectral function with the quasiparticle pole subtracted off and the 
* quasiparticle broadening removed.
* input: w - energy (omega)
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       omp - plasma frequency omega_p
*       se - real part of the on shell self energy
*       width - absolute value of the imaginargy part of the 
*            on shell self energy plus the core hole broadening
*       se2 - real part of the self energy at energy w
*       xise - imaginary part of the self energy at energy w
*       z1 - real part of the renormalization constant
*       z1i - imaginary part of the renormalization constant
      implicit none
      integer it1,it2
      double precision w,etot,emain,etothi,etotlo,z1m
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      z1m=sqrt(z1**2+z1i**2)
*      etot=-z1*(width-xise)+z1i*(w+se-se2)
      etot=-(width-xise)
      etot=etot/((w+se-se2)**2+(width-xise)**2)
*      it1=int(xmkesat)
*      it2=int(2*xmkesat)
*      if(it1.eq.it2.and.it1.gt.5) then
*        etot=-sqrt(z1**2+z1i**2)*(width-xise)
*        etothi=etot/((w+omp*1.d-3+se-se2)**2+(width-xise)**2)
*        etotlo=etot/((w-omp*1.d-3+se-se2)**2+(width-xise)**2)
*        etot=(etothi+etotlo)/2.d0
*      endif
*      emain=0.d0
      emain=-z1i/(w*pi*z1m)
      emain=emain*exp(-(w/(2*omp))**2)
      xmkesat=etot/(pi*z1m)-emain
*      it1=int(xmkesat)
*      it2=int(2*xmkesat)
*      if(it1.eq.it2.and.it1.gt.5) then
*        write(6,*) 'nan error: xmkesat'
*        write(6,*) z1,z1i,width,xise,w,se,se2,etot,xmkesat
*      endif
      return
500   format(1x,5(e12.5,1x))
      end

      double precision function xmkgwext(w)
* Find the extrinsic satellite, with full broadening and all quasiparticle
* contributions.
* input: w - energy (omega)
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       se - real part of the on shell self energy
*       se2 - real part of the self energy at energy w
*       xise - imaginary part of the self energy at energy w
      implicit none
      integer it1,it2
      double precision w,etot,emain
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      etot=xise/(pi*((w+se-se2)**2+xise**2))
      emain=0.d0
      xmkgwext=etot-emain
      it1=int(xmkgwext)
      it2=int(2*xmkgwext)
      if(it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'nan error: xmkgwext'
        write(6,*) z1,z1i,width,xise,w,se,se2,etot,xmkgwext
      endif
      return
500   format(1x,5(e12.5,1x))
      end

      double precision function xintxsat(q)
* the integrand of the interference spectral function
* input: q - momentum to be integrated over
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       omp - plasma frequency omega_p
*       ek - bare photoelectron kinetic energy = pk**2/2
*       brd - global broadening parameter to stabilize logarithms
*       ac2 - additional accuracy parameter
*       wp2 - omega prime, an additional energy variable
      implicit none
      integer numcal,maxns,nsing
      double precision wwq,q,xk,vpp2,wdisp,xfact,xloren,
     2       eps,ac2,wp2,dw1,abr,rlr,xsing,error,tol
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      common /ff/ ac2,wp2
      external vpp2,wdisp
      wwq=wdisp(q)
      tol=2.d-1*omp
      if (ek-wp2.ge.0.d0) then
        xk=sqrt(2.d0*(ek-wp2))
        xfact=log(((wwq-q**2/2+xk*q)**2+tol**2)/
     2            ((wwq-q**2/2-xk*q)**2+tol**2))/2.d0
        xloren=ac2/(pi*((wp2-wwq)**2+ac2**2))
        xintxsat=q*vpp2(q)*xloren*xfact/(wwq*xk)
      else
        xk=sqrt(-2.d0*(ek-wp2))
        xfact=datan(xk*q/(wwq-q**2/2))
        xloren=ac2/(pi*((wp2-wwq)**2+ac2**2))
        xintxsat=q*vpp2(q)*xloren*xfact/(wwq*xk)
      endif
      return
      end

      double precision function xintisat(q)
* the integrand of the intrinsic spectral function
* input: q - momentum to be integrated over
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ac2 - additional accuracy parameter
*       wp2 - omega prime, an additional energy variable
      implicit none
      integer numcal,maxns,nsing
      double precision wwq,q,xk,vpp2,wdisp,xfact,xloren,
     3       ac2,wp2,dw1,abr,rlr,xsing,error
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external vpp2,wdisp
      wwq=wdisp(q)
      xloren=ac2/(((wp2-wwq)**2+ac2**2)*pi)
      xintisat=q**2*vpp2(q)*xloren/wwq**2
      return
      end

      double precision function xmkxsat(w,width)
* generates the interference spectral function
* input: w - energy (omega)
*        width - additional broadening
* input from common blocks
*       omp - plasma frequency omega_p
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp2 - omega prime, an additional energy variable
      implicit none
      integer nsing,numcal,maxns,i
      double precision xintxsat,grater,abr,rlr,xsing,error,
     2       qk,wp,acc,pi,ef,xmu,qf,omp,ompl,wt,w,pk,width,
     3       ekp,ek,ac2,wp2,q2,qwidth,qmin,qmax,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external xintxsat,grater
      wp2=w
      ac2=width
      rlr=acc
      abr=omp*acc
      nsing=0
      q2=sqrt(max(2*(w-ompl),ac2))
      qwidth=10.d0*ac2/q2
      qmin=max(0.d0,q2-qwidth)
      qmax=q2+qwidth
*      do i=1,1000
*      enddo
*      write(6,*) 'xmkxsat'
      xmkxsat=grater(xintxsat,qmin,q2,
     2             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkxsat=xmkxsat+grater(xintxsat,q2,qmax,
     2             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkxsat=xmkxsat/(2.d0*pi)**2
      return
      end

      double precision function xmkisat(w,width)
* generates the intrinsic spectral function
* input: w - energy (omega)
*        width - additional broadening
* input from common blocks
*       omp - plasma frequency omega_p
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp2 - omega prime, an additional energy variable
      implicit none
      integer nsing,numcal,maxns
      double precision xintisat,grater,abr,rlr,xsing,error,
     2       qk,wp,acc,pi,ef,xmu,qf,omp,ompl,wt,z,w,pk,width,
     3       ekp,ek,ac2,wp2,q2,qmin,qmax,qwidth,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external xintisat,grater
      wp2=w
      ac2=width
      rlr=acc
      abr=omp*acc
      nsing=0
      if ((w-ompl).gt.ac2) then
        q2=dsqrt(2.d0*(w-ompl))
      else
        q2=dsqrt(2.d0*ac2)
      endif
      qwidth=10.d0*min(q2,ac2/q2)
      qmin=max(0.d0,q2-qwidth)
      qmax=q2+qwidth
      xmkisat=grater(xintisat,0.d0,q2,
     2             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkisat=xmkisat+grater(xintisat,q2,qmax,
     2             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkisat=xmkisat/(2.d0*pi**2)
      return
      end

*      double precision function xmkisat(w,width)
**     generates the intrinsic spectral function
*      implicit none
*      double precision w,wdisp,qdisp,dwdq,q,vpp2,width
*      double precision pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
*      common /convsf/ pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
*      external wdisp,dwdq,qdisp,vpp2
*      xmkisat=0.d0
*      if (w.gt.omp) then
*        q=qdisp(w)
*        xmkisat=q**2*vpp2(q)/(2.d0*pi**2*w**2*dwdq(q))
*      endif
*      return
*      end
*
*      double precision function xmkxsat(w)
**     generates the interference spectral function
*      implicit none
*      double precision w,q,wdisp,qdisp,dwdq,vpp2,xk,xfact,
*     2                 eta
*      complex*16 cfac,coni
*      parameter (coni=(0,1))
*      double precision pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
*      common /convsf/ pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
*      external wdisp,dwdq,qdisp,vpp2
*      xmkxsat=0.d0
*      eta=1.d-2*omp
*      if (w.gt.omp.and.ek.gt.w) then
*        q=qdisp(w)
*        xk=dsqrt(2.d0*(ek-w))
*        xfact=log(((w-q**2/2+xk*q)**2+(omp*acc)**2)/
*     2            ((w-q**2/2-xk*q)**2+(omp*acc)**2))/2.d0
**        cfac=log((w-q**2/2+xk*q-coni*eta)/(w-q**2/2-xk*q-coni*eta))
**        xfact=dble(cfac)
*        xmkxsat=q*vpp2(q)*xfact/((2.d0*pi)**2*xk*w*dwdq(q))
*      endif
*      return
*      end

      double precision function xmkak(w)
* This function returns the interference contribution to the 
* quasiparticle.
* input: w - energy (omega)
* input from common blocks
*       omp - plasma frequency omega_p
*       ek - bare photoelectron kinetic energy = pk**2/2
*       acc - global accuracy parameter 
*       wmax - highest allowed energy
* common block control of subprograms
*       pkq - momentum variable used in the integrand but kept 
*             fixed throughout the integration
      implicit none
      integer i,j,jj,nsing,numcal,maxns,it1,it2
      double precision w,qmax,xintak,
     2                 abr,rlr,error,xsing,grater
      integer nnpts
      double precision pkq
      common /funct2/ pkq
      double precision wmin,wmax,wmin1,wmax1
      common /limits/ wmin,wmax,wmin1,wmax1
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external xintak
      if (w.gt.0.d0) then
        rlr=acc
        abr=dsqrt(omp)*acc
        qmax=dsqrt(2*wmax)
        pkq=dsqrt(2*ek)
        nsing=0
*        write(6,*) 'xmkak'
        xmkak=grater(xintak,abr,qmax,abr,rlr,nsing,xsing,
     2                 error,numcal,maxns)
      else
        xmkak=0.d0
      endif
      return
      end

      double precision function xintak(q)
* Integrand for the function xmkak.
* input: q - momentum to be integrated over.
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       omp - plasma frequency omega_p
*       pkq - momentum variable used in the integrand but kept 
*             fixed throughout the integration
      implicit none
      double precision wq,q,pkq,xlog,wdisp,vpp2,eps
      common /funct2/ pkq
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external vpp2,wdisp
      wq=wdisp(q)
      eps=1.d-1
      xlog=(wq+q**2/2.d0+pkq*q)**2+(ompl*eps)**2
      xlog=xlog/((wq+q**2/2.d0-pkq*q)**2+(ompl*eps)**2)
      xlog=log(xlog)/2.d0
      xintak=q*vpp2(q)*xlog/(wq*pkq*4.d0*pi**2)
      return
      end
c  Computes spectral function window centered on resonance.  
c  Spectral function "deconvoluted" to remove broadening of main peak
      subroutine mkspectf(rs,pk,gammach,xreduc,wpts,spectf,
     2                    weights,isattype,npl,nplmax,
     3                    plengy,plwt,plbrd,brpole) 
* input: rs - radius that contains one electron
*        pk - photoelectron momentum
*        gammach - core hole lifetime broadening
*        xreduc - add-hoc interference reduction
*        isattype - approximation to use when computing satellite.
*        npl - number of poles in epsilon^{-1}
*        nplmax - maximum number of poles for dimensioning arrays
*        plengy - energy of each pole in epsilon^{-1}
*        plwt - weight of each pole in epsilon^{-1}
*        plbrd - broadening of each pole in epsilon^{-1}
*        brpole - true if poles are to be calculated with broadening
* output: wpts - array of energy points on which the spectral
*                function is specified.  array is of length npts
*         spectf - array of values of the spectral function.
*                array is of length npts
*         weights - array containing the total spectral weight of 
*                various contributions to the spectral function
*         weights(1) - extrinsic quasiparticle weight
*         weights(2) - asymmetry of quasiparticle
*         weights(3) - interference quasiparticle weight
*         weights(4) - extrinsic satellite weight
*         weights(5) - interference satellite weight
*         weights(6) - intrinsic satellite weight
*         weights(7) - weight of the extrinsic satellite clipped 
*                      to include only that part in the satellite region
*         weights(8) - weight of the extrinsic satellite clipped
*                      to include that part in the vicinity of the
*                      quasiparticle
      implicit none
      logical brpole
      integer npts,i,ii,j,jj,k,iemax,iixmax,ishift,isetold,
     2        iset,iset2,iswitch,isattype,iqph,iqpl,
     3        isetd,isetd2,iswd,iswd2,jwidth,ipl
*             npts - number of grid points for spectral function
*             i,ii,j,jj,k - counters and dummy variables
*             iixmax - grid point of highest value of the
*                   interference spectral function
*             ishift - number of grid points to shift the 
*                   extrinsic spectral function
*             iset,iset2,isetd,isetd2,isetold - keeping track
*                   of trigger conditions to separate the 
*                   extrinsic satellite into "satellite" and
*                   "main peak" regions
*             iswd,iswd2,iswitch - keeping track of grid points
*                   at which the above triggering conditions are met
*             iqph,iqpl - array points that bound the quasiparticle
*             jwidth - extra region tacked onto integration for 
*                   lorentzian broadening of the self energy
*             ipl - counter for loops over poles
*      parameter(npts=80)
      parameter(npts=112)
      integer npl,nplmax
      double precision rs,pk,gammach,xreduc,
     2                 spectf(8,npts),wpts(npts)
      double precision plengy(nplmax),plwt(nplmax),plbrd(nplmax)
      double precision conc,xa,expa,ak,aangstrom,eV,dsat,d2sat,
     2                 specttmp(npts),wswitch,swidth,wd,wd2,
     3                 dsatold,d2satold,dsath,wsearch,esfhi,esflo
*             conc - electron concentration
*             xa - dimensionless coupling constant of the electron gas
*             expa - the base of the natural logarithm raised to the
*                  power xa (e**xa)
*             ak - interference contribution to the quasiparticle
*             aangstrom - convert from angstrom to Bohr units of length
*             eV - convert from electron volts to Hartree units of energy
*             dsat,d2sat,dsatold,d2satold,dsath - finite difference 
*                  derivatives of the extrinsic satellite, to find
*                  the triggers to separate the "main peak" from the
*                  "satellite" structure of the extrinsic satellite 
*                  spectral function
*             specttmp - interpolated spectral function for shifting
*                  the extrinsic spectral function without numerical
*                  jitter
*             wsearch - energy used for finding the shift in the
*                  extrinsic spectral function
*             wswitch,wd,wdold - keep track of energy at which the trigger
*                  for separating the extrinsic spectral function
*                  is met
*             esfhi,esflo - used to bracket the high point of the
*                  extrinsic spectral function for its shift
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      double precision xfact1,xfact2,wemax,wshift,wshift2,
     3                 whi,wlo,whi2,wlo2,ehi,elo,emax,ehi1,ehi2,
     4                 elo1,elo2,delta
*             xfact1,xfact2 - dummy variables for keeping track of 
*                  lengthly computations
*             wemax - the energy at the maximum value of the 
*                  extrinsic spectral function
*             wshift,wshift2 - amounts to shift the extrinsic 
*                  spectral function
*             whi,wlo,whi2,wlo2 - used to bracket the search region for the
*                  maximum of the extrinsic spectral function
*             ehi,elo,ehi1,ehi2,elo1,elo2 - values of the extrinsic 
*                  spectral function at the above bracket points
*             emax - the maximum value of the extrinsic spectral function
*             delta - a small value
      double precision sefr(npts),sefi(npts),
     2                 sef2r(npts),sef2i(npts),sumr,sumi,brpl,
     3                 wh,wl,w2,wh2,wl2,dw2
*             sefr - array of the real part of the self energy 
*                  values on the energy grid wpts
*             sefi - array of the imaginary part of the self energy 
*                  values on the energy grid wpts
*             sef2r,sef2i - self energy after lorentzian broadening
*             sumr,sumi - running integration counters
*             brpl - the broadening of a given pole in epsilon^{-1}
*             wh,wl,wh2,wl2 - high and low values of energy intervals,
*                  used for integration bounds over that interval
*             w2 - an energy variable
*             dw2 - an energy interval
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of ipl'th pole in epsilon^{-1}
*       wt - weight of ipl'th pole in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       qpk - photoelectron momentum
*       acc - global accuracy parameter
*       brd - width of pole in epsilon^(-1)
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
      double precision wmin,wmax,wmin1,wmax1,dw,dw1,w,wlim(0:npts)
      common /limits/ wmin,wmax,wmin1,wmax1
*       wmin,wmax,wmin1,wmax1 - a bunch of extreme values of the energy
*       dw,dw1 - energy intervals
*       w - energy
*       wlim - energy values separating the energy gridpoints
      double precision se,xibeta,ce,width,xnn,xaa,se2,xise,
     2                 z0,z0i,zm,z1,z1i,z1m,xrz,xiz,x2a,x2n,
     3                 sef0,se0,yr,yi,zzr,zzi,zzm,qpengy,qpwidth
*       se - real part of the on shell self energy
*       xibeta - imaginary part of the on shell self energy, primarily 
*             used only to find width
*       ce - core energy
*       width - quasiparticle broadening
*       xnn - on shell energy derivative of the real part of the self
*             energy
*       xaa - on shell energy derivative of the imaginary part of the self
*             energy
*       se2 - real part of the self energy calculated at energy w
*       xise - imaginary part of the self energy calculated at energy w
*       z0 - approximation for the real part of the renormalization constant
*       z0i - approximation for the imaginary part of the 
*             renormalization constant
*       zm - magnitude of the above approximation for the 
*             renormalization constant
*       z1,z1i,z1m - a better approximation for the renormalization constant 
*       xrz,xiz - variables for intermediate steps in the calculation of
*             the renormalization constant
*       x2n,x2a - the real and imaginary parts, respectively, of the
*             second energy derivative of the on shell self energy.
*             These variables are not currently used in this version
*             of the code, but if needed, here they are.
*       sef0 - self energy at the fermi level
*       se0 - approximation for the self energy
*       yr - real part of second derivative of self energy
*       yi - imaginary part of second derivative of self energy
*       zzr,zzi,zzm - refined approximation for renormalization constant
*       qpengy - refined approximation for quasiparticle energy
*       qpwidth - refined approximation for quasiparticle width
      double complex zzc,yyc,xxc,discr,epc
*       zzc - complex renormalization constant
*       yyc - second derivative of self energy
*       xxc - 1-derivative of self energy
*       discr - square root of complex discriminant
*       epc - complex energy of quasiparticle pole
      double precision sxnn,sxaa,sse,sxise
*       for sums over the poles in epsilon^{-1} when computing the 
*       self energy and its derivatives
      common /energies/ se,ce,width,z1,z1i,se2,xise
      double precision xsat,xmain,emain,esat,xisat,esat2
*       the values of the interference satellite, interference 
*       quasiparticle, extrinsic quasiparticle, extrinsic satellite,
*       intrinsic satellite, and quasiparticle asymmetry spectral 
*       functions, respectively
      double precision xwidth, xiwidth
*       xwidth - broadening of the interference satellite
*       xiwidth - broadening of the intrinsic satellite
      double precision wtemain,wtesat,wtxmain,wtxsat,wtisat,
     2                 wtmesat,wtmemain,wtesat2,weights(8)
*       the total spectral weights of the extrinsic quasiparticle,
*       extrinsic satellite, interference quasiparticle, 
*       interference satellite, intrinsic satellite, clipped 
*       "satellite" part of the extrinsic satellite, clipped
*       "quasiparticle" part of the extrinsic satellite, and quasipartcle
*       asymmetry terms in the spectral function, respectively, in 
*       addition to the array of spectral weights
      double precision swtcorr,satwt,swtfac
*       swtcorr - the weight of the negative regions of the satellite
*       satwt - the satellite weight (neglecting factors of expa, as above)
*       swtfac - renormalization factor to keep satellite weight the 
*          same after negative parts are chopped off
      double precision grater,abr,rlr,xsing,error
      integer nsing,numcal,maxns
*       Grater is an integration function.  See its description for
*       the meaning of these variables.
      double precision xmkesat,beta,xmkak,xmkxsat,xmkisat,
     2                 xmkgwext,xmkssat,xmkpsat,xmkxg,xmkig
*       These are functions called in the execution of this subroutine.
      integer iwrite,jcount
      common /flag/ iwrite,jcount
*       These are used to flag a particular run of this subroutine
*       to write intermediate results to a file.
      integer lowq
*       lowq - If not equal to zero, calculate contributions
*          to self energy from below Fermi level.
      common /belowqf/ lowq
      external beta,grater,xmkesat,xmkxsat,xmkisat,xmkak,xmkgwext,
     2         xmkssat,xmkxg,xmkig
      double precision w3,q3,q4
      common /fqso2/ q3
      common /fw/ w
c     Josh - comment out extraneous externals
c      double precision wintr1,winti1,wintr2,winti2,xdumr,xdumi
c      external wintr1,winti1,wintr2,winti2
c      double precision srrpaw1,sirpaw1,rpadsfr,rpadsfi
c      external srrpaw1,sirpaw1,rpadsfr,rpadsfi
c      double precision rpapolr,rpapoli
c     external rpapolr,rpapoli
c      double precision xtestr,xtesti,qq,wplas,dq
c      external xtestr,xtesti,wplas
c      double precision rpasigr,rpasigi
c      external rpasigr,rpasigi
c     End Josh
      double precision exchange
      external exchange
      integer iw
      integer ijkwrite
      common /morewrite/ ijkwrite
* test functions

* compute initial values
      ijkwrite=0
      acc=1.d-4
      pi=dacos(-1.d0)
      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
      ef=qf*qf/2.d0
      conc=3.d0/(4.d0*pi*(rs**3))
      omp=dsqrt(4.d0*pi*conc)
      qpk=qf
      ek=ef
      ekp=ef
      adisp=2.d0*ef/3.d0
      sef0=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(0.d0,sse,sxise)
        else
          call renergies(0.d0,sse)
        endif
*        call ppset(rs,pi,qf,ef,omp)
        sef0=sef0+sse*wt
      enddo
      sef0=sef0+exchange(qf)
*      call ppset(rs,pi,qf,ef,omp)
      xmu=ef+sef0
      ekp=xmu

* first estimate for energies and self energies
      qpk=pk
      ek=pk*pk/2.d0
      ekp=ek
      se0=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(0.d0,sse,sxise)
          call dbrsigma(0.d0,sxnn,sxaa)
        else
          call renergies(0.d0,sse)
          call drenergies(0.d0,sxnn)
          call dienergies(0.d0,sxaa)
        endif
*        call ppset(rs,pi,qf,ef,omp)
        se0=se0+sse*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se0=se0+exchange(pk)
*      call ppset(rs,pi,qf,ef,omp)
      xrz=1.d0-xnn
      xiz=-xaa
      z0=xrz/(xrz**2+xiz**2)
      z0i=-xiz/(xrz**2+xiz**2)
      zm=dsqrt(z0**2+z0i**2)
      ekp=ek+sef0+z0*(se0-sef0)

* refined estimate for self energies, using the self energy at the 
* Fermi level to increase self consistency
      se=0.d0
      xibeta=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(-sef0,sse,sxise)
          call dbrsigma(-sef0,sxnn,sxaa)
        else
          call renergies(-sef0,sse)
          call xienergies(-sef0,sxise)
          call drenergies(-sef0,sxnn)
          call dienergies(-sef0,sxaa)
        endif
        se=se+sse*wt
        xibeta=xibeta+sxise*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se=se+exchange(pk)
      width=dabs(xibeta)+gammach

* calculate renormalization constant
      xa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        xa=xa+3*wt*(omp/ompl)**2/(8*dsqrt(2.d0*ompl))
      enddo
      expa=exp(-xa)
      xrz=1.d0-xnn
      xiz=-xaa
      z1=xrz/(xrz**2+xiz**2)
      z1i=-xiz/(xrz**2+xiz**2)
      z1m=sqrt(z1**2+z1i**2)

* Find energy gridpoints
      dw=omp/30.d0
      iqph=54
      iqpl=53
      wpts(iqph)=dw*1.d-2
      wpts(iqpl)=-dw*1.d-2
      wpts(iqph+1)=dw*2.d-2
      wpts(iqpl-1)=-dw*2.d-2
      do i=1,30
        wpts(i+1+iqph)=i*dw
        wpts(iqpl-1-i)=-i*dw
      enddo
      do i=1,3
        wpts(i+31+iqph)=wpts(31+iqph)+i*dw
        wpts(iqpl-31-i)=wpts(iqpl-31)-i*dw
      enddo
      do i=1,3
        wpts(i+34+iqph)=wpts(34+iqph)+(2*i)*dw
        wpts(iqpl-34-i)=wpts(iqpl-33)-(2*i)*dw
      enddo
      do i=1,3
        wpts(i+37+iqph)=wpts(37+iqph)+(4*i)*dw
        wpts(iqpl-37-i)=wpts(iqpl-36)-(4*i)*dw
      enddo
      do i=1,12
        wpts(i+40+iqph)=wpts(40+iqph)+(10*i)*dw
        wpts(iqpl-40-i)=wpts(iqpl-39)-(10*i)*dw
      enddo
      do i=1,6
        wpts(i+52+iqph)=wpts(52+iqph)+(30*i)*dw
      enddo
      do i=1,npts-1
        wlim(i)=(wpts(i)+wpts(i+1))/2.d0
      enddo
      wlim(0)=2.d0*wpts(1)-wpts(2)
      wlim(npts)=2.d0*wpts(npts)-wpts(npts-1)
      wmin=wpts(1)
      wmax=wpts(npts)
      wmax1=ekp+wmax
      wmin1=ekp+wmin

* compute self energies on energy gridpoints
      do i=1,npts
        sefr(i)=0.d0
        sefi(i)=0.d0
      enddo
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        do i=1,npts
          w=wpts(i)+ekp
          dw1=wlim(i)-wlim(i-1)
*          dw1=wpts(i)-wpts(i-1)
  
*          call renergies(w-ekp-sef0,sse)
*          call xienergies(w-ekp-sef0,sxise)
          if (brpole) then
            call brsigma(w-ekp-se,sse,sxise)
          else
            call renergies(w-ekp-se,sse)
            call xienergies(w-ekp-se,sxise)
          endif
          sefr(i)=sefr(i)+sse*wt
          sefi(i)=sefi(i)+sxise*wt
        enddo
*        call ppset(rs,pi,qf,ef,omp)
      enddo
      se2=0.d0
      xise=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(-se,sse,sxise)
          call dbrsigma(-se,sxnn,sxaa)
        else
          call renergies(-se,sse)
          call xienergies(-se,sxise)
          call drenergies(-se,sxnn)
          call dienergies(-se,sxaa)
        endif
        se2=se2+sse*wt
        xise=xise+sxise*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se=se2+exchange(pk)
      width=dabs(xise)+gammach


      do i=1,npts
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        sefr(i)=sefr(i)+exchange(pk)
        sefi(i)=dabs(sefi(i))+gammach
        if (iwrite.eq.jcount) write(24,500) wpts(i),sefr(i),sefi(i)
      enddo
*      se=(sefr(iqph)+sefr(iqpl))/2
*      width=(sefi(iqph)+sefi(iqpl))/2
*      xnn=(sefr(iqph)-sefr(iqpl))/(wpts(iqph)-wpts(iqpl))
*      xaa=(sefi(iqph)-sefi(iqpl))/(wpts(iqph)-wpts(iqpl))
      xrz=1.d0-xnn
      xiz=-xaa
      z1=xrz/(xrz**2+xiz**2)
      z1i=-xiz/(xrz**2+xiz**2)
      z1m=sqrt(z1**2+z1i**2)
      qpengy=ekp+width*z1i
      qpwidth=width*z1
      zzr=z1
      zzi=z1i
      zzm=z1m

*      write(68,*) pk/qf,se,width,z1,z1i

* correct for endpoint effects
      ak=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        ak=ak+xmkak(ekp)*xreduc*wt
*        call ppset(rs,pi,qf,ef,omp)
      enddo
      wtemain=(datan(wlim(0)/width)+pi/2.d0)/pi
     2        +(pi/2.d0-datan(wlim(npts)/width))/pi
      wtxmain=2.d0*wtemain*zm*z1*ak
      wtemain=wtemain*z1*expa
      wtesat=0.d0
      wtesat2=0.d0
      wtxsat=0.d0
      wtisat=0.d0
      do i=1,npts
        do ii=1,8
          spectf(ii,i)=0.d0
        enddo
      enddo

* compute spectral function
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)

        se2=sefr(i)
        xise=sefi(i)

* This form for the extrinsic quasiparticle is calculated at the
* correct quasiparticle pole, for much better cancelation of the
* quasiparticle structure from the xmkgwext function, below.
        emain=z1*(datan((wlim(i)-qpengy+ekp)/qpwidth)
     2        -datan((wlim(i-1)-qpengy+ekp)/qpwidth))/(pi*dw1)
     3        -z1i*dlog((qpwidth**2+(wlim(i)-qpengy+ekp)**2)/
     4         (qpwidth**2+(wlim(i-1)-qpengy+ekp)**2))/(2.d0*pi*dw1)
     5         *exp(-((w-qpengy)/(2*omp))**2)
        xmain=2.d0*zm*ak*emain
        wtemain=wtemain+emain*expa*dw1
        wtxmain=wtxmain+xmain*expa*dw1
* Depending on the value of isattype, different approximations
* can be used for the satellite.
        if (isattype.eq.1) then
          esat=xmkgwext(w-ekp)-emain
        elseif (isattype.eq.2) then
          esat=(xise-width-(w-ekp)*xaa)/(pi*(w-ekp)**2)
        elseif (isattype.eq.3) then
* generate full extrinsic spectral function - do not use for spectroscopy!
          esat=xmkgwext(w-ekp)
        else
          esat=xmkesat(w-ekp)
        endif
        xsat=0.d0
        xisat=0.d0
        do ipl=1,npl
          call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
          xwidth=max(5.d0*dw,brd)
          xiwidth=max(2.d0*dw,brd)
          if (isattype.eq.3) then
            xwidth=xwidth+width
            xiwidth=xiwidth+width
          endif
          xsat=xsat+xmkxsat(w-ekp,xwidth)*xreduc*wt
          xisat=xisat+xmkisat(w-ekp,xiwidth)*wt
*          call ppset(rs,pi,qf,ef,omp)
        enddo
        wtxsat=wtxsat+xsat*dw1*expa
        wtesat=wtesat+esat*dw1*expa
        wtisat=wtisat+xisat*dw1*expa
        spectf(1,i)=emain
        spectf(2,i)=esat
        spectf(3,i)=xmain
        spectf(4,i)=xsat
        spectf(5,i)=xisat
        spectf(6,i)=(esat+xisat-2.d0*xsat)
        if (isattype.eq.3) then
          spectf(6,i)=spectf(6,i)+xmain
        endif
      enddo
      spectf(2,iqpl)=(spectf(2,iqpl)+spectf(2,iqph))/2
      spectf(2,iqph)=spectf(2,iqpl)

*     separate quasiparticle-like structure from satellite in non-delta
*     function part of extrinsic spectral function.
      iset=0
      iset2=0
      dsatold=0.d0
      d2satold=0.d0
      isetd=0
      isetd2=0
      do ii=2,npts-1
        isetold=iset2
        i=npts+1-ii
        w=wpts(i)+ekp
        dsat=(spectf(2,i)-spectf(2,i-1))/(wpts(i)-wpts(i-1))
        dsath=(spectf(2,i+1)-spectf(2,i))/(wpts(i+1)-wpts(i))
        d2sat=(dsath-dsat)/(wlim(i)-wlim(i-1))
        if (dsat.gt.0.d0.and.spectf(2,i).gt.0.d0
     2      .and.iset.eq.0.d0) iemax=i
        if (spectf(5,i).lt.spectf(5,i+1)
     2      .and.spectf(5,i+1).gt.spectf(5,i+2)) iixmax=i
        if (dsat.gt.0.d0.and.spectf(2,i).gt.0.d0) iset=1
        if (beta(0.d0).gt.0.d0) then
          if (dsat.lt.0.d0.and.dsatold.ge.0.d0
     2        .and.iset.eq.1.and.isetd.eq.0) then
             isetd=1
             iswd=i
             wd=w
          endif
          if (d2sat.gt.0.d0.and.d2satold.le.0.d0
     2        .and.iset.eq.1.and.isetd2.eq.0) then
             isetd2=1
             iswd2=i
             wd2=w
          endif
        else
          if (dsat.lt.0.d0.and.dsatold.ge.0.d0
     2        .and.iset.eq.1.and.isetd.eq.0.and.w.gt.0.d0) then
             isetd=1
             iswd=i
             wd=w
          endif
          if (d2sat.gt.0.d0.and.d2satold.le.0.d0
     2        .and.iset.eq.1.and.isetd2.eq.0) then
             isetd2=1
             iswd2=i
             wd2=w
          endif
        endif
      enddo
      if (isetd.eq.1) then
        iswitch=iswd
        wswitch=wd
      else
        iswitch=iswd2
        wswitch=wd2
      endif
*      swidth=(omp/4.d0)*spectf(2,iswitch)/spectf(2,iemax)
      swidth=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        if (swidth.gt.0.d0) then
          spectf(7,i)=spectf(2,i)/(1.d0+exp((w-wswitch)/swidth))
          spectf(8,i)=spectf(2,i)/(1.d0+exp((wswitch-w)/swidth))
        elseif (i.ge.iswitch) then
          spectf(8,i)=spectf(2,i)
        else
          spectf(7,i)=spectf(2,i)
        endif
      enddo

** find the energy at which the extrinsic satellite reaches its 
** maximum value.
*      wemax=wpts(iemax)+ekp
*      whi=wpts(iemax+1)+ekp
*      wlo=wpts(iemax-1)+ekp
*      emax=spectf(2,iemax)
*      do i=1,8
*        whi2=(whi+wemax)/2.d0
*        wlo2=(wlo+wemax)/2.d0
*        if (whi2-ekp.ne.0.d0) then
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          ehi=xmkesat(whi2-ekp)
*        else
*          delta=(whi2-wemax)/1.d3
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          ehi1=xmkesat(whi2-ekp+delta)
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          ehi2=xmkesat(whi2-ekp-delta)
*          ehi=(ehi1+ehi2)/2.d0
*        endif
*        if (wlo2-ekp.ne.0.d0) then
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          elo=xmkesat(wlo2-ekp)
*        else
*          delta=(wemax-wlo2)/1.d3
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          elo1=xmkesat(wlo2-ekp+delta)
*          se2=0.d0
*          xise=0.d0
*          do ipl=1,npl
*            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
*            call brsigma(whi2-ekp-sef0,sse,sxise)
**            call ppset(rs,pi,qf,ef,omp)
*            se2=se2+sse*wt
*            xise=xise+sxise*wt
*          enddo
*          se2=se2+exchange(pk)
*          xise=dabs(xise)+gammach
*          elo2=xmkesat(wlo2-ekp-delta)
*          elo=(elo1+elo2)/2.d0
*        endif
*        if (ehi.gt.emax.and.ehi.gt.elo) then
*          wlo=wemax
*          wemax=whi2
*          emax=ehi
*        elseif(elo.gt.emax.and.elo.ge.ehi) then
*          whi=wemax
*          wemax=wlo2
*          emax=elo
*        else
*          whi=whi2
*          wlo=wlo2
*        endif
*      enddo
*      wemax=wemax-ekp
*
** if the extrinsic spectral function peaks beyond the plasma frequency,
** shift the spectral function so it peaks at the plasma frequency.
**      if(wemax.lt.omp) then
*      if(.false.) then
*        wshift=omp-wemax
*        ishift=53-iemax
**       here, wpts(53)=omp
*        wshift2=wpts(iemax)-wemax
*        if (wshift.lt.0.d0) then
*          ishift=ishift-1
*          wshift2=wpts(iemax+1)-wemax
*        endif
**        ishift=int(wshift/dw)
**        wshift2=wshift-ishift*dw
*        do i=1,ishift
*          specttmp(i)=0.d0
*        enddo
*        do i=1,npts
*          wsearch=wpts(i)-wshift
*          if (wpts(1).gt.wsearch) then
*            specttmp(i)=0.d0
*          else
*            do j=1,npts
*              if (wpts(j).gt.wsearch) then
*                whi=wpts(j)
*                wlo=wpts(j-1)
*                esfhi=spectf(2,j)
*                esflo=spectf(2,j-1)
*                wshift2=wpts(j)-wsearch
*                specttmp(i)=esflo+(esfhi-esflo)*wshift2/(whi-wlo)
*                goto 100
*              endif
*            enddo
*          endif
*100       continue
*        enddo
*      endif


* Find the weights for the clipped parts of the extrinsic satellite
      wtmesat=0.d0
      wtmemain=0.d0
      do i=1,npts
        wtmesat=wtmesat+spectf(8,i)*expa*dw
        wtmemain=wtmemain+spectf(7,i)*expa*dw
      enddo
** Weight the satellite terms by the power of the renormalization
** constant reflecting how many powers of the extrinsic Green's function
** appear in the original expressions
*      wtxsat=wtxsat*(zm+wtmemain)**2
*      wtisat=wtisat*(zm+wtmemain)
      do i=1,npts
*        spectf(4,i)=spectf(4,i)*(zm+wtmemain)**2
*        spectf(5,i)=spectf(5,i)*(zm+wtmemain)
        spectf(6,i)=spectf(2,i)-2.d0*spectf(4,i)+spectf(5,i)
      enddo
      
* eliminate regions of negative spectral weight, keep track of 
* the total satellite weight before correction and weight of negative
* regions removed for later renormalization.
      swtcorr=0.d0
      satwt=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)
        satwt=satwt+spectf(6,i)*dw1
        if (spectf(6,i).lt.0.d0) then
          swtcorr=swtcorr+spectf(6,i)*dw1
          spectf(6,i)=0.d0
          spectf(4,i)=(spectf(2,i)+spectf(5,i))/2.d0
        endif
      enddo
* renormalize to keep satellite weight the same.
      swtfac=satwt/(satwt-swtcorr)
      if (swtfac.lt.0.d0) swtfac=0.d0
      wtxsat=0.d0
      wtisat=0.d0
      wtesat=0.d0
      wtmesat=0.d0
      wtmemain=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)
        spectf(4,i)=(spectf(2,i)+spectf(5,i)-spectf(6,i)*swtfac)/2.d0
        wtxsat=wtxsat+spectf(4,i)*dw1*expa
        wtesat=wtesat+spectf(2,i)*dw1*expa
        wtisat=wtisat+spectf(5,i)*dw1*expa
        wtmesat=wtmesat+spectf(8,i)*expa*dw
        wtmemain=wtmemain+spectf(7,i)*expa*dw
      enddo
        
* Write out array of weights
      weights(1)=z1*expa
      write(66,*) z1,expa
      weights(2)=z1i*expa
      weights(3)=2.d0*z1*zm*ak*xreduc*expa
      weights(4)=wtesat
      weights(5)=wtxsat
      weights(6)=wtisat
      weights(7)=wtmesat
      weights(8)=wtmemain

*      write(6,*) xxc,yyc,discr,abs(discr),
*     2           -2.d0*width*(0,1)*yyc,zzc,abs(zzc),z1,z1i
*      z1=zzr
*      z1i=zzi
*      z1m=zzm

      if (iwrite.eq.jcount) then
        do i=1,npts
          write(12,500) wpts(i),spectf(1,i),spectf(3,i),
     2                  spectf(7,i),spectf(8,i)
          write(13,500) wpts(i),spectf(2,i),spectf(4,i),
     2                  spectf(5,i),
     3                  spectf(2,i)+spectf(5,i)-2.d0*spectf(4,i)
*     2                  spectf(5,i),spectf(6,i)
        enddo
      endif

      return
 500  format(1x,5(e12.5,1x))
 700  format(1x,7(f10.5,1x))
      end
      double precision function wdisp(q)
* dispertion relation
* input: q - momentum or wavenumber
* input from common blocks
*       ompl - zero q energy of mode
*       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      wdisp=sqrt(ompl**2+adisp*q**2+q**4/4.d0)
      return
      end

      double precision function dwdq(q)
* the derivative of the dispertion relation 
* with respect to the plasmon momentum q
* input: q - momentum or wavenumber
* input from common blocks
*       ompl - zero q energy of mode
*       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp
      dwdq=(q**3+2.d0*adisp*q)/(2.d0*wdisp(q))
      return
      end

      double precision function d2wdq2(q)
* the second derivative of the dispertion relation 
* with respect to the plasmon momentum q
* input: q - momentum or wavenumber
* input from common blocks
*       omp - plasma frequency omega_p
*       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q,wdisp,dwdq
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp,dwdq
      d2wdq2=(3.d0*q**2+2.d0*adisp)*wdisp(q)
     2       -(q**3+2.d0*adisp*q)*dwdq(q)
      d2wdq2=d2wdq2/(2.d0*wdisp(q)**2)
      return
      end

      double precision function qdisp(w)
* The inverse dispertion relation
* input: w - energy (omega)
* input from common blocks
*       ompl - zero q energy of mode
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
      implicit none
      double precision w,x,y
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      x=adisp**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*adisp+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          qdisp=dsqrt(y)
          return
        endif
      endif
      qdisp=0.d0
      return
      end

      double precision function qthresh(ompl,AA,ef,qf)
* Find the photoelectron momentum corresponding to the onset of plasmon
* losses.
* input: ompl - zero q energy of mode
*        AA - dispersion parameter omega(q)**2=ompl**2+AA*q**2+q**4/4
*        ef - fermi energy
*        qf - fermi momentum
      implicit none
      integer i,nrts,nrts2,nflip,nqa,nqb
      double precision AA,ompl,ef,qf,a,b,c,d,qthresh1,qthresh2,q01,
     2                 q02,xfact,test1,test2,test3,qh,ek,q1a,q2a,
     3                 q3a,q0a,q1b,q2b,q3b,q0b,w,x,y
      complex*16 rt1,rt2,rt3,rtt1,rtt2,rtt3,rt(3),rtt(3),qst
      a=1.d0
      b=-3.d0*AA
      c=3.d0*AA**2-27.d0*ompl**2/4.d0
      d=-AA**3
      call croots(a,b,c,d,rt1,rt2,rt3,nrts)
      rt(1)=rt1
      rt(2)=rt2
      rt(3)=rt3
      if (nrts.eq.1) then
 10     continue
          nflip=0
          do i=1,2
            if (dimag(rt(i)).lt.dimag(rt(i+1))) then
              qst=rt(i)
              rt(i)=rt(i+1)
              rt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 10
        qthresh1=dble(rt(2))
      else
        qthresh1=max(dble(rt1),dble(rt2),dble(rt3))
      endif
*      write(6,*) rt1,rt2,rt3
*      write(6,*) qthresh1
      if (qthresh1.gt.0.d0) then
        qthresh1=dsqrt(qthresh1)
      else
        qthresh1=0.d0
      endif
*      write(6,*) qthresh1
      a=1.d0
      b=1.5d0*qf+AA/qf
      c=qf**2+2.d0*AA
      d=qf**3/4.d0+AA*qf+ompl**2/qf
      call croots(a,b,c,d,rt1,rt2,rt3,nrts)
      rt(1)=rt1
      rt(2)=rt2
      rt(3)=rt3
      b=-b
      d=-d
      call croots(a,b,c,d,rtt1,rtt2,rtt3,nrts2)
      rtt(1)=rtt1
      rtt(2)=rtt2
      rtt(3)=rtt3
      if (nrts.eq.1) then
 11     continue
          nflip=0
          do i=1,2
            if (dimag(rt(i)).lt.dimag(rt(i+1))) then
              qst=rt(i)
              rt(i)=rt(i+1)
              rt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 11
        q01=dble(rt(2))
      else
        xfact=sqrt(AA**2+(dble(rt(1))**2/2.d0)**2-ompl**2)-AA
        test1=dble(rt(1))-qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rt(2))**2/2.d0)**2-ompl**2)-AA
        test2=dble(rt(2))-qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rt(3))**2/2.d0)**2-ompl**2)-AA
        test3=dble(rt(3))-qf-sqrt(2.d0*xfact)
        if(test1.lt.test2.and.test1.lt.test3) then
          q01=dble(rt(1))
        elseif(test2.lt.test3) then
          q01=dble(rt(2))
        else
          q01=dble(rt(3))
        endif
      endif
      if (nrts2.eq.1) then
 12     continue
          nflip=0
          do i=1,2
            if (dimag(rtt(i)).lt.dimag(rtt(i+1))) then
              qst=rtt(i)
              rtt(i)=rtt(i+1)
              rtt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 12
        q02=dble(rtt(2))
      else
        xfact=sqrt(AA**2+(dble(rtt(1))**2/2.d0)**2-ompl**2)-AA
        test1=dble(rtt(1))+qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rtt(2))**2/2.d0)**2-ompl**2)-AA
        test2=dble(rtt(2))+qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rtt(3))**2/2.d0)**2-ompl**2)-AA
        test3=dble(rtt(3))+qf-sqrt(2.d0*xfact)
        if(test1.lt.test2.and.test1.lt.test3) then
          q02=dble(rt(1))
        elseif(test2.lt.test3) then
          q02=dble(rt(2))
        else
          q02=dble(rt(3))
        endif
      endif
*      write(6,*) q01,q02
      qthresh2=min(abs(q01),abs(q02))
*      write(6,*) qthresh1,qthresh2

      qh=1000.d0*qf
      ek=qthresh1**2/2.d0
      call qlimits(ek,qthresh1,ompl,AA,qh,nqa,q1a,q2a,q3a)
      q0a=0.d0
      w=ek-ef
      x=AA**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*AA+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          q0a=dsqrt(y)
        endif
      endif
      ek=qthresh2**2/2.d0
      call qlimits(ek,qthresh2,ompl,AA,qh,nqb,q1b,q2b,q3b)
      q0b=0.d0
      w=ek-ef
      x=AA**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*AA+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          q0b=dsqrt(y)
        endif
      endif

      if (nqa.eq.0) then
        qthresh=qthresh1
      elseif(abs(q1a-q2a).lt.abs(q1b-q0b)) then
        qthresh=qthresh1
      else
        qthresh=qthresh2
      endif
*      write(6,*) qthresh
      return
      end

      double precision function vpp2(q)
* the square of the coupling potential
* input: q - momentum or wavenumber
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ompl - zero q energy of mode
      implicit none
      double precision q,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp
      vpp2=2*pi*omp**2/(q**2*wdisp(q))
      return
      end

      subroutine qlimits(w,pk,omp,Ap,qh,nq,q1,q2,q3)
* finds the limiting q values of the inequalities
* omega(q)+(q-k)^2/2-omega < 0 
* omega(q)+(q+k)^2/2-omega > 0 and
* omega(q)-(q-k)^2/2+omega > 0
* for omega(q)^2 = omp^2+Ap q^2+q^4/4
* qh is the upper allowable value for q (usually set by omega(q)+E_f-omega < 0
* for the first two inequalities.
* input: energy w in atomic hartrees.
*        momentum pk in atomic units.
*        resonance frequency omp in hartrees.
*        dispersion parameter Ap.
*        upper limit qh.
* output: number of limiting q values (either 1 or 3).
*         limiting q values q1, q2,and q3
*         q1 and q2 bracket a region of allowed q values given by
*         the first two inequalities listed above,
*         q3 is the upper bound of the third inequality
      implicit none
      integer i,j,nar,nq
      double precision w,pk,omp,Ap,qh,q1,q2,q3,a,b,c,d
      double precision dev1,dev2,dev3,qa1,qa2,qa3,wdisp
      complex*16 aa1,aa2,aa3
      external wdisp
*      write(6,*) 'qlimits'
* find q for which omega(q)+(q-k)^2/2-omega = 0
      a=pk
      b=w+Ap-3.d0*pk**2/2.d0
      c=pk**3-2.d0*w*pk
      d=omp**2-w**2+w*pk**2-pk**4/4.d0
      call croots(a,b,c,d,aa1,aa2,aa3,nar)
* One of these roots is the solution of
* omega(q)-(q-k)^2/2+omega = 0.  This is q3.  The other two roots, if present,
* solve omega(q)+(q-k)^2/2-omega = 0
* It must hold that q>=0.  Fortunately, the symmetry of the inequalities
* means that solving omega(q)+(q+k)^2/2-omega = 0 gives roots which are the 
* negatives of the roots from the cubic we just solved.  It suffices to
* take the absolute value of the roots (if they are real) to give our
* limiting values of q.
      if (nar.eq.3) then
        qa1=dble(aa1)
        qa2=dble(aa2)
        qa3=dble(aa3)
        dev1=dabs(wdisp(qa1)+(qa1-pk)**2/2.d0-w)
        dev2=dabs(wdisp(qa2)+(qa2-pk)**2/2.d0-w)
        dev3=dabs(wdisp(qa3)+(qa3-pk)**2/2.d0-w)
        if (dev1.gt.dev2.and.dev1.gt.dev3) then
          q1=min(dabs(qa2),dabs(qa3))
          q2=max(dabs(qa2),dabs(qa3))
          q3=dabs(qa1)
        elseif(dev2.gt.dev3) then
          q1=min(dabs(qa1),dabs(qa3))
          q2=max(dabs(qa1),dabs(qa3))
          q3=dabs(qa2)
        else
          q1=min(dabs(qa1),dabs(qa2))
          q2=max(dabs(qa1),dabs(qa2))
          q3=dabs(qa3)
        endif
        q1=min(q1,qh)
        q2=min(q2,qh)
        nq=3
      else
* The equation omega(q)+(q-k)^2/2-omega = 0 has no real solutions.  The
* one real root of the cubic is q3. 
        q1=0.d0
        q2=0.d0
        qa1=dabs(dimag(aa1))
        qa2=dabs(dimag(aa2))
        qa3=dabs(dimag(aa3))
        if (qa1.lt.qa2.and.qa1.lt.qa3) then
          q3=dabs(dble(aa1))
        elseif (qa2.lt.qa3) then
          q3=dabs(dble(aa2))
        else
          q3=dabs(dble(aa3))
        endif
        nq=1
      endif
      return
      end
      subroutine renergies(w,rbeta)
* Calculates the real part of the photelectron self energy due to 
* pole ipl in the inverse dielectric function epsilon^{-1}.
* input: w - energy (omega)
* output: rbeta - the real part of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - width of pole in epsilon^(-1)
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,rbeta,beta,rw1beta,grater,wmax,wmin,
     2                 abr,rlr,xsing,error,exchange,qmax
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision rseint1,rseint2,rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,rw1beta,grater,exchange,rseint1,rseint2,rseint3
      integer ijkwrite
      common /morewrite/ ijkwrite
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
*     rlr for plasmon pole should be 10^-7 to eliminate numerical 
*     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        rbeta=grater(rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        rbeta=rbeta+grater(rseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,pk-qf,pk+qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      elseif (pk.lt.qf) then
        rbeta=grater(rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,qf-pk,pk+qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        if (lowq.ne.0) rbeta=rbeta+grater(rseint3,0.d0,
     2             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        rbeta=grater(rseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,0.d0,2.d0*qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      endif
      rbeta=-rbeta*omp**2/(2.d0*pi*pk)
*      rbeta=rbeta+exchange(pk)
      return
      end

      subroutine brsigma(w,rbeta,xibeta)
* Calculates the broadened photelectron self energy due to 
* pole ipl in the inverse dielectric function epsilon^{-1}.
* input: w - energy (omega)
* output: rbeta - the real part of the self energy
*         xibeta - the imaginary part of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - linewidth of pole in epsilon^{-1}
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns,nq,nqq
      parameter (nqq=5)
      double precision w,rbeta,xibeta,grater,qdisp,wmax,wmin,
     2                 abr,rlr,xsing(20),error
      double precision qmax,qlimh,qliml,qh,q0,q1,q2,q3,qsing(nqq)
      double precision sig1r,sig2r,sig3r,sig4r,sig5r,
     2                 sig6r,sig7r,sig8r,sig9r,sig10r,
     3                 sig1i,sig2i,sig3i,sig4i,sig5i,
     4                 sig6i,sig7i,sig8i,sig9i,sig10i
      double precision fqlogr1,fqlogr2,fqlogr3,fqlogr4,
     2                 fqlogi1,fqlogi2,fqlogi3,fqlogi4,
     3                 fqatnr1,fqatnr2,fqatnr3,fqatnr4,
     4                 fqatni1,fqatni2,fqatni3,fqatni4
      double complex cfact,csigma
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer iwrite,jj
      common /flag/ iwrite,jj
      integer lowq
      common /belowqf/ lowq
      external grater,qdisp,
     2         fqlogr1,fqlogr2,fqlogr3,fqlogr4,
     3         fqlogi1,fqlogi2,fqlogi3,fqlogi4,
     4         fqatnr1,fqatnr2,fqatnr3,fqatnr4,
     5         fqatni1,fqatni2,fqatni3,fqatni4
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
*     rlr for plasmon pole should be 10^-7 to eliminate numerical 
*     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      qlimh=pk+qf
      qliml=abs(pk-qf)
      qh=qdisp(max(wp-ef,ompl))
      q0=qdisp(max(ef-wp,ompl))
      call qlimits(wp,pk,ompl,adisp,qmax,nq,q1,q2,q3)
      qsing(1)=q0
      qsing(2)=q1
      qsing(3)=q2
      qsing(4)=q3
      qsing(5)=qh
      sig1r=0.d0
      sig1i=0.d0
      sig2r=0.d0
      sig2i=0.d0
      sig3r=0.d0
      sig3i=0.d0
      sig4r=0.d0
      sig4i=0.d0
      sig5r=0.d0
      sig5i=0.d0
      sig6r=0.d0
      sig6i=0.d0
      sig7r=0.d0
      sig7i=0.d0
      sig8r=0.d0
      sig8i=0.d0
      sig9r=0.d0
      sig9i=0.d0
      sig10r=0.d0
      sig10i=0.d0
      call findsing(qlimh,qmax,nqq,qsing,nsing,xsing)
      sig1r=grater(fqlogr1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig1i=grater(fqlogi1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig2r=grater(fqatnr1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig2i=grater(fqatni1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      call findsing(qliml,qlimh,nqq,qsing,nsing,xsing)
      sig3r=grater(fqlogr2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig3i=grater(fqlogi2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig4r=grater(fqatnr2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      sig4i=grater(fqatni2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      if (lowq.ne.0) then
        sig7r=grater(fqlogr3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        sig7i=grater(fqlogi3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        sig8r=grater(fqatnr3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        sig8i=grater(fqatni3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      endif
      call findsing(0.d0,qliml,nqq,qsing,nsing,xsing)
      if (pk.gt.qf) then
        sig5r=grater(fqlogr1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig5i=grater(fqlogi1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig6r=grater(fqatnr1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig6i=grater(fqatni1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      endif
      if (pk.lt.qf.and.lowq.ne.0) then
        sig9r=grater(fqlogr4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig9i=grater(fqlogi4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig10r=grater(fqatnr4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        sig10i=grater(fqatni4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      endif
      rbeta=(sig1r+sig3r+sig5r+sig7r+sig9r)*omp**2/(4.d0*pi*pk)
     2     +(sig2r+sig4r+sig6r+sig8r+sig10r)*omp**2/(2.d0*pi*pk)
      xibeta=(sig1i+sig3i+sig5i+sig7i+sig9i)*omp**2/(4.d0*pi*pk)
     2     -(sig2i+sig4i+sig6i+sig8i+sig10i)*omp**2/(2.d0*pi*pk)
      cfact=(1,0)-(brd/ompl)*(0,1)
      csigma=((1,0)*rbeta+(0,1)*xibeta)*cfact
      rbeta=dble(csigma)
      xibeta=dimag(csigma)
*      rbeta=rbeta+exchange(pk)
      return
      end

      double precision function fqlogr1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr1=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi1=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=(ef-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr2=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=(ef-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi2=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/((ef-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr3=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/((ef-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi3=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr4=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi4=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqatni1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp-wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatni1=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp-wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatnr1=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp-wq-ef)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatni2=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp-wq-ef)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatnr2=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-ef)/brd)
      fqatni3=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-ef)/brd)
      fqatnr3=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-(pk+q)**2/2.d0)/brd)
      fqatni4=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-(pk+q)**2/2.d0)/brd)
      fqatnr4=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      subroutine dbrsigma(w,drbeta,dibeta)
* Calculates the derivative of the broadened photelectron self energy 
* with respect to frequency w due to 
* pole ipl in the inverse dielectric function epsilon^{-1}.
* input: w - energy (omega)
* output: drbeta - the real part of the self energy derivative
*         dibeta - the imaginary part of the self energy derivative
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - linewidth of pole in epsilon^{-1}
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns,nq,nqq
      parameter (nqq=5)
      double precision w,drbeta,dibeta,grater,qdisp,wmax,wmin,
     2                 abr,rlr,xsing(20),error
      double precision qmax,qlimh,qliml,qh,q0,q1,q2,q3,qsing(nqq)
      double precision dsig1r,dsig2r,dsig3r,dsig4r,dsig5r,
     2                 dsig6r,dsig7r,dsig8r,dsig9r,dsig10r,
     3                 dsig1i,dsig2i,dsig3i,dsig4i,dsig5i,
     4                 dsig6i,dsig7i,dsig8i,dsig9i,dsig10i
      double precision fqlogr1,fqlogr2,fqlogr3,fqlogr4,
     2                 fqlogi1,fqlogi2,fqlogi3,fqlogi4,
     3                 fqatnr1,fqatnr2,fqatnr3,fqatnr4,
     4                 fqatni1,fqatni2,fqatni3,fqatni4
      double complex cfact,csigma
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer iwrite,jj
      common /flag/ iwrite,jj
      integer lowq
      common /belowqf/ lowq
      external grater,qdisp,
     2         dqlogr1,dqlogr2,dqlogr3,dqlogr4,
     3         dqlogi1,dqlogi2,dqlogi3,dqlogi4,
     4         dqatnr1,dqatnr2,dqatnr3,dqatnr4,
     5         dqatni1,dqatni2,dqatni3,dqatni4
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
*     rlr for plasmon pole should be 10^-7 to eliminate numerical 
*     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      qlimh=pk+qf
      qliml=abs(pk-qf)
      qh=qdisp(max(wp-ef,ompl))
      q0=qdisp(max(ef-wp,ompl))
      call qlimits(wp,pk,ompl,adisp,qmax,nq,q1,q2,q3)
      qsing(1)=q0
      qsing(2)=q1
      qsing(3)=q2
      qsing(4)=q3
      qsing(5)=qh
      dsig1r=0.d0
      dsig1i=0.d0
      dsig2r=0.d0
      dsig2i=0.d0
      dsig3r=0.d0
      dsig3i=0.d0
      dsig4r=0.d0
      dsig4i=0.d0
      dsig5r=0.d0
      dsig5i=0.d0
      dsig6r=0.d0
      dsig6i=0.d0
      dsig7r=0.d0
      dsig7i=0.d0
      dsig8r=0.d0
      dsig8i=0.d0
      dsig9r=0.d0
      dsig9i=0.d0
      dsig10r=0.d0
      dsig10i=0.d0
      call findsing(qlimh,qmax,nqq,qsing,nsing,xsing)
      dsig1r=grater(dqlogr1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig1i=grater(dqlogi1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig2r=grater(dqatnr1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig2i=grater(dqatni1,qlimh,qmax,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      call findsing(qliml,qlimh,nqq,qsing,nsing,xsing)
      dsig3r=grater(dqlogr2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig3i=grater(dqlogi2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig4r=grater(dqatnr2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      dsig4i=grater(dqatni2,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      if (lowq.ne.0) then
        dsig7r=grater(dqlogr3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        dsig7i=grater(dqlogi3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        dsig8r=grater(dqatnr3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
        dsig8i=grater(dqatni3,qliml,qlimh,abr,rlr,nsing,xsing,
     2           error,numcal,maxns)
      endif
      call findsing(0.d0,qliml,nqq,qsing,nsing,xsing)
      if (pk.gt.qf) then
        dsig5r=grater(dqlogr1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig5i=grater(dqlogi1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig6r=grater(dqatnr1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig6i=grater(dqatni1,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      endif
      if (pk.lt.qf.and.lowq.ne.0) then
        dsig9r=grater(dqlogr4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig9i=grater(dqlogi4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig10r=grater(dqatnr4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        dsig10i=grater(dqatni4,0.d0,qliml,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      endif
      drbeta=(dsig1r+dsig3r+dsig5r+dsig7r+dsig9r)*omp**2/(2.d0*pi*pk)
     2     +(dsig2r+dsig4r+dsig6r+dsig8r+dsig10r)*omp**2/(2.d0*pi*pk)
      dibeta=(dsig1i+dsig3i+dsig5i+dsig7i+dsig9i)*omp**2/(2.d0*pi*pk)
     2     -(dsig2i+dsig4i+dsig6i+dsig8i+dsig10i)*omp**2/(2.d0*pi*pk)
      cfact=(1,0)-(brd/ompl)*(0,1)
      csigma=((1,0)*drbeta+(0,1)*dibeta)*cfact
      drbeta=dble(csigma)
      dibeta=dimag(csigma)
*      rbeta=rbeta+exchange(pk)
      return
      end

      double precision function dqlogr1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr1=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi1=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogr2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr2=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi2=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogr3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr3=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi3=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end


      double precision function dqlogr4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr4=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi4=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqatni1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni1=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr1(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr1=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni2=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr2(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr2=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni3=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr3(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr3=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni4=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr4(q)
* Integrand for one of the intergals in calculating the 
* self energy in subroutine brsigma.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
*      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr4=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      subroutine findsing(ql,qh,nqq,qsing,nsing,xsing)
* finds which values in the array qsing are between ql and qh, and puts
* them in the array xsing.  Returns number of elements in xsing as nsing
      implicit none
      integer i,j,k,nqq,nsing,it
      double precision ql,qh,qsing(nqq),xsing(20),store
      integer iwrite,jj
      common /flag/ iwrite,jj
      nsing=0
      do i=1,nqq
        if (qsing(i).gt.ql.and.qsing(i).lt.qh) then
          nsing=nsing+1
          xsing(nsing)=qsing(i)
        elseif (qsing(i).lt.ql.and.qsing(i).gt.qh) then
          nsing=nsing+1
          xsing(nsing)=qsing(i)
        endif
      enddo
      if (nsing.le.1) goto 40
      j=2
 20   k=j
 30   continue
      if (xsing(k-1).gt.xsing(k)) then
        store=xsing(k-1)
        xsing(k-1)=xsing(k)
        xsing(k)=store
        k=k-1
        if (k.gt.1) goto 30
      endif
      j=j+1
      if (j.le.nsing) goto 20
 40   return
      end

      double precision function exchange(qk)
* compute the HF exchange potential for the free electron gas
* input: qk - photoelectron momentum
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       qf - Fermi momentum
      implicit none
      integer i,ii,j,jj
      double precision qk
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      if (qk.eq.qf) then
        exchange=-qf/pi
      else
        exchange=-(1.d0/pi)*(qf+
     2  ((qf**2-qk**2)/(2.d0*qk))*log(dabs((qk+qf)/(qk-qf))))
      endif
      return
      end

      subroutine xienergies(w,xibeta)
* Calculates the imaginary part of the photelectron self energy.
* input: w - energy (omega)
* output: xibeta - the imaginary part of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,xibeta,beta,xiw1beta,grater,wmax,wmin,
     2                 abr,rlr,xsing,error
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,xiw1beta,grater
      xibeta=-beta(w)*pi
      return
      end

      subroutine drenergies(w,drbeta)
* Calculates the real part of the derivative of the photelectron 
* self energy with respect to omega.
* input: w - energy (omega)
* output: drbeta - the imaginary part of the derivative of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,drbeta,beta,grater,wmax,wmin,
     2                 abr,rlr,xsing,error,qmax
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision rseint1,rseint2,rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,grater,drseint1,drseint2,drseint3
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)
*     rlr for plasmon pole should be 10^-7 to eliminate numerical 
*     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        drbeta=grater(drseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        drbeta=drbeta+grater(drseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,pk-qf,pk+qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
      elseif (pk.lt.qf) then
        drbeta=grater(drseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,qf-pk,pk+qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        if (lowq.ne.0) drbeta=drbeta+grater(drseint3,0.d0,
     2             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        drbeta=grater(drseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,0.d0,2.d0*qf,abr,rlr,nsing,
     2             xsing,error,numcal,maxns)
      endif
      drbeta=drbeta*omp**2/(2.d0*pi*pk)
      return
      end

      subroutine dienergies(w,dibeta)
* Calculates the imaginary part of the derivative of the photelectron 
* self energy with respect to omega.
* input: w - energy (omega)
* output: dibeta - the imaginary part of the derivative of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=ompl+adisp*q**2+q**4/4
      implicit none
      integer i,npts,nnpts
      double precision w,dibeta
      integer it1,it2,nq
      double precision qh,q1,q2,q3,q0,wh,wq1,wq2,wq3,wq0,qmax
* the q's are limiting values of momenta, the w's are energies
* corresponding to these limiting momenta
      double precision qdisp,A,wdisp,test1,test2
      double precision dq0dw,dq1dw,dq2dw,dq3dw,dqhdw,xfact
* the dqdw's are derivatives of the momentum limits with
* changing energy
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      A=adisp
      q1=0.d0
      q2=0.d0
      q3=0.d0
*      xfact=0.d0
      dibeta=0.d0
* Must find the momenta which limit the 
* final integration (done analytically).
* Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
* Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
* and omega(q)+omega-(k-q)^2/2=0 for q>0.
      qmax=1.d6*qf
      call qlimits(w+ekp,pk,ompl,A,qmax,nq,q1,q2,q3)
* find derivatives of q's
      wq1=wdisp(q1)
      if (w+ekp-ef.gt.ompl) then
        dqhdw=(w+ekp-ef)/(qh*dsqrt(A**2+(w+ekp-ef)**2-ompl**2))
        dq0dw=0.d0
      elseif(ef-w-ekp.gt.ompl) then
        dq0dw=-(ef-w-ekp)/(q0*dsqrt(A**2+(ef-w-ekp)**2-ompl**2))
        dqhdw=0.d0
      else
        dqhdw=0.d0
        dq0dw=0.d0
      endif
      test1=(pk+q1)**2/2.d0-(w+ekp)+wq1
      test2=(pk-q1)**2/2.d0-(w+ekp)+wq1
      if (q1.ge.qh) then
        q1=qh
        dq1dw=dqhdw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq1dw=wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)
      else
        dq1dw=wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)
      endif
      wq2=wdisp(q2)
      test1=(pk+q2)**2/2.d0-(w+ekp)+wq2
      test2=(pk-q2)**2/2.d0-(w+ekp)+wq2
      if (q2.ge.qh) then
        q2=qh
        dq2dw=dqhdw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq2dw=wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)
      else
        dq2dw=wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)
      endif
      wq3=wdisp(q3)
      test1=(pk+q3)**2/2.d0-(w+ekp)-wq3
      test2=(pk-q3)**2/2.d0-(w+ekp)-wq3
      if (q3.ge.q0) then
        q3=q0
        dq3dw=dq0dw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq3dw=wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)
      else
        dq3dw=wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)
      endif
* find derivative of imaginary part of self energy
* find contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        xfact=A*q1*(1.d0/wq1+1.d0/ompl)+q1**3/(2.d0*wq1)
        xfact=2.d0/q1-xfact/(ompl+wq1+A*q1**2/(2.d0*ompl))
        dibeta=dibeta+omp**2/(4.d0*pk*ompl)*dq1dw*xfact
        xfact=A*q2*(1.d0/wq2+1.d0/ompl)+q2**3/(2.d0*wq2)
        xfact=2.d0/q2-xfact/(ompl+wq2+A*q2**2/(2.d0*ompl))
        dibeta=dibeta-omp**2/(4.d0*pk*ompl)*dq2dw*xfact
      endif
* find contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        xfact=A*q0*(1.d0/wq0+1.d0/ompl)+q0**3/(2.d0*wq0)
        xfact=2.d0/q0-xfact/(ompl+wq0+A*q0**2/(2.d0*ompl))
        dibeta=dibeta+omp**2/(4.d0*pk*ompl)*dq0dw*xfact
        xfact=A*q3*(1.d0/wq3+1.d0/ompl)+q3**3/(2.d0*wq3)
        xfact=2.d0/q3-xfact/(ompl+wq3+A*q3**2/(2.d0*ompl))
        dibeta=dibeta-omp**2/(4.d0*pk*ompl)*dq3dw*xfact
      endif
      dibeta=dibeta
* Test that dibeta is a number (stop for nan results).
      it1=int(dibeta)
      it2=int(2.d0*dibeta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'dienergies ',dibeta
        stop
      endif
      dibeta=dibeta
      return
      end

      subroutine d2renergies(w,d2rbeta)
* Calculates the real part of the second derivative of the photelectron 
* self energy with respect to omega.
* input: w - energy (omega)
* output: d2rbeta - the real part of the second derivative of the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,d2rbeta,beta,grater,wmax,wmin,
     2                 abr,rlr,xsing,error,qmax
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision d2rseint1,d2rseint2,d2rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,grater,d2rseint1,d2rseint2,d2rseint3
      wp=w+ekp
      wmin=-100*(pi*beta(0.d0)+ompl)+w
      wmax=100*(pi*beta(0.d0)+ompl)+w
      qmax=1.d2*dsqrt(ompl)
*     rlr for plasmon pole should be 10^-7 to eliminate numerical 
*     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        d2rbeta=grater(d2rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,pk-qf,pk+qf,abr,rlr,nsing,
     2             xsing,error,numcal,maxns)
      elseif (pk.lt.qf) then
        d2rbeta=grater(d2rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,qf-pk,pk+qf,abr,rlr,nsing,
     2             xsing,error,numcal,maxns)
        if (lowq.ne.0) d2rbeta=d2rbeta+grater(d2rseint3,0.d0,
     2             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        d2rbeta=grater(d2rseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,
     2             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,0.d0,2.d0*qf,abr,rlr,nsing,
     2             xsing,error,numcal,maxns)
      endif
      d2rbeta=d2rbeta*omp**2/(2.d0*pi*pk)
      return
      end

      subroutine d2ienergies(w,d2ibeta)
* Calculates the imaginary part of the second derivative of the photelectron 
* self energy with respect to omega.
* input: w - energy (omega)
* output: d2ibeta - the imaginary part of the second derivative of 
*                   the self energy
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      integer i,npts,nnpts
      double precision w,d2ibeta
      integer it1,it2,nq
      double precision qh,q1,q2,q3,q0,wh,wq1,wq2,wq3,wq0,qmax
* the q's are limiting values of momenta, the w's are energies
* corresponding to these limiting momenta
      double precision qdisp,A,wdisp,test1,test2,www,d2fact,xfact
      double precision dq0dw,dq1dw,dq2dw,dq3dw,dqhdw,dwqdw
* the dqdw's are derivatives of the momentum limits with
* changing energy
      double precision d2q0dw2,d2q1dw2,d2q2dw2,d2q3dw2,d2qhdw2
* the d2qdw2's are second derivatives of the momentum limits with
* changing energy (but you've already figured that out by now, right?)
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp,d2fact
      A=adisp
      q1=0.d0
      q2=0.d0
      q3=0.d0
      d2ibeta=0.d0
* Must find the momenta which limit the 
* final integration (done analytically).
* Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
* Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
* and omega(q)+omega-(k-q)^2/2=0 for q>0.
      qmax=1.d6*qf
      call qlimits(w+ekp,pk,ompl,A,qmax,nq,q1,q2,q3)
* find derivatives of q's
      wq1=wdisp(q1)
      if (w+ekp-ef.gt.ompl) then
        www=A**2+(w+ekp-ef)**2-ompl**2
        dqhdw=(w+ekp-ef)/(qh*dsqrt(www))
        dq0dw=0.d0
        d2qhdw2=1/(qh*dsqrt(www))-(w+ekp-ef)**2/(qh*dsqrt(www)**3)
     2          -(w+ekp-ef)*dqhdw/(qh**2*dsqrt(www))
        d2q0dw2=0.d0
      elseif(ef-w-ekp.gt.ompl) then
        www=A**2+(w+ekp-ef)**2-ompl**2
        dq0dw=(w+ekp-ef)/(q0*dsqrt(www))
        dqhdw=0.d0
        d2q0dw2=1/(q0*dsqrt(www))-(w+ekp-ef)**2/(q0*dsqrt(www)**3)
     2          -(w+ekp-ef)*dq0dw/(q0**2*dsqrt(www))
        d2qhdw2=0.d0
      else
        dqhdw=0.d0
        dq0dw=0.d0
        d2qhdw2=0.d0
        d2q0dw2=0.d0
      endif
      test1=(pk+q1)**2/2.d0-(w+ekp)+wq1
      test2=(pk-q1)**2/2.d0-(w+ekp)+wq1
      if (q1.ge.qh) then
        q1=qh
        dq1dw=dqhdw
        d2q1dw2=d2qhdw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq1dw=wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)
        dwqdw=(A*q1+q1**3/2.d0)/wq1*dq1dw
        d2q1dw2=dwqdw/((q1+pk)*wq1+A*q1+q1**3/2.d0)
     2  -wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)
     3  *((wq1+A+3*q1**2/2.d0)*dq1dw+(q1+pk)*dwqdw)
      else
        dq1dw=wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)
        dwqdw=(A*q1+q1**3/2.d0)/wq1*dq1dw
        d2q1dw2=dwqdw/((q1-pk)*wq1+A*q1+q1**3/2.d0)
     2  -wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)
     3  *((wq1+A+3*q1**2/2.d0)*dq1dw+(q1-pk)*dwqdw)
      endif
      wq2=wdisp(q2)
      test1=(pk+q2)**2/2.d0-(w+ekp)+wq2
      test2=(pk-q2)**2/2.d0-(w+ekp)+wq2
      if (q2.ge.qh) then
        q2=qh
        dq2dw=dqhdw
        d2q2dw2=d2qhdw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq2dw=wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)
        dwqdw=(A*q2+q2**3/2.d0)/wq2*dq2dw
        d2q2dw2=dwqdw/((q2+pk)*wq2+A*q2+q2**3/2.d0)
     2  -wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)
     3  *((wq2+A+3*q2**2/2.d0)*dq2dw+(q2+pk)*dwqdw)
      else
        dq2dw=wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)
        dwqdw=(A*q2+q2**3/2.d0)/wq2*dq2dw
        d2q2dw2=dwqdw/((q2-pk)*wq2+A*q2+q2**3/2.d0)
     2  -wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)
     3  *((wq2+A+3*q2**2/2.d0)*dq2dw+(q2-pk)*dwqdw)
      endif
      wq3=wdisp(q3)
      test1=(pk+q3)**2/2.d0-(w+ekp)-wq3
      test2=(pk-q3)**2/2.d0-(w+ekp)-wq3
      if (q3.ge.q0) then
        q3=q0
        dq3dw=dq0dw
        d2q3dw2=d2q0dw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq3dw=wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)
        dwqdw=(A*q3+q3**3/2.d0)/wq3*dq3dw
        d2q3dw2=dwqdw/((q3+pk)*wq3-A*q3-q3**3/2.d0)
     2  -wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)
     3  *((wq3-A-3*q3**2/2.d0)*dq3dw+(q3+pk)*dwqdw)
      else
        dq3dw=wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)
        dwqdw=(A*q3+q3**3/2.d0)/wq3*dq3dw
        d2q3dw2=dwqdw/((q3-pk)*wq3-A*q3-q3**3/2.d0)
     2  -wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)
     3  *((wq3-A-3*q3**2/2.d0)*dq3dw+(q3-pk)*dwqdw)
      endif
* find second derivative of imaginary part of self energy
* find contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        d2ibeta=d2ibeta+d2fact(q1,wq1,dq1dw,d2q1dw2)
        d2ibeta=d2ibeta-d2fact(q2,wq2,dq2dw,d2q2dw2)
      endif
* find contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        d2ibeta=d2ibeta+d2fact(q0,wq0,dq0dw,d2q0dw2)
        d2ibeta=d2ibeta-d2fact(q3,wq3,dq3dw,d2q3dw2)
      endif
* Test that d2ibeta is a number (stop for nan results).
      it1=int(d2ibeta)
      it2=int(2.d0*d2ibeta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2ienergies ',d2ibeta
        stop
      endif
      d2ibeta=d2ibeta
      return
      end

      double precision function d2fact(q,wq,dqdw,d2qdw2)
* the terms in a sum used to calculate the imaginary part of
* the second derivative of the self energy in the subroutine
* d2ienergies.
* input: q - momentum limit
*        wq - energy corresponding to q
*        dqdw - derivative of momentum limit with respect to energy
*        d2qdw2 - second derviative of momentum limit
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      double precision q,wq,dqdw,d2qdw2,xfact1,xfact2,dwqdw
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,A
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,A
      dwqdw=(A*q+q**3/2.d0)/wq*dqdw
      xfact1=((A*(1.d0/wq+1.d0/ompl)+3.d0*q**2/(2.d0*wq))*dqdw
     2  -dwqdw*(q**3+2.d0*A*q)/(2*wq**2))
     3  /(ompl+wq+A*q**2/(2.d0*ompl))
      xfact1=-xfact1+(dwqdw+A*q*dqdw/ompl)
     2  *(A*q*(1.d0/wq+1.d0/ompl)+q**3/(2.d0*wq))
     3  /(ompl+wq+A*q**2/(2.d0*ompl))**2
      xfact1=xfact1-2.d0*dqdw/q**2
      xfact2=A*q*(1.d0/wq+1.d0/ompl)+q**3/(2.d0*wq)
      xfact2=2.d0/q-xfact2/(ompl+wq+A*q**2/(2.d0*ompl))
      d2fact=omp**2/(4.d0*pk*ompl)*(dqdw*xfact1+d2qdw2*xfact2)
      return
      end

      double precision function rseint1(q)
* Integrand for one of the intergals in calculating the real part
* of the self energy in subroutine renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      integer ijkwrite
      common /morewrite/ ijkwrite
      wq=wdisp(q)
      xlog=((pk+q)**2/2.d0-wp+wq)**2+(acc*ompl)**2
      xlog=xlog/(((pk-q)**2/2.d0-wp+wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint1=xlog/(wq*dsqrt(q**2+ompl*acc))
** NAN detector
*      it1=int(rseint1)
*      it2=int(2.d0*rseint1)
*      if (it1.eq.it2.and.it1.gt.5) then
*        write(6,*) 'rseint1 ',rseint1,q,wq,xlog,wp,pk,acc,ompl
*        stop
*      endif
      return
      end

      double precision function rseint2(q)
* Integrand for one of the intergals in calculating the real part
* of the self energy in subroutine renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xlog=1.d0
      if (lowq.ne.0) then
        xlog=(ef-wp-wq)**2+(acc*ompl)**2
        xlog=xlog/(((pk-q)**2/2.d0-wp-wq)**2+(acc*ompl)**2)
      endif
      xlog=xlog*(((pk+q)**2/2.d0-wp+wq)**2+(acc*ompl)**2)
      xlog=xlog/((ef-wp+wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint2=xlog/(wq*dsqrt(q**2+ompl*acc))
** NAN detector
*      it1=int(rseint2)
*      it2=int(2.d0*rseint2)
*      if (it1.eq.it2.and.it1.gt.5) then
*        write(6,*) 'rseint2 ',rseint2,q,wq,xlog,wp,pk,acc,ompl
*        stop
*      endif
      return
      end

      double precision function rseint3(q)
* Integrand for one of the intergals in calculating the real part
* of the self energy in subroutine renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xlog=((pk+q)**2/2.d0-wp-wq)**2+(acc*ompl)**2
      xlog=xlog/(((pk-q)**2/2.d0-wp-wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint3=xlog/(wq*dsqrt(q**2+ompl*acc))
** NAN detector
*      it1=int(rseint3)
*      it2=int(2.d0*rseint3)
*      if (it1.eq.it2.and.it1.gt.5) then
*        write(6,*) 'rseint3 ',rseint3,q,wq,xlog,wp,pk,acc,ompl
*        stop
*      endif
      return
      end

      double precision function drseint1(q)
* Integrand for one of the intergals in calculating the derivative 
* of the real part of the self energy in subroutine drenergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp+wq
      xfact=xnum1/(xnum1**2+brd**2)
      xnum2=(pk-q)**2/2.d0-wp+wq
      xfact=xfact-xnum2/(xnum2**2+brd**2)
      drseint1=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint1)
      it2=int(2.d0*drseint1)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint1 ',drseint1,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function drseint2(q)
* Integrand for one of the intergals in calculating the derivative 
* of the real part of the self energy in subroutine drenergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2,xnum3,xnum4
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xfact=0.d0
      if (lowq.ne.0) then
        xnum1=ef-wp-wq
        xfact=xnum1/(xnum1**2+brd**2)
        xnum2=(pk-q)**2/2.d0-wp-wq
        xfact=xfact-xnum2/(xnum2**2+brd**2)
      endif
      xnum3=(pk+q)**2/2.d0-wp+wq
      xfact=xfact+xnum3/(xnum3**2+brd**2)
      xnum4=ef-wp+wq
      xfact=xfact-xnum4/(xnum4**2+brd**2)
      drseint2=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint2)
      it2=int(2.d0*drseint2)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint2 ',drseint2,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function drseint3(q)
* Integrand for one of the intergals in calculating the derivative 
* of the real part of the self energy in subroutine drenergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp-wq
      xfact=xnum1/(xnum1**2+brd**2)
      xnum2=(pk-q)**2/2.d0-wp-wq
      xfact=xfact-xnum2/sqrt(xnum2**2+brd**2)
      drseint3=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint3)
      it2=int(2.d0*drseint3)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint3 ',drseint3,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function d2rseint1(q)
* Integrand for one of the intergals in calculating the second derivative 
* of the real part of the self energy in subroutine d2renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp+wq
      xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
      xnum2=(pk-q)**2/2.d0-wp+wq
      xfact=xfact-(xnum2**2-brd**2)/(xnum2**2+brd**2)**2
      d2rseint1=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint1)
      it2=int(2.d0*d2rseint1)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint1 ',d2rseint1,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function d2rseint2(q)
* Integrand for one of the intergals in calculating the second derivative 
* of the real part of the self energy in subroutine d2renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2,xnum3,xnum4
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xfact=0.d0
      if (lowq.ne.0) then
        xnum1=ef-wp-wq
        xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
        xnum2=(pk-q)**2/2.d0-wp-wq
        xfact=xfact-(xnum2**2-brd**2)/(xnum2**2+brd**2)**2
      endif
      xnum3=(pk+q)**2/2.d0-wp+wq
      xfact=xfact+(xnum3**2-brd**2)/(xnum3**2+brd**2)**2
      xnum4=ef-wp+wq
      xfact=xfact-(xnum4**2-brd**2)/(xnum4**2+brd**2)**2
      d2rseint2=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint2)
      it2=int(2.d0*d2rseint2)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint2 ',d2rseint2
        stop
      endif
      return
      end

      double precision function d2rseint3(q)
* Integrand for one of the intergals in calculating the second derivative 
* of the real part of the self energy in subroutine d2renergies.
* input: q - momentum (variable to be integrated over)
* input from common blocks
*       ompl - energy of pole ipl in epsilon^{-1}
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
* common block control of subprograms
*       ac2 - additional accuracy parameter
*       wp - omega prime, an additional energy variable to be held
*            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp-wq
      xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
      xnum2=(pk-q)**2/2.d0-wp-wq
      xfact=xfact-(xnum2**2-brd**2)/sqrt(xnum2**2+brd**2)**2
      d2rseint3=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint3)
      it2=int(2.d0*d2rseint3)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint3 ',d2rseint3
        stop
      endif
      return
      end

      double precision function beta(w)
* the extrinsic beta function
* beta(k,w)=(1/pi)*|Im(self energy(k,w+(k^2)/2))|*theta((w+(k^2)/2)-xmu)
* input: w - energy (omega)
* input from common blocks
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       xmu - chemical potential = Fermi energy + self consistent 
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole ipl in epsilon^{-1}
*       wt - weight of pole ipl in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of 
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       pk - photoelectron momentum
*       acc - global accuracy parameter 
*       brd - global broadening parameter to stabilize logarithms
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      integer it1,it2,nq,i
      double precision w,q1,q2,q3,wth,wq0,wq1,wq2,wq3
      double precision q0,qh,qdisp,A,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
*      common /temp/ q1,q2
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      A=adisp
      q1=0.d0
      q2=0.d0
* Must find the momenta which limit the
* final integration (done analytically).
      beta=0.d0
* Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
* Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
* and omega(q)+omega-(k-q)^2/2=0.
      call qlimits(w+ekp,pk,ompl,A,qh,nq,q1,q2,q3)
* Calculate beta
* Calculate contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        beta=beta+omp**2/(4*pi*pk*ompl)
     2       *dlog(q2**2/(ompl+wq2+A*q2**2/(2.d0*ompl))
     2       *(ompl+wq1+A*q1**2/(2.d0*ompl))/q1**2)
*        beta=beta+ompl/(4*pi*pk)
*     2       *dlog(q2**2/(ompl+wq2+A*q2**2/(2.d0*ompl)))
      endif
* Calculate contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        beta=beta-omp**2/(4*pi*pk*ompl)
     2       *dlog(q0**2/(ompl+wq0+A*q0**2/(2.d0*ompl))
     2       *(ompl+wq3+A*q3**2/(2.d0*ompl))/q3**2)
      endif
* Test that beta is a number (stop for nan results).
      it1=int(beta)
      it2=int(2.d0*beta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'beta ',beta,q1,q2,q3,q0
        write(6,*) 'beta ',pk,w,ekp,acc,ompl
        write(6,*) nq
        write(6,*) q1,q2
        stop
      endif
      return
      end
      
* Program so2conv, by Luke Campbell (2002)

      subroutine so2conv (ispec, ipr6, ipse, ipsk, wsigk, cen, cfname)

* This program computes many body effects on an XAS spectrum given in
* one or more FEFF output files of the form chi.dat,chipnnnn.dat,
* xmu.dat, or feffnnnn.dat.  The files read are from the most recent FEFF 
* from the directory in which this (sub) program is run, as determined
* by the files feff.inp and list.dat.  The files are overwritten
* with the computed spectra corrected for many body effects.

      implicit none
      integer npts,npts2,j,nsfpts,nqpts,maxfiles,ifiles,ifeffct
*      parameter(npts=400,npts2=401,nsfpts=80,nqpts=66,maxfiles=200)
      parameter(npts=400,npts2=401,nsfpts=112,nqpts=66,maxfiles=200)
*       npts - The number of points of the uniform energy grid used for 
*              the spectral function.
*       npts2 - The maximum number of data points in the data files
*       j - Used to store the number of data points in a data file.
*       nsfpts - Number of points of the minimal energy grid used for 
*              the spectral function.
*       nqpts - Number of points in the minimal momentum grid used for
*              the spectral function.
*       maxfiles - The maximum expected number of FEFF files to be read.
*       ifiles - The actual number of FEFF files specified.
*       ifeffct - Number of data points in a feffnnnn.dat file.
      integer i,ii,jj,ik,ikk,ipl
*       i,ii,jj,ik,ikk,ipl - counters for loops.
      integer last,lfirst,isearch,ifirst(maxfiles),ilast(maxfiles)
*       last - Collum number of last non-blank character in text string.
*       lfirst - Collum number of first non-blank character in text string.
*       isearch - Another counter for loops, used in text string searches.
*       ifirst - Collum number of start of file name.
*       ilast - Collum number of end of file name.
      integer ispec,ipr6
*       ispec - Identifies type of spectroscopy.  0 for EXAFS, 1 for XANES.
*       ipr6 - Identifes which files need to be convoluted.
      integer iwrite,iwrite2,intout,iwrconv,iwrcfile,iout,ipw,
     2        ipwc,ipse,ipsk
      double precision cen,den,den2
      character*12 cfname
*       iwrite - Write the spectral function and self energy 
*              to a file at this point on the minimal grid.
*       iwrite2 - Write the spectral function to a file at this point 
*              on the uniform grid.
*       intout - Set to one to write out the running integration of 
*              the convolution with the spectral function.
*       iwrconv - The data point at which the running integration is
*              to be written.
*       iwrcfile - The file for which the running integration is to be 
*              written.
*       ipw - Set greater than zero to write a file of integrated spectral
*              weights.
*       ipwc - Set greater than zero to write a file of integrated spectra
*              weights with the extrinsic spectral function separated into
*              quasiparticle and satellite terms.
*       ipse - Set greater to zero to write a file of the on shell
*              self energy.
*       ipsk - Set greater than zero to write the spectral function 
*              and self energy at a specified momentum.
*       iout - Flags whether the file specfunct.dat needs to
*              be created.
*       cen -  closest energy point to write out running convolution
*       den -  difference in gridpoint energy and cen (above)
*       cfname - name of file shown in running convolution (NULL if none)
      double precision wsigk
*       wsigk - Momentum at which to write self energy file if 
*              ipsk specified.
      integer iedge,ios,itype(maxfiles),iasym,iasymt,isattype,isattypt,
     2        icut,nleg,npi
*       iedge - Pointer to the Fermi energy in xmu.dat.
*       ios - Greater than 1 upon failure to open a file.
*       itype - Identifies type of file read.
*       iasym - Set to 1 to include quasiparticle phase as an
*              asymmetric 1/omega term to the extrinsic satellite
*              rather than as a complex spectral weight.  This is
*              necessary when convoluting with a real valued
*              function or one whose imaginary part is not known as
*              is the case with the xmu.dat file.
*       iasymt - Check if spectral function read from file uses
*              the same method of calculating the quasiparticle phase.
*       isattype - Indicates what approximations to use for 
*              satellite.
*       isattypt - Check if spectral function read from file uses
*              the same approximation.
*       icut - Set to 0 to avoid cutting off the spectral function at
*              energies where excitations are energetically forbidden
*              in the convolution.
*       nleg - Number of legs in the scattering path, not used for 
*              any calculations.
*       npi - Number of 2 pi jumps in phase.
      integer ipath,nlegs
*       Data from file list.dat.  Only ipath is used.
*       ipath - Path index.
      double precision sig2,ampratio,reff
*       Data from file list.dat.  Not used in this program except 
*       as dummy variables.
      character*80 cfirst,cinp,cblankl
      character*12 cfile(maxfiles)
      character*4 cpath
*       cfirst - Lines from FEFF files.
*       cinp - Lines from so2conv.inp.
*       cblankl - A long line of blanks.
*       cfile - The names of the FEFF files to read.
*       cpath - text containing path index
      character*20 cblanks
*       cblanks - A short line of blanks.
      character*10 vintx,mutx
*       Text containing the interstitial potential and the chemical
*       potential, respectively.
      character*9 gchtx,kftx
*       Text containing the core hole lifetime and the Fermi level,
*       respectively.
      character*5 rstx
*       Text containing the interstitial electron density parameter R_s.
      double precision aangstrom,eV
*       Conversion constants.
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      double precision dw,w,wp,dk,dk2,dp,delta
*       Various energy and momenta and their intervals.
      double precision rs,vint,deg,Rnn,gammach,conc,brpl1,brpl,cmu,
     2                 ckf,rst,gcht,brplt
*       rs - Electron density parameter R_s.
*       vint - Interstitial potential.
*       deg - Path degeneracy, number of identical scattering paths.
*       Rnn - Half length of scattering path.
*       gammach - Core hole lifetime.
*       conc - Electron concentration.
*       brpl1 - plasmon broadening in units of plasma frequency
*       brpl - plasmon broadening
*       cmu - Chemical potential from FEFF file.
*       ckf - Fermi level from FEFF file.
*       rst - Test to see if rs matches old value.
*       gcht - Test to see if gammach matches old value.
*       brplt - Test to see if brpl1 matches old value.
      double precision sef0,se0,se,ce,width,z1,z1i,se2,xise,
     2                 zkk,dsig,xpkg(npts2),seg(npts2),ekpg(npts2),
     3                 sse,sxise
*       sef0 - Self energy at the Fermi level Sigma(E_F,k_F)
*       se0 - Non-self consistent on shell self energy.
*       se - Real part of on shell self energy.
*       ce - Electron gas core hole self energy.
*       width - Quasiparticle lifetime.
*       z1 - Real part of renormalization constant.
*       z1i - Imaginary part of renormalization constant.
*       se2 - Real part of self energy (not neccessarily on shell).
*       xise - Imaginary part of self energy (not neccessarily on shell).
*       zkk - Momentum derivative renormalization constant.
*       dsig - Difference in self energies, to estimate zkk.
*       xpkg - Array of momenta, to estimate zkk.
*       seg - Array of self energies, to estimate zkk.
*       ekpg - Array of energies, to estimate zkk.
*       sse - self energy due to given pole in epsilon^{-1}
      common /energies/ se,ce,width,z1,z1i,se2,xise
      double precision xktest,xktold
*       Used for checking the location of the fermi level in xmu.dat.
      double precision xreduc,xfact1
*       xreduc - Ad hoc overall uniform amplitude reduction.  Obsolete.
*       xfact1 - Dummy variable for storing intermediate calculations.
      double precision pk(npts2),pgrid(nqpts),epts(nsfpts),wpts(npts),
     2                 pthresh
*       pk - Array of momentum values at data points from FEFF files.
*       pgrid - Array of momentum values on minimal momentum grid.
*       epts - Array of energy values on minimal energy grid.
*       wpts - Array of energy variables on uniform energy grid.
*       pthresh - Estimate of momentum threshold for plasmon creation.
      double precision spectf(8,nsfpts),cspec(npts),weights(8)
*       spectf - Spectral function computed on minimal grid.
*       cspec - Spectral function interpolated onto uniform grid.
*       Weights - Spectral weights of the components of the spectral function.
      double precision dphsum,s02sum,xnorm,ww
*       Used in finite element integration for averaging over
*            data points on uniform momentum grid to get output 
*            for momentum points in feffnnnn.dat file.
*       dphsum - Phase shift integration.
*       s02sum - Amplitude reduction integration.
*       xnorm - Normalization.
*       ww - Weighting or importance function.
      double precision chirr,chiii,so2mag,phaseshft,
     2                 phshftold,phrmu,phrmu0,
     3                 xchir,xchii,phchir,phchii,xmu2,xmu02,
     4                 phmu,phmu0,xchi2,phchi,rmu2,rmu02,ximu2
*       chirr - Real part of chi.
*       chiii - Imaginary part of chi.
*       so2mag - Amplitude reduction.
*       phaseshift - Phaseshift.
*       phaseshiftold - Used to check for jumps of 2 pi in phase.
*       xchir - Intermediate real part of chi.
*       xchii - Intermediate imaginary part of chi.
*       phchir - Phase of convolution of real part of chi.
*       phchii - Phase of convolution of imaginary part of chi.
*       xmu2 - Convolution of XANES spectrum with spectral function.
*       xmu02 - Convolution of background with spectral function.
*       xchi2 - Convolution of XAFS signal from xmu.dat.
*       phmu - 'Phase' of XANES spectrum, should be zero, but the 
*           subroutine needs a variable to put it in.
*       phmu0 - 'Phase' of background, should be zero, but the 
*           subroutine needs a variable to put it in.
*       phchi - 'Phase' of XAFS signal, should be zero, but the 
*           subroutine needs a variable to put it in.
      double precision emsf(nqpts,nsfpts),essf(nqpts,nsfpts),
     2                 xmsf(nqpts,nsfpts),xssf(nqpts,nsfpts),
     3                 xissf(nqpts,nsfpts),escsf(nqpts,nsfpts),
     4                 engrid(nqpts,nsfpts),wgts(nqpts,8),
     5                 sfinfo(nqpts,8)
*       emsf - Array of extrinsic quasiparticle spectral function
*            values for file specfunct.dat.
*       essf - Array of extrinsic satellite spectral function
*            values for file specfunct.dat.
*       xmsf - Array of interference quasiparticle spectral function
*            values for file specfunct.dat.
*       xssf - Array of interference satellite spectral function
*            values for file specfunct.dat.
*       xissf - Array of intrinsic satellite spectral function
*            values for file specfunct.dat.
*       essf - Array of clipped extrinsic satellite spectral function
*            values for file specfunct.dat.
*       engrid - Energy grid for spectral function in file specfunct.dat.
*       wgts - Array of spectral weights in file specfunct.dat.
*       sfinfo - Array of miscellaneous information in file specfunct.dat.
      double precision xk(npts2),chi(npts2),xmag(npts2),phase(npts2),
     2                 phm2kr(npts2),epts2(npts2),chir(npts2),
     3                 chii(npts2),e1(npts2),xmu(npts2),xmu0(npts2),
     5                 xk2(npts2),caph(npts2),xmfeff(npts2),
     6                 phfeff(npts2),redfac(npts2),xlam(npts2),
     7                 realck2(npts2),caph2(npts2),xmfeff2(npts2),
     8                 phfeff2(npts2),redfac2(npts2),xlam2(npts2)
      double precision rmu(npts2),rmu0(npts2),ximu(npts2)
*       xk - EXAFS wavenumber from chi.dat or chipnnnn.dat.
*       chi - EXAFS signal from chi.dat or chipnnnn.dat.
*       xmag - Magnitude of chi from chi.dat or chipnnnn.dat.
*       phase - Phase of chi from chi.dat or chipnnnn.dat.
*       phm2kr - PHase Minus 2 K R, phase of chi with dominant 
*           k dependant oscillation removed from chi.dat or chipnnnn.dat.
*       epts2 - Energy from files, measured from edge.
*       chir - Real part of EXAFS signal.
*       chii - Imaginary part of EXAFS signal.
*       e1 - Energy from xmu.dat.
*       xmu - XANES signal from xmu.dat.
*       xmu0 - Embedded atom background signal.
*       xk2 - EXAFS wavenumber from feffnnnn.dat
*       caph - Central atom phase shift interpolated onto uniform 
*           momentum grid.
*       xmfeff - Magnitude of f_(eff) interpolated onto uniform
*           momentum grid.
*       phfeff - Phase of f_(eff) interpolated onto uniform 
*           momentum grid.
*       redfac - Reduction factor on uniform momentum grid.
*       xlam - Mean free path of photoelectron, interpolated
*           onto uniform momentum grid.
*       realck2 - Real part of complex momentum from feffnnnn.dat.
*       caph2 - Central atom phase shift on feffnnnn.dat grid.
*       xmfeff2 - Magnitude of f_(eff) from feffnnnn.dat.
*       phfeff2 - Phase of f_(eff) from feffnnnn.dat.
*       redfac2 - Reduction factor feffnnnn.dat grid
*       xlam2 - Mean free path of photoelectron from feffnnnn.dat.
      double precision s02list(npts2),phlist(npts2)
*       s02list - Array of amplitude reduction values on uniform 
*           momentum grid.
*       phlist - Array of phase shifts on uniform momentum grid.
      double precision emain,esat,xmain,xsat,xisat,sat,eclip,esatr
*       emain - Extrinsic quasiparticle spectral function.
*       esat - Extrinsic satellite spectral function.
*       xmain - Interference quasiparticle spectral function.
*       xsat - Interference satellite spectral function.
*       xisat - Intrinsic spectral function.
*       sat - Total satellite spectral function.
*       eclip - Anomalous extrinsic satellite in region of quasiparticle.
*       esatr - Everything left in the extrinsic satellite after eclip is
*         subtracted off.
      integer npl,nplmax,nplt
      parameter (nplmax=5000)
*       npl - number of poles in epsilon^{-1}
*       nplmax - maximum number of poles in epsilon^{-1} for array dimensioning
      double precision plengy(nplmax),plwt(nplmax),oscstr(nplmax),
     2                 plbrd(nplmax),epswt,
     3                 plengyt(nplmax),plwtt(nplmax),plbrdt(nplmax)
*       plengy - energy of poles in epsilon^{-1}
*       plwt - weight of poles in epsilon^{-1}, such that sum(plwt)=1
*       oscstr - oscilator strength f of resonances epsilon^{-1}
*       plbrd - broadening of poles in epsilon^{-1}
*       epswt - sum of weight of poles, to enforce sum(plwt)=1 {a consequence of
*               int(0,oo) d_omega omega epsilon^{-1}(q,omega)=-(pi/2)*omp}
      logical brpole
*       brpole - true if pole broadening is to be calculated
      double precision pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
*       pi - ratio of circumference to diameter of a circle in
*            euclidian geometry
*       ef - Fermi energy
*       fmu - chemical potential = Fermi energy + self consistent
*             on shell self energy at the Fermi level
*       qf - Fermi momentum
*       omp - plasma frequency omega_p
*       ompl - energy of pole in epsilon^{-1}
*       wt - weight of pole in epsilon^{-1}
*       ekp - photoelectron energy = bare kinetic energy + real part of
*             on shell self energy
*       ek - bare photoelectron kinetic energy = pk**2/2
*       qpk - photoelectron momentum
*       acc - global accuracy parameter
*       brd - width of pole in epsilon^(-1)
*       adisp - dispersion parameter for dispersion relation,
*               w(q)**2=omp+adisp*q**2+q**4/4
      common /convsf/ pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
*       Common block one is used in most functions and subroutines.
      common /flag/ iwrite,jj
*       Common block flag tells subroutine  mkspectf when to print out
*       interesting information.
      double precision qthresh,beta,exchange
*       Functions which may be called.
      integer lowq,lowqt
*       lowq - Set equal to 0 to avoid calculating self energy from 
*              contributions below Fermi level.  When considering 
*              interference effects, lowq should be set equal to zero.
*              When finding quasiparticle properties, lowq should not
*              be zero.
*       lowqt - Value of lowq from previous run.
      common /belowqf/ lowq
      external qthresh,beta,exchange


      integer elnes,ipmin,ipmax,ipstep  !KJ added these variables 2-06
	character*12 f1,f2  !KJ dummy filenames  2-06
	integer kilast(10) !KJ stupid fix for stupid setup of this program 2-06
	integer ip !KJ local index 2-06
c     elnes : we're in EELS mode if elnes=1
c     ipmin,ipmax,ipstep : range of polarization components to calculate

c  !KJ Next section added to read ELNES variables 2-06     
c     read eels.inp
      elnes=0
      open(file='eels.inp',unit=3,status='old',err=900)
        read(3,*,err=900,end=900) 
        read(3,20,err=900,end=900) elnes
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,20,err=900,end=900) ipmin,ipstep,ipmax
  20  format (20i4)
      close(3)
      goto 901
900   continue
      elnes=0
901   continue
      if(elnes.eq.0) then
        ipstep=1
        ipmax=1
        ipmin=1
      endif
      do i=1,10
      kilast(i)=9
      enddo
      kilast(1)=7
c  !KJ end my changes



c  !KJ big loop for all polarization types 2-06 :
      do ip=ipmin,ipmax,ipstep
      if(ip.eq.1) then
	  f1(1:12)='chi.dat     '
	  f2(1:12)='xmu.dat     '
	elseif(ip.eq.10) then
	  f1(1:12)='chi10.dat   '
	  f2(1:12)='xmu10.dat   '
	elseif(ip.gt.1.and.ip.lt.10) then
	  f1(1:4)='chi0'
	  f1(5:5)= char(48+ip)
	  f1(6:12)='.dat   '
	  f2(1:4)='xmu0'
	  f2(5:5)= char(48+ip)
	  f2(6:12)='.dat   '
	else
	  stop 'crazy ip in ff2xmu'
	endif
c  !KJ end my changes


      isattype=0
      iwrite=0
      iwrite2=0
      iwrconv=0
      iwrcfile=0
      ipw=0
      ipwc=0
*      ipse=0
*      ipsk=0
      lowq=0
      icut=0
      brpole=.true.
*      brpole=.false.
*      wsigk=1.d0
      if (abs(ispec).eq.1) then
        cfile(1)= f2  !KJ 2-06    'xmu.dat     '
        ilast(1)= kilast(ip)  !KJ 2-06     7
        itype(1)=2
        cfile(2)= f1  !KJ 2-06    'chi.dat     '
        ilast(2)=7
        itype(2)=1
        ifiles=2
      elseif (abs(ispec).eq.0) then
c       Josh - Added convolution of xmu.dat for path expansion
        cfile(1)= f2  !KJ 2-06    'xmu.dat     '
        ilast(1)= kilast(ip)  !KJ 2-06   7
        itype(1)=2
c       Josh END
        cfile(2)= f1  !KJ 2-06    'chi.dat     '
        ilast(2)= kilast(ip)  !KJ 2-06   7
        itype(2)=1
        ik=2
        if (ipr6.eq.2.or.ipr6.eq.3) then
          open(unit=1,file='list.dat',status='old',iostat=ios)
          if (ios.gt.0) goto 10
 2        read(1,'(A)') cinp
          if (cinp(6:14).ne.'---------') goto 2
          read(1,'(A)') cinp
 3          read(1,'(A)',end=10) cinp
            ik=ik+1
            lfirst=1
            last=1
 4          if (cinp(lfirst:lfirst).eq.' ') then
              lfirst=lfirst+1
              last=lfirst
              goto 4
            endif
 5          if (cinp(last:last).ne.' ') then
              last=last+1
              goto 5
            endif
            if (last-lfirst.eq.3) then
              cpath='0'//cinp(lfirst:last-1)
            elseif (last-lfirst.eq.2) then
              cpath='00'//cinp(lfirst:last-1)
            elseif (last-lfirst.eq.1) then
              cpath='000'//cinp(lfirst:last-1)
            else
              cpath=cinp(lfirst:last-1)
            endif
*            open(unit=23,file='temp',status='unknown')
*            close(unit=23)
*            open(unit=23,file='temp',status='old')
*            read(23,*) cpath
*            close(unit=23,status='delete')
*            last=4
* 4          if (cpath(last:last).eq.' ') then
*              last=last-1
*              goto 4
*            endif
*            if (last.eq.3) then
*              cpath='0'//cpath(:last)
*            elseif (last.eq.2) then
*              cpath='00'//cpath(:last)
*            elseif (last.eq.1) then
*              cpath='000'//cpath(:last)
*            endif
            if (ipr6.eq.2) then
              cfile(ik+1)='chip'//cpath//'.dat'
              itype(ik+1)=1
            elseif (ipr6.eq.3) then
              cfile(ik+1)='feff'//cpath//'.dat'
              itype(ik+1)=3
            endif
            ilast(ik+1)=12
          goto 3
 10       close(unit=1)
        endif
        ifiles=ik+1
      endif
      do i=1,ifiles
        if (cfile(i).eq.cfname) iwrcfile=i
      enddo

      if (ipsk.ne.0) then
        iwrite=0
        iwrite2=0
      elseif (iwrite.ne.0) then
        iwrite2=0
      endif
        
* Open diagnostic and auxilliary output files.
      if (iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0) then
        open(unit=12,status='unknown',file='qpsf.dat')
        open(unit=13,status='unknown',file='satsf.dat')
      endif
      if (ipsk.ne.0.or.iwrite.ne.0) 
     2  open(unit=24,status='unknown',file='sigma.dat')
      if (ipwc.ne.0) open(unit=14,status='unknown',file='weightscl.dat')
      if (ipw.ne.0) open(unit=15,status='unknown',file='weights.dat')
      if (ipse.ne.0) 
     2 open(unit=18,status='unknown',file='selfenergy.dat')

      do ikk=1,ifiles
* Open input and output files.
        open(unit=21,status='old',
     2       file=cfile(ikk)(:ilast(ikk)),iostat=ios)
        open(unit=16,status='unknown',
     2       file='mbheader.temp')
*        if (itype(ikk).eq.1) then
*          open(unit=17,status='unknown',
*     2       file=cfile(ikk)(:ilast(ikk)-4)//'so2.dat')
*        endif
        if (ios.gt.0) goto 640

 12     read(21,'(A)') cfirst
        last=80
* Trim off excess white space.
 13     if (cfirst(last:last).eq.' ') then
          last=last-1
          goto 13
        endif
* Josh - Check that convolution has not been done on this file.
        if(cfirst(1:last).eq.'# Convoluted with A(omega).') then
          print '(A)', 'WARNING: '// cfile(ikk) 
          print '(A)', 'has already been convoluted ' //
     2          'in a previous calculation.'
          print '(A)', 'Rerun module 6 if you still '//
     2          'wish to proceed.'
          stop
        end if
* Find physical parameters of the material needed for computation.
        do isearch=1,last-7
          if(cfirst(isearch:isearch+6).eq.'Gam_ch=') then
            gchtx=cfirst(isearch+7:isearch+15)
          elseif(cfirst(isearch:isearch+6).eq.'Rs_int=') then
            rstx=cfirst(isearch+8:isearch+12)
          elseif(cfirst(isearch:isearch+4).eq.'Vint=') then
            vintx=cfirst(isearch+5:isearch+14)
          elseif(cfirst(isearch:isearch+2).eq.'Mu=') then
            mutx=cfirst(isearch+3:isearch+12)
          elseif(cfirst(isearch:isearch+2).eq.'kf=') then
            kftx=cfirst(isearch+3:isearch+11)
          endif
        enddo
* Write header lines.
        if ((iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0).and.ikk.eq.1) 
     2  then
          write(12,'(A)') cfirst(1:last)
          write(13,'(A)') cfirst(1:last)
        endif
        if (ikk.eq.1.and.(ipsk.ne.0.or.iwrite.ne.0))
     2    write(24,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipwc.ne.0) write(14,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipw.ne.0) write(15,'(A)') cfirst(1:last)
        write(16,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipse.ne.0) write(18,'(A)') cfirst(1:last)
        if (cfirst(6:14).ne.'---------') goto 12
* feffnnnn.dat files have more information lines after the line with
* dashes than the other file types.
        if (itype(ikk).eq.3) then
 14       read(21,'(A)') cfirst
          last=80
 15        if (cfirst(last:last).eq.' ') then
            last=last-1
            goto 15
          endif
* If first character is a number, keep it.  Otherwise it might be a
* comment character, so dump it.
          if (cfirst(1:1).eq.'1'.or.cfirst(1:1).eq.'2'
     2        .or.cfirst(1:1).eq.'3'.or.cfirst(1:1).eq.'4'
     3        .or.cfirst(1:1).eq.'5'.or.cfirst(1:1).eq.'6'
     4        .or.cfirst(1:1).eq.'7'.or.cfirst(1:1).eq.'8'
     5        .or.cfirst(1:1).eq.'9'.or.cfirst(1:1).eq.'0') then
            lfirst=1
          else
            lfirst=2
          endif
* Search for scattering half length, read into variable Rnn
          do isearch=1,last
            if(cfirst(isearch:isearch+3).eq.'reff') then
              open(unit=23,file='temp',status='unknown')
              write(23,*) cfirst(lfirst:last)
              close(unit=23)
              open(unit=23,file='temp',status='old')
              read(23,*) nleg,deg,Rnn
              close(unit=23,status='delete')
              Rnn=Rnn*aangstrom
* @# indicates end of comment lines.
            elseif(cfirst(isearch:isearch+1).eq.'@#') then
              goto 16
            endif
          enddo
          write(16,'(A)') cfirst(1:last)
          goto 14
        else 
          read(21,'(A)') cfirst
        endif
        last=80
 16     if (cfirst(last:last).eq.' ') then
          last=last-1
          goto 16
        endif
* Write column labels.
        if (ikk.eq.1) then
          if (iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0) then
            write(12,'(A)') '#  omega       delta ext    interf.'   
     2       // '      main ext     sat ext'
            write(13,'(A)') '#  omega       extrinsic    interf.'    
     2       // '      intrinsic    total'
          endif
          if (iwrite.ne.0.or.ipsk.ne.0)
     2    write(24,'(A)') '#   omega      Re(Sigma)    -Im(Sigma)'
          if (ipwc.ne.0)
     2    write(14,'(A)') '#   energy     ext q.p.   inter q.p. ext sat'
     3     // '    inter sat  intrin sat'
          if (ipw.ne.0)
     2    write(15,'(A)') '#   energy     ext q.p.   inter q.p. ext sat'
     3     // '    inter sat  intrin sat'
          if (ipse.ne.0)
     2    write(18,'(A)') '#   energy      Re(Sigma)   -Im(Sigma)'
     3     //'        Re(Z)        Im(Z)'
        endif
        write(16,'(A)') cfirst(1:last)
*        if (itype(ikk).eq.1) then
*          write(17,'(A)') '#     k          So2        phase shift'
*        endif
* Find numerical values of core hole lifetime, Wigner-Seitz radius, 
* interstitial potential, chemical potential, and Fermi momentum.
        open(unit=23,file='temp',status='unknown')
        write(23,*) gchtx
        write(23,*) rstx
        write(23,*) vintx
        write(23,*) mutx
        write(23,*) kftx
        close(unit=23)
        open(unit=23,file='temp',status='old')
        read(23,*) gammach
        read(23,*) rs
        read(23,*) vint
        read(23,*) cmu
        read(23,*) ckf
        close(unit=23,status='delete')
        gammach=(gammach/2.d0)*eV
        cmu=(cmu-vint)*eV
        vint=vint*eV
        ckf=ckf/aangstrom
        xreduc=1.d0
* Compute important material properties and constants.
        pi=dacos(-1.d0)
        qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
        ef=qf*qf/2.d0
        conc=3.d0/(4.d0*pi*(rs**3))
        omp=dsqrt(4.d0*pi*conc)
        adisp=2.d0*ef/3.d0
        ekp=ef
        qpk=qf
        acc=1.d-4

* Read pole expansion for epsilon^{-1}
        call rdeps(omp,nplmax,npl,plengy,oscstr,plbrd)
        epswt=0.d0
        do ipl=1,npl
          plwt(ipl)=dabs(oscstr(ipl)*plengy(ipl)**2/(omp**2))
          epswt=epswt+plwt(ipl)
        enddo
        open (unit=88,file='apl.dat', status='unknown')
        do ipl=1,npl
          write(88,'(5f10.5)') plengy(ipl)/eV,oscstr(ipl)*plengy(ipl)
        enddo
        close (88)

        sef0=0.d0
        do ipl=1,npl
          call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
          if (brpole) then
            call brsigma(0.d0,sse,sxise)
          else 
            call renergies(0.d0,sse)
          endif
          sef0=sef0+sse*wt
        enddo
        sef0=sef0+exchange(qf)

* Estimate threshold for plasmon creation.
        pthresh=qthresh(omp,adisp,ef,qf)
* Find minimal momentum grid.
        do i=1,10
          dp=(pthresh-qf)/10.d0
          pgrid(i)=qf+i*dp
        enddo
        do i=1,30
          dp=0.25d0*pthresh/30.d0
          pgrid(i+10)=pgrid(10)+i*dp
        enddo
        do i=1,10
          dp=0.75d0*pthresh/10.d0
          pgrid(i+40)=pgrid(40)+i*dp
        enddo
        do i=1,10
          dp=2.d0*pthresh/10.d0
          pgrid(i+50)=pgrid(50)+i*dp
        enddo
        pgrid(61)=5.d0*pthresh
        pgrid(62)=7.d0*pthresh
        pgrid(63)=1.d1*pthresh
        pgrid(64)=3.d1*pthresh
        pgrid(65)=1.d2*pthresh
        pgrid(66)=3.d2*pthresh
  
* Read information from files.
        if (itype(ikk).eq.1) then
          iasym=0
          do i=1,npts2
            epts2(i)=0.d0
            xk(i)=0.d0
            pk(i)=0.d0
            chi(i)=0.d0
            xmag(i)=0.d0
            phase(i)=0.d0
            phm2kr(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=25) xk(i),chi(i),xmag(i),phase(i)
            xk(i)=xk(i)/aangstrom
            chir(i)=xmag(i)*cos(phase(i))
            chii(i)=xmag(i)*sin(phase(i))
            j=i
          enddo
        elseif (itype(ikk).eq.2) then
          iasym=1
          do i=1,npts2
            e1(i)=0.d0
            epts2(i)=0.d0
            xk(i)=0.d0
            pk(i)=0.d0
            xmu(i)=0.d0
            xmu0(i)=0.d0
            chi(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=25) e1(i),epts2(i),xk(i),xmu(i),xmu0(i),chi(i)
            e1(i)=e1(i)*eV
            epts2(i)=epts2(i)*eV
            xktest=dabs(xk(i))
            if (xktest.lt.xktold) iedge=i
            xktold=xktest
            xk(i)=xk(i)/aangstrom
            j=i
          enddo
          dw=epts2(j)-epts2(j-1)
        elseif (itype(ikk).eq.3) then
          iasym=0
          do i=1,npts2
            epts2(i)=0.d0
            xk(i)=0.05d0*(i-1)/aangstrom
            pk(i)=0.d0
            chi(i)=0.d0
            xmag(i)=0.d0
            phase(i)=0.d0
            phm2kr(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=17) xk2(i),caph2(i),xmfeff2(i),phfeff2(i),
     2          redfac2(i),xlam2(i),realck2(i)
            xk2(i)=xk2(i)/aangstrom
            xmfeff2(i)=xmfeff2(i)*aangstrom
            xlam2(i)=xlam2(i)*aangstrom
            realck2(i)=realck2(i)/aangstrom
            j=i
          enddo
 17       ifeffct=j
* Interpolate data from feffnnnn.dat grid to the uniform grid used
* in chi.dat.
          do i=1,npts2
            do jj=1,j-1
              if (xk(i).ge.xk2(jj).and.xk(i).lt.xk2(jj+1)) then
                caph(i)=caph2(jj)+(caph2(jj+1)-caph2(jj))
     2                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                xmfeff(i)=xmfeff2(jj)+(xmfeff2(jj+1)-xmfeff2(jj))
     2                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                phfeff(i)=phfeff2(jj)+(phfeff2(jj+1)-phfeff2(jj))
     2                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                redfac(i)=redfac2(jj)+(redfac2(jj+1)-redfac2(jj))
     2                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                xlam(i)=xlam2(jj)+(xlam2(jj+1)-xlam2(jj))
     2                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                goto 18
              endif
            enddo
            if (xk(i).eq.xk2(jj)) then
              caph(i)=caph2(jj)
              xmfeff(i)=xmfeff2(jj)
              phfeff(i)=phfeff2(jj)
              redfac(i)=redfac2(jj)
              xlam(i)=xlam2(jj)
            endif
* Compute exafs signal Chi from feffnnnn.dat info.
 18         xmag(i)=(deg*xmfeff(i)*redfac(i)*exp(-2.d0*Rnn/xlam(i)))
     2              /(xk(i)*Rnn**2)
            phm2kr(i)=phfeff(i)+caph(i)
            phase(i)=phm2kr(i)+2.d0*xk(i)*Rnn
            chir(i)=xmag(i)*cos(phase(i))
            chii(i)=xmag(i)*sin(phase(i))
          enddo
          xmag(1)=xmag(2)+(xk(1)-xk(2))*(xmag(3)-xmag(2))/(xk(3)-xk(2))
          chir(1)=xmag(1)*cos(phase(1))
          chii(1)=xmag(1)*sin(phase(1))
          j=npts2
        endif
 25     continue

        close(unit=21)
        close(unit=16)
        open(unit=21,status='old',file='mbheader.temp',iostat=ios)
        if (ios.gt.0) goto 640
        open(unit=16,status='unknown',
     2       file=cfile(ikk)(:ilast(ikk)))
*	Josh - added to header so that we can check if a previous convolution
*	has been performed.
        write(16,'(A)') '# Convoluted with A(omega).'
 27       read(21,'(A)',end=35) cfirst
          last=80
 30       if (cfirst(last:last).eq.' ') then
            last=last-1
            goto 30
          endif
          write(16,'(A)') cfirst(1:last)
          goto 27
 35     close(unit=21,status='delete')
          
* Need to find photoelectron momentum.
* First, find kinetic energy.
        do i=1,j
          if (xk(i).ge.0.d0) then
            ekp=xk(i)**2/2.d0+cmu
          else
            ekp=-xk(i)**2/2.d0+cmu
          endif
          ekpg(i)=ekp
          if (itype(ikk).eq.1.or.itype(ikk).eq.3) epts2(i)=ekp
* Make arrays of 0 order self energies and 0 order estimates for momentum.
          if (ekp.ge.0.d0) then
            qpk=sqrt(qf**2+2.d0*(ekp-fmu))
            xpkg(i)=qpk
            se0=0.d0
            do ipl=1,npl
              call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
              if (brpole) then
                call brsigma(0.d0,sse,sxise)
              else
                call renergies(0.d0,sse)
              endif
              se0=se0+sse*wt
            enddo
            seg(i)=se0+exchange(qpk)
          endif
        enddo
* Estimate momentum derivative renormalization constant.
        do i=1,j
          if (ekpg(i).ge.0) then
            if (i.ne.1.and.i.ne.j) then
              dsig=seg(i+1)-seg(i-1)
              dk2=xpkg(i+1)**2/2.d0-xpkg(i-1)**2/2.d0
            elseif (i.eq.1) then
              dsig=seg(i+1)-seg(i)
              dk2=xpkg(i+1)**2/2.d0-xpkg(i)**2/2.d0
            else
              dsig=seg(i)-seg(i-1)
              dk2=xpkg(i)**2/2.d0-xpkg(i-1)**2/2.d0
            endif
            zkk=1.d0/(1.d0+dsig/dk2)
* Find array of photoelectrom momenta.
            pk(i)=sqrt(xpkg(i)**2-2.d0*zkk*(seg(i)-sef0))
          endif
        enddo
* if running convolution is written, choose momentum point closest to desired energy
        if (iwrcfile.ne.0) then
          den=abs(ekpg(1)-cen)
          iwrconv=1
          do i=2,j
            den2=abs(ekpg(i)-cen)
            if (den2.lt.den) then
              den=den2
              iwrconv=i
            endif
          enddo
        endif
  
* Pad out the ends of the arrays so the convolution does not screw up the
* end data points.
        if (itype(ikk).eq.1.or.itype(ikk).eq.3) then
          dw=epts2(j)-epts2(j-1)
          do i=j+1,npts2
            epts2(i)=epts2(i-1)+dw
          enddo
        elseif (itype(ikk).eq.2) then
          dw=epts2(j)-epts2(j-1)
          do i=j,npts2
            e1(i)=e1(i-1)+dw
            epts2(i)=epts2(i-1)+dw
            xmu0(i)=xmu0(j)
            xmu(i)=xmu0(j)
            chi(i)=0.d0
          enddo
          call mkrmu(xmu,xmu0,rmu,epts2,npts2)
*          call mkrmu(xmu0,xmu0,rmu0,epts2,npts2)
          do i=1,npts2
            ximu(i)=xmu(i)-xmu0(i)
          enddo
        endif

* Find the spectral function and associated data.
        iout=0
* If specfunct.dat does not exist, or if it is for a material with 
* different electron gas properties, set flag 'iout' to recompute 
* spectral function.  Otherwise, read in spectral function.
        open(unit=23,file='specfunct.dat',status='old',
     2       access='sequential',form='unformatted',iostat=ios)
        if (ios.gt.0) then
          iout=1
        else
          read(23) rst,gcht,iasymt,isattypt,lowqt,nplt
          read(23) (plengyt(ii), ii=1,nplmax)
          read(23) (plbrdt(ii), ii=1,nplmax)
          read(23) (plwtt(ii), ii=1,nplmax)
          read(23) ((sfinfo(jj,ii), jj=1,nqpts), ii=1,8)
          read(23) ((wgts(jj,ii), jj=1,nqpts), ii=1,8)
          read(23) ((emsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((essf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xmsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xssf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xissf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((escsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((engrid(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        endif
        close(unit=23)
        if (rst.ne.rs) iout=1
        if (gcht.ne.gammach) iout=1
        if (iasymt.ne.iasym) iout=1
        if (lowqt.ne.lowq) iout=1
        if (isattypt.ne.isattype) iout=1
        if (nplt.ne.npl) iout=1
        do jj=1,npl
          if (plengy(jj).ne.plengyt(jj)) iout=1
          if (plbrd(jj).ne.plbrdt(jj)) iout=1
          if (plwt(jj).ne.plwtt(jj)) iout=1
        enddo
        do jj=1,nqpts
          if (real(sfinfo(jj,1)).ne.real(pgrid(jj))) iout=1
        enddo
        if (iwrite.ne.0.and.ikk.eq.1) iout=1

        npi=0
        if (iout.eq.1) then
          write(6,*) 'computing spectral function'
        endif
        if(ipsk.ne.0.and.ikk.eq.1) then
          jj=0
          call mkspectf(rs,wsigk*qf,gammach,xreduc,
     2                  epts,spectf,weights,isattype,
     3                  npl,nplmax,plengy,plwt,plbrd,brpole)

*          do ipl=1,npl
*            write(6,*) ipl,plengy(ipl),plbrd(ipl),plwt(ipl),
*     2        plengy(ipl)+sef0-se,sef0,se
*          enddo
  
        endif
        do jj=1,nqpts
          if (iout.eq.1) then
* Spectral function has not yet been computed, or was computed for the
* wrong material.  Compute spectral function.
            call mkspectf(rs,pgrid(jj),gammach,xreduc,
     2                    epts,spectf,weights,isattype,
     3                  npl,nplmax,plengy,plwt,plbrd,brpole)
            write(6,'(I3,A10)') int(real(jj)/real(nqpts)*100.0),
     &                '% computed'
            sfinfo(jj,1)=pgrid(jj)
            sfinfo(jj,2)=ekp
            sfinfo(jj,3)=ek
            sfinfo(jj,4)=se
            sfinfo(jj,5)=ce
            sfinfo(jj,6)=width
            sfinfo(jj,7)=z1
            sfinfo(jj,8)=z1i
            do ii=1,8
              wgts(jj,ii)=weights(ii)
            enddo
            do ii=1,nsfpts
              emsf(jj,ii)=spectf(1,ii)
              essf(jj,ii)=spectf(2,ii)
              xmsf(jj,ii)=spectf(3,ii)
              xssf(jj,ii)=spectf(4,ii)
              xissf(jj,ii)=spectf(5,ii)
              escsf(jj,ii)=spectf(8,ii)
              engrid(jj,ii)=epts(ii)
            enddo
          else
* Correct spectral function has already been computed.
            ekp=sfinfo(jj,2)
            ek=sfinfo(jj,3)
            se=sfinfo(jj,4)
            ce=sfinfo(jj,5)
            width=sfinfo(jj,6)
            z1=sfinfo(jj,7)
            z1i=sfinfo(jj,8)
            do ii=1,8
              weights(ii)=wgts(jj,ii)
            enddo
          endif
        enddo
        write(6,*) 'convoluting file '
     2           //cfile(ikk)(:ilast(ikk))
* Interpolate spectral function onto uniform momentum grid.
        do jj=1,j
          do i=1,nqpts-1
            if (pk(jj).ge.pgrid(i).and.pk(jj).lt.pgrid(i+1)) then
              delta=(pk(jj)-pgrid(i))/(pgrid(i+1)-pgrid(i))
              do ii=1,nsfpts
                epts(ii)=engrid(i,ii)
     2               +delta*(engrid(i+1,ii)-engrid(i,ii))
                spectf(1,ii)=emsf(i,ii)+delta*(emsf(i+1,ii)-emsf(i,ii))
                spectf(2,ii)=essf(i,ii)+delta*(essf(i+1,ii)-essf(i,ii))
                spectf(3,ii)=xmsf(i,ii)+delta*(xmsf(i+1,ii)-xmsf(i,ii))
                spectf(4,ii)=xssf(i,ii)+delta*(xssf(i+1,ii)-xssf(i,ii)) 
                spectf(5,ii)=xissf(i,ii)
     2               +delta*(xissf(i+1,ii)-xissf(i,ii))
                spectf(6,ii)=essf(i,ii)+xissf(i,ii)-2.d0*xssf(i,ii)
     2               +delta*(essf(i+1,ii)+xissf(i+1,ii)
     3                     -2.d0*xssf(i+1,ii)
     4                     -(essf(i,ii)+xissf(i,ii)-2.d0*xssf(i,ii)))
                spectf(7,ii)=essf(i,ii)-escsf(i,ii)
     2               +delta*(essf(i+1,ii)-escsf(i+1,ii)
     3                     -(essf(i,ii)-escsf(i,ii)))
                spectf(8,ii)=escsf(i,ii)
     2               +delta*(escsf(i+1,ii)-escsf(i,ii))
              enddo
              do ii=1,8
                weights(ii)=wgts(i,ii)+delta*(wgts(i+1,ii)-wgts(i,ii))
              enddo
              se=sfinfo(i,4)+delta*(sfinfo(i+1,4)-sfinfo(i,4))
              ce=sfinfo(i,5)+delta*(sfinfo(i+1,5)-sfinfo(i,5))
              width=sfinfo(i,6)+delta*(sfinfo(i+1,6)-sfinfo(i,6))
              z1=sfinfo(i,7)+delta*(sfinfo(i+1,7)-sfinfo(i,7))
              z1i=sfinfo(i,8)+delta*(sfinfo(i+1,8)-sfinfo(i,8))
            endif
          enddo
* Take special care with the endpoints.
          if (pk(jj).ge.pgrid(nqpts)) then
            do ii=1,nsfpts
              epts(ii)=engrid(nqpts,ii)
              spectf(1,ii)=emsf(nqpts,ii)
              spectf(2,ii)=essf(nqpts,ii)
              spectf(3,ii)=xmsf(nqpts,ii)
              spectf(4,ii)=xssf(nqpts,ii)
              spectf(5,ii)=xissf(nqpts,ii)
              spectf(6,ii)=essf(nqpts,ii)+xissf(nqpts,ii)
     2                     -2.d0*xssf(nqpts,ii)
              spectf(7,ii)=essf(nqpts,ii)-escsf(nqpts,ii)
              spectf(8,ii)=escsf(nqpts,ii)
            enddo
            do ii=1,8
              weights(ii)=wgts(nqpts,ii)
            enddo
            se=sfinfo(nqpts,4)
            ce=sfinfo(nqpts,5)
            width=sfinfo(nqpts,6)
            z1=sfinfo(nqpts,7)
            z1i=sfinfo(nqpts,8)
          elseif (pk(jj).lt.pgrid(1)) then
            do ii=1,nsfpts
              epts(ii)=engrid(nqpts,ii)
              spectf(1,ii)=emsf(1,ii)
              spectf(2,ii)=essf(1,ii)
              spectf(3,ii)=xmsf(1,ii)
              spectf(4,ii)=xssf(1,ii)
              spectf(5,ii)=xissf(1,ii)
              spectf(6,ii)=essf(1,ii)+xissf(1,ii)-2.d0*xssf(1,ii)
              spectf(7,ii)=essf(1,ii)-escsf(1,ii)
              spectf(8,ii)=escsf(1,ii)
            enddo
            do ii=1,8
              weights(ii)=wgts(1,ii)
            enddo
            se=sfinfo(1,4)
            ce=sfinfo(1,5)
            width=sfinfo(1,6)
            z1=sfinfo(1,7)
            z1i=sfinfo(1,8)
          endif
* Interpolate spectral function onto uniform energy grid.
          call interpsf(npts,epts,wpts,spectf,cspec)
* Write information to diagnostic and supplementary files.
          if (ikk.eq.1.and.ipwc.ne.0) then
            write(14,701) epts2(jj),weights(1)+weights(8),
     2                  weights(3),weights(7),weights(5),weights(6)
          endif
          if (ikk.eq.1.and.ipw.ne.0) then
            write(15,701) epts2(jj),weights(1),weights(3)/2.d0,
     2                  weights(4),weights(5),weights(6)
          endif
          if (ikk.eq.1.and.ipse.ne.0) then
            write(18,500) epts2(jj),se,width,z1,z1i
          endif
          do i=1,nsfpts
            w=epts(i)
            emain=spectf(1,i)
            esat =spectf(2,i)
            xmain=spectf(3,i)
            xsat =spectf(4,i)
            xisat =spectf(5,i)
            sat  =spectf(6,i)
            eclip=spectf(7,i)
            esatr=spectf(8,i)
            if (jj.eq.iwrite2.and.ikk.eq.1) then
              write(12,500) w,emain,xmain,eclip,esatr
              write(13,500) w,esat,-2.d0*xsat,xisat,
     2                      (esat+xisat-2.d0*xsat)
            endif
          enddo
          wp=epts2(jj)
* Begin convolution.
          if (itype(ikk).eq.1.or.itype(ikk).eq.3) then
            if (jj.eq.iwrconv.and.ikk.eq.iwrcfile) then
              intout=1
              open (unit=28,file='realint.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu,gammach,npts2,epts2,chir,npts,
     2                wpts,cspec,weights,xchir,phchir,iasym,
     3                1,intout,omp)
            if (jj.eq.iwrconv.and.ikk.eq.iwrcfile) then
              intout=1
              open (unit=28,file='imagint.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu,gammach,npts2,epts2,chii,npts,
     2                wpts,cspec,weights,xchii,phchii,iasym,
     3                1,intout,omp)
* Find real and imaginary many body EXAFS signal chi.
            chirr=xchir*cos(phchir)-xchii*sin(phchii)
            chiii=xchii*cos(phchii)+xchir*sin(phchir)
* Compute phase shift, remove jumps of 2 pi.
            phaseshft=datan2(chiii,chirr)
            if(phaseshft-phshftold.gt.5.d0) then
              npi=npi+2
            elseif(phaseshft-phshftold.lt.-5.d0) then
              npi=npi-2
            endif
            phshftold=phaseshft
            phaseshft=phaseshft-pi*npi
            if (itype(ikk).eq.1) then
* Write output file of many body EXAFS signal.
              write(16,630) xk(jj)*aangstrom,
     2                chiii,dsqrt(chirr**2+chiii**2),
     3                phaseshft,
     4                phaseshft+phm2kr(jj)-phase(jj)
            endif
 630        format (1x, f10.4, 3x, 4(1pe13.6,1x))
* Find amplitude reduction and phase shift.
            so2mag=dsqrt(chirr**2+chiii**2)/xmag(jj)
            phaseshft=datan2(chiii,chirr)
            if(phaseshft-phshftold.gt.5.d0) then
              npi=npi+2
            elseif(phaseshft-phshftold.lt.-5.d0) then
              npi=npi-2
            endif
            phshftold=phaseshft
            phaseshft=phaseshft-pi*npi-phase(jj)
            s02list(jj)=so2mag
            phlist(jj)=phaseshft
*            if (itype(ikk).eq.1) then
** Write output file of amplitude reduction and phase shift.
*              write(17,500) xk(jj)*aangstrom,
*     2                so2mag,phaseshft
*            endif
          elseif (itype(ikk).eq.2) then
* Find many body convolution on absorption signal, embedded atom background,
* and fine structure.
*            call sfconv(wp,cmu+vint,gammach,npts2,epts2,chi,npts,
*     2                wpts,cspec,weights,xchi2,phchi,iasym,1,0,omp)
            if (iasym.eq.0) then
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='imagint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,ximu,npts,
     2            wpts,cspec,weights,ximu2,phmu,iasym,icut,intout,omp)
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='realint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,rmu,npts,
     2            wpts,cspec,weights,rmu2,phrmu,iasym,icut,intout,omp)
            else
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='realint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,xmu,npts,
     2            wpts,cspec,weights,xmu2,phrmu,iasym,icut,intout,omp)
            endif
            if (jj.eq.iwrconv) then
              intout=1
              open (unit=28,file='xmu0int.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu+vint,gammach,npts2,epts2,xmu0,npts,
     2            wpts,cspec,weights,xmu02,phmu0,iasym,icut,intout,omp)
            if (iasym.eq.0) then
              xmu2=ximu2*cos(phmu)+rmu2*sin(phrmu)+xmu02
            endif
            write(16,800)
     2                e1(jj)/eV,epts2(jj)/eV,xk(jj)*aangstrom,
!     2                xmu2,xmu02,(xmu2-xmu02)/xmu02    !KJ this is the original instruction  3-06
     2                xmu2,xmu02,(xmu2-xmu02)   !KJ removed normalization of chi   3-06  
          endif
        enddo
        if (itype(ikk).eq.3) then
          dk=0.05d0
* Average over nearby points on uniform momentum grid to find
* values for coarser grid in feffnnnn.dat files.  Use a trianguar
* weighting function to ensure each point contributes its full weight
* to the average.
          do jj=1,ifeffct
            s02sum=0.d0
            dphsum=0.d0
            xnorm=0.d0
            do i=1,j
              if (xk2(jj).eq.xk(i)) then
                ww=1.d0
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              elseif (xk(i).gt.xk2(jj-1).and.xk(i).le.xk2(jj)
     2                .and.xk2(jj-1).ne.xk2(jj)) then
                ww=(xk(i)-xk2(jj-1))/(xk2(jj)-xk2(jj-1))
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              elseif(xk(i).gt.xk2(jj).and.xk(i).lt.xk2(jj+1)
     2                .and.xk2(jj+1).ne.xk2(jj)) then
                ww=(xk2(jj+1)-xk(i))/(xk2(jj+1)-xk2(jj))
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              else
                ww=0.d0
              endif
            enddo
* Find many body values of phase and reduction.
            redfac2(jj)=redfac2(jj)*s02sum/xnorm
            caph2(jj)=caph2(jj)+dphsum/xnorm
* Write convoluted feffnnnnc.dat.
            write(16,400) xk2(jj)*aangstrom,caph2(jj),
     2          xmfeff2(jj)/aangstrom,phfeff2(jj),redfac2(jj),
     2          xlam2(jj)/aangstrom,realck2(jj)*aangstrom
  400       format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
     1                            2(1pe11.4,1x))
          enddo
        endif
* Write out the spectral function to a data file.
        open(unit=23,file='specfunct.dat',status='unknown',
     2       access='sequential',form='unformatted',iostat=ios)
        write(23) rs,gammach,iasym,isattype,lowq,npl
        write(23) (plengy(ii), ii=1,nplmax)
        write(23) (plbrd(ii), ii=1,nplmax)
        write(23) (plwt(ii), ii=1,nplmax)
        write(23) ((sfinfo(jj,ii), jj=1,nqpts), ii=1,8)
        write(23) ((wgts(jj,ii), jj=1,nqpts), ii=1,8)
        write(23) ((emsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((essf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xmsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xssf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xissf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((escsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((engrid(jj,ii), jj=1,nqpts), ii=1,nsfpts)
* Do a little housecleaning.
        close(unit=23)
        close(unit=21)
        close(unit=16)
        close(unit=17)
 640    continue
      enddo
 500  format(1x,5(e12.5,1x))
 700  format(1x,7(f10.5,1x))
 701  format(1x,e10.5,1x,6(f10.5,1x))
 800  format(1x,2f11.3,f8.3,1p,3e13.5)
* End of program.


      enddo  !KJ end of loop do ip=ipmin,ipmax,ipstep  2-06


c     sub-program so2conv
c     stop
      return

      end
      subroutine res02 ( mso2conv, ispec, ipr6, ipse, ipsk, wsigk, cen,
     2                   cfname)

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

        integer  mso2conv, ipse, ipsk, ispec, ipr6
        double precision  wsigk, cen
        character*12 cfname

c     Local stuff
      character*512 slog
      character*80 head(nheadx)
      dimension lhead(nheadx)

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)


c     read s02.inp
      open (file='s02.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,20)  mso2conv, ipse, ipsk
        read (3,10)  slog
        read (3,30)  wsigk, cen
        read (3,10)  slog
        read (3,20)  ispec, ipr6
        read (3,10)  slog
        read (3,10)  cfname
      close(3)

      return
      end
      subroutine rdeps (omp,nplmax,npl,plengy,oscstr,plbrd)

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      integer  npl,nplmax,ipl
      double precision  plengy(nplmax),oscstr(nplmax),plbrd(nplmax)

c     Local stuff
      character*512 slog
      character*80 head(nheadx),line
      dimension lhead(nheadx)

*     initialize
      do ipl=1,nplmax
        plengy(ipl)=0.d0
        oscstr(ipl)=0.d0
        plbrd(ipl)=0.d0
      enddo

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)

c     read s02.inp
      open (file='exc.dat', unit=3, status='old',iostat=ios)
      if (ios.eq.0) then
        npl=0
 100      read (3,'(a)',end=200)  line
          istart=1
 110      if (line(istart:istart).eq.' ') then
            istart=istart+1
            goto 110
          endif
          if (line(istart:istart).ne.'#') then
            npl=npl+1
            read (line,*) plengy(npl),plbrd(npl),oscstr(npl)
          endif
          goto 100
 200    continue
        do ipl=1,npl
          plengy(ipl)=plengy(ipl)/hart
          plbrd(ipl)=plbrd(ipl)/hart
        enddo
      else
        npl=1
        plengy(1)=omp
        plbrd(1)=0.001d0*omp
        oscstr(1)=1.d0
        open (file='exc.dat', unit=3, status='unknown')
        write(3,30) plengy(1)*hart,plbrd(1)*hart,oscstr(1)
      endif
      close(3)
      


*      write(6,*) '***',npl,'*',plengy(1),'*',plbrd(1),'*',oscstr(1),'***'


      return
      end
      subroutine plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
* Sets calculation parameters for a given pole in epsilon^{-1}
*        ipl - pole to set parameters for
*        nplmax - array dimensioning (maximum number of poles)
*        plengy - energy of poles (array)
*        plwt - weight of poles (array)
*        ompl - energy of selected pole
*        wt - weight of selected pole
*        adisp - dispersion parameter for selected pole
*        conc - electron density that would give plasma frequency
*               equal to pole energy
*        rs - Wigner Seitz radius that would give plasma frequency
*               equal to pole energy
*        ef - Fermi energy that would give plasma frequency
*               equal to pole energy
*        qf - Fermi momentum that would give plasma frequency
*               equal to pole energy
      implicit none
      integer ipl,nplmax
      double precision plengy(nplmax),plwt(nplmax),plbrd(nplmax)
      double precision ompl,wt,brd
*      double precision pi,adisp,conc,rs,ef,qf
      ompl=plengy(ipl)
      wt=plwt(ipl)
      brd=plbrd(ipl)
*      adisp=(0.75d0*pi/ompl**2)**(2.d0/3.d0)/3.d0
*      rs=(3.d0/ompl**2)**(1.d0/3.d0)
*      conc=3.d0/(4.d0*pi*(rs**3))
*      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
*      ef=qf*qf/2.d0
      return
      end

      subroutine ppset(rs,pi,qf,ef,omp)
* Sets electron gas parameters for a given Wigner-Seitz radius rs.
*        rs - Wigner Seitz radius that would give plasma frequency
*               equal to pole energy
*        omp - plasma frequency
*        ef - Fermi energy that would give plasma frequency
*               equal to pole energy
*        qf - Fermi momentum that would give plasma frequency
*               equal to pole energy
      implicit none
      double precision rs,qf,ef,conc,omp,pi
      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
      ef=qf*qf/2.d0
      conc=3.d0/(4.d0*pi*(rs**3))
      omp=dsqrt(4.d0*pi*conc)
      return
      end

c///////////////////////////////////////////////////////////////////////
c Distribution:  COMMON 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: COMMON 1.0
c      COMMON 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
      i = istrln(fname)
      j = istrln(mod)
      write(slog,100)  fname(1:i), mod(1:j)
      call wlog(slog)

  100 format (' Error opening file, ', a, 
     2        ' in module ', a)

      call wlog(' Fatal error')
      call par_stop('CHOPEN')
      end
      subroutine fixdsp (dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect code.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc0(251), dpc0(251)
      dimension dgcx(nrptx), dpcx(nrptx)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.
      do 100  i = 251, 1, -1
         if ( abs(dgc0(i)) .ge. 1.0d-11 .or. 
     1        abs(dpc0(i)) .ge. 1.0d-11 )  then
            imax = i
            goto 16
         endif
  100 continue
      call wlog(' Should never see this line from sub fixdsp')
   16 continue
c     jmax is the first point where both dpc and dgc are zero in
c     the original grid
      jmax = imax + 1
      if (jmax.gt.251) jmax = 251

      delta = dxorg
      do 10  j = 1, jmax
         xorg(j) = xxx(j)
   10 continue
      rmax = rrr(jmax)

c     How far out do we go in the new grid?  To the last new grid
c     point before jmax.  Everything will be zero beyond jmax.
      delta = dxnew
      jnew = jjj(rmax)
      do 20  j = 1, jnew
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of rmax
      do 30  j = 1, jnew
         call terp (xorg, dgc0,  jmax, 3, xnew(j), dgcx(j))
         call terp (xorg, dpc0,  jmax, 3, xnew(j), dpcx(j))
   30 continue

c     and zero the arrays past rmax
      do 32  j = jnew+1, nrptx
         dgcx(j) = 0
         dpcx(j) = 0
   32 continue

      return
      end
      subroutine fixdsx (iph, dxorg, dxnew, dgc, dpc, dgcn, dpcn)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect and phase codes.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.

      delta = dxorg
      do 10  j = 1, 251
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      do 20  j = 1, nrptx
         xnew(j) = xxx(j)
   20 continue

      do 200 iorb = 1, 30
         imax = 0
         do 100  i = 251, 1, -1
            if ( abs(dgc(i,iorb,iph)) .ge. 1.0d-11 .or. 
     1           abs(dpc(i,iorb,iph)) .ge. 1.0d-11 )  then
               imax = i
               goto 16
            endif
  100    continue
   16    continue
         if (imax .eq. 0) then
            jnew = 0
            goto 35
         endif
c        jmax is the first point where both dpc and dgc are zero in
c        the original grid
         jmax = imax + 1
         if (jmax .gt. 251) jmax = 251

         delta = dxorg
         rmax = rrr(jmax)

c        How far out do we go in the new grid?  To the last new grid
c        point before jmax.  Everything will be zero beyond jmax.
         delta = dxnew
         jnew = jjj(rmax)

c        interpolate to new grid using x, only inside of rmax
         do 30  j = 1, jnew
            call terp(xorg,dgc(1,iorb,iph),jmax,3, xnew(j),dgcn(j,iorb))
            call terp(xorg,dpc(1,iorb,iph),jmax,3, xnew(j),dpcn(j,iorb))
   30    continue

c        and zero the arrays past rmax
   35    do 40  j = jnew+1, nrptx
            dgcn(j,iorb) = 0
            dpcn(j,iorb) = 0
   40    continue
  200 continue

      return
      end
      subroutine fixvar (rmt, edens, vtot, dmag,
     1                   vint, rhoint, dxorg, dxnew, jumprm,
     2                   vjump, ri, vtotph, rhoph, dmagx)

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}


      dimension edens(251), vtot (251), dmag(251)
      dimension vtotph(nrptx), rhoph(nrptx), dmagx(nrptx)
      dimension ri(nrptx)
      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     PHASE needs
c     vtot = total potential including gs xcorr, no r**2
c     edens = rho, charge density, no factor of 4*pi, no r**2
c     From overlapping, vtot = potential only, ok as is
c                       edens = density*4*pi, so fix this here.
c     ri = r grid through imt+1

c     Only values inside the muffin tin are used, except that XCPOT
c     (in PHASE) uses values at imt+1 and requires these to be the
c     interstitial values.  So set the last part of the arrays to
c     interstitial values...

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
      
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

      delta = dxorg
      jmtorg = jjj(rmt)
      jriorg = jmtorg + 1
      jrior1 = jriorg + 1
      do 10  j = 1, jrior1
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      jmtnew = jjj(rmt)
      jrinew = jmtnew + 1
      jrine1 = jrinew + 1
      do 20  j = 1, jrine1
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of muffintin
c     jri (first interstitial point) must be set to interstitial value
      do 30  j = 1, jrinew
         call terp (xorg, vtot,  jriorg, 3, xnew(j), vtotph(j))
         call terp (xorg, edens, jrior1, 3, xnew(j), rhoph(j))
         call terp (xorg, dmag,  jrior1, 3, xnew(j), dmagx(j))
   30 continue

      if (jumprm .eq. 1) then
         xmt = log(rmt)
         call terp (xorg, vtot,  jriorg, 3, xmt, vmt)
         vjump = vint - vmt
      endif
      if (jumprm .gt. 0) then
         do 90  j = 1, jrinew
            vtotph(j) = vtotph(j) + vjump
   90    continue
      endif

      delta = dxnew
      do 180  j = 1, nrptx
         ri(j) = rrr(j)
  180 continue
      do 190  j = 1, jrinew
         rhoph(j) = rhoph(j)/(4*pi)
  190 continue
      do 200  j = jrinew+1, nrptx
         vtotph(j) = vint
         rhoph(j) = rhoint/(4*pi)
c fix later : need to calculate interstitial dmint
c      want interpolation beyond mt also
         dmagx(j) = 0.0d0
  200 continue

      return
      end
      subroutine getorb (iz, ihole, xion, iunf, norb, norbco, iorb,
     1                  iholep, nqn, nk, xnel, xnval, xmag)
c     Gets orbital data for chosen element.  Input is:
c       iz - atomic number of desired element,
c       ihole - index of core-hole orbital
c       xion  - ionicity (usually zero)
c     other arguments are output.
c       norb - total number of orbitals
c       norbco - number of core orbitals
c       iorb - index of orbital for making projections (last occupied)
c       iholep - index of core hole orbital in compacted list
c       nqn - principal quantum number for each orbital
c       nk - quantum number kappa for each orbital
c       xnel - occupation for each orbital
c       xnval - valence occupation for each orbital
c       xmag - spin magnetization for each orbital
c     Feel free to change occupation numbers for element of interest.
c     ival(i) is necessary only for partly nonlocal exchange model.
c     iocc(i) and ival(i) can be fractional
c     But you have to keep the sum of iocc(i) equal to nuclear charge.
c     Also ival(i) should be equal to iocc(i) or zero. 
c     Otherwise you have to change this subroutine or contact authors 
c     for help.

      implicit double precision (a-h, o-z)

c     Written by Steven Zabinsky, July 1989
c     modified (20 aug 1989)  table increased to at no 99
c     Recipe for final state configuration is changed. Valence
c     electron occupations are added. ala 17.1.1996

c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1

      dimension nqn(30), nk(30), xnel(30), xnval(30), xmag(30)
      dimension kappa (29)
      real iocc, ival, ispn
      dimension iocc (100, 29), ival (100, 29), ispn (100, 29)
      dimension nnum (29), iorb(-4:3)
      character*512 slog

c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/

c     principal quantum number (energy eigenvalue)
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/

c     occupation of each level for z = 1, 99
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 2,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 4,i),i=1,29)  /0,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 4,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 5,i),i=1,29)  /0,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 5,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

c     data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (iocc( 6,i),i=1,29)  /2,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
c     data (ival( 6,i),i=1,29)  /0,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (ival( 6,i),i=1,29)  /0,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 6,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 7,i),i=1,29)  /0,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 7,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 8,i),i=1,29)  /0,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 8,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 9,i),i=1,29)  /0,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 9,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(10,i),i=1,29)  /0,0,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(10,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(11,i),i=1,29)  /0,0,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(11,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(12,i),i=1,29)  /2,2,2,4,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(12,i),i=1,29)  /0,0,0,0,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(12,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(13,i),i=1,29)  /0,0,0,0,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(13,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(14,i),i=1,29)  /0,0,0,0,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(14,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(15,i),i=1,29)  /0,0,0,0,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(15,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(16,i),i=1,29)  /0,0,0,0,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(16,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(17,i),i=1,29)  /0,0,0,0,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(17,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(18,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(18,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(19,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(19,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(20,i),i=1,29)  /0,0,0,0,0,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(20,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(21,i),i=1,29)  /0,0,0,0,0,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(21,i),i=1,29)  /0,0,0,0,0,  0,0,1,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(22,i),i=1,29)  /0,0,0,0,0,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(22,i),i=1,29)  /0,0,0,0,0,  0,0,2,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(23,i),i=1,29)  /0,0,0,0,0,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(23,i),i=1,29)  /0,0,0,0,0,  0,0,3,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(24,i),i=1,29)  /0,0,0,0,0,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(24,i),i=1,29)  /0,0,0,0,0,  0,0,4,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(26,i),i=1,29)  /0,0,0,0,0,  0,0,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(26,i),i=1,29)  /0,0,0,0,0,  0,0,2,2,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(27,i),i=1,29)  /0,0,0,0,0,  0,0,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(27,i),i=1,29)  /0,0,0,0,0,  0,0,0,3,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(28,i),i=1,29)  /0,0,0,0,0,  0,0,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(28,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(29,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(29,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(30,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(30,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(31,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(31,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(32,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(32,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,1,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,3,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ival(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ispn(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ispn(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,1,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           2,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           2,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           2,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           3,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           4,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           4,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           5,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           6,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           2,0,0,0,0,  0,0,0,2,0,  0,0,0,0/

      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           7,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ival(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ispn(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  0,0,0,0,0,  0,0,0,0/

      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (ival(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  0,0,0,2,0,  0,0,0,0/
      data (ispn(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,3,  0,0,0,0,0,  0,0,0,0/

      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (ival(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,0,0,0,4,  0,0,0,2,0,  0,0,0,0/
      data (ispn(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  0,0,0,0,0,  0,0,0,0/

      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (ival(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,2,0,  0,0,0,0/
      data (ispn(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,0,0,  0,0,0,0/

      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (ival(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  2,0,0,2,0,  0,0,0,0/
      data (ispn(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  2,0,0,0,0,  0,0,0,0/

      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (ival(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  3,0,0,2,0,  0,0,0,0/
      data (ispn(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  3,0,0,0,0,  0,0,0,0/

      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (ival(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  5,0,0,1,0,  0,0,0,0/
      data (ispn(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  2,0,0,0,0,  0,0,0,0/

      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (ival(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,1,0,  0,0,0,0/
      data (ispn(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  1,0,0,0,0,  0,0,0,0/

      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (ival(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,0,  0,0,0,0/
      data (ispn(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (ival(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,1,  0,0,0,0/
      data (ispn(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (ival(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  0,0,0,0/
      data (ispn(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (ival(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  1,0,0,0/
      data (ispn(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (ival(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  2,0,0,0/
      data (ispn(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (ival(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  3,0,0,0/
      data (ispn(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (ival(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,0/
      data (ispn(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (ival(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,1/
      data (ispn(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (ival(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,0,0,2/
      data (ispn(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (ival(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,1,0,2/
      data (ispn(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,1,0,0/

      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (ival(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,2,0,2/
      data (ispn(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,2,0,0/

      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (ival(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,2,  4,1,0,2/
      data (ispn(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,0,  0,0,0,0/

      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (ival(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,0,0,2,  4,1,0,2/
      data (ispn(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1.5,0,0,0,  0,0,0,0/

      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (ival(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,2,  4,1,0,2/
      data (ispn(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,0,  0,0,0,0/

      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (ival(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,0,0,2,  4,0,0,2/
      data (ispn(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,0,0,0,  0,0,0,0/

      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/
      data (ival(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,1,0,2,  4,0,0,2/
      data (ispn(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,1,0,0,  0,0,0,0/

      data (iocc(96,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,2,2,2,  4,0,0,2/
      data (ival(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,2,0,2,  4,0,0,2/
      data (ispn(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,2,0,0,  0,0,0,0/

      data (iocc(97,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,3,2,2,  4,0,0,2/
      data (ival(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,3,0,2,  4,0,0,2/
      data (ispn(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,3,0,0,  0,0,0,0/

      data (iocc(98,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,4,2,2,  4,0,0,2/
      data (ival(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,4,0,2,  4,0,0,2/
      data (ispn(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1,4,0,0,  0,0,0,0/

      data (iocc(99,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,5,2,2,  4,0,0,2/
      data (ival(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,5,0,2,  4,0,0,2/
      data (ispn(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,4,0,0,  0,0,0,0/

      data (iocc(100,i),i=1,29) /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,6,2,2,  4,0,0,2/
      data (ival(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,6,0,2,  4,0,0,2/
      data (ispn(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,3,0,0,  0,0,0,0/

      if (iz .lt. 1  .or.  iz .gt. 99)  then
    8    format(' Atomic number ', i5, ' not available.')
         write(slog,8)  iz
         call wlog(slog)
         call par_stop('GETORB-0')
      endif

      ion = nint(xion)
      delion=xion-ion
      index = iz - ion
      ilast = 0
      iscr = 0
      iion = 0
      iholep = ihole

c     find last occupied orbital (ilast) and iion for delion.ge.0
      do 30 i=29,1,-1
         if (iion.eq.0 .and. iocc(index,i).gt.delion) iion=i
         if (ilast.eq.0 .and. iocc(index,i).gt.0) ilast=i
 30   continue

      if (ihole.gt.0) then
         if ( iocc(index,ihole) .lt. 1 ) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
         endif
      endif
      if (ihole.eq.ilast) then 
         if ( iocc(index,ihole)-delion.lt.1) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
        endif
      endif

c        the recipe for final state atomic configuration is changed
c        from iz+1 prescription, since sometimes it changed occupation
c        numbers in more than two orbitals. This could be consistent
c        only with s02=0.0. New recipe remedy this deficiency.

c     find where to put screening electron
      index1 = index + 1
      do 10  i = 1, 29
 10   if (iscr.eq.0 .and. (iocc(index1,i)-iocc(index,i)).gt.0.5) iscr=i
c     special case of hydrogen like ion
c     if (index.eq.1) iscr=2

c     find where to add or subtract charge delion (iion).
c     if (delion .ge. 0) then
c        removal of electron charge
c        iion is already found
      if (delion .lt. 0) then
c        addition of electron charge
         iion = iscr
c        except special cases
         if (ihole.ne.0 .and.
     1       iocc(index,iscr)+1-delion.gt.2*abs(kappa(iscr))) then
             iion = ilast
             if (ilast.eq.iscr .or. iocc(index,ilast)-delion.gt.
     1                          2*abs(kappa(ilast)) ) iion = ilast + 1
         endif
      endif

      norb = 0
      do 19 i=-4, 3
 19   iorb(i) = 0
      do 20  i = 1, 29
         if (iocc(index,i).gt.0 .or. (i.eq.iscr .and. ihole.gt.0)
     1       .or. (i.eq.iion .and. iocc(index,i)-delion.gt.0))  then
            if (i.ne.ihole .or. iocc(index,i).ge.1) then
               norb = norb + 1
               nqn(norb) = nnum(i)
               nk(norb)  = kappa(i)
               xnel(norb) = iocc(index,i)
               if (i.eq.ihole) then
                  xnel(norb) = xnel(norb) - 1
                  iholep = norb
               endif
               if (i.eq.iscr .and. ihole.gt.0)  xnel(norb)=xnel(norb)+1
               xnval(norb)= ival(index,i)
               if ((kappa(i).eq.-4 .or. kappa(i).eq.3) .and. iunf.eq.0)
     1           xnval(norb) = 0
               xmag(norb) = ispn(index,i)
               iorb(nk(norb)) = norb
               if (i.eq.ihole .and. xnval(norb).ge.1)
     1                         xnval(norb) = xnval(norb) - 1
               if (i.eq.iscr .and. ihole.gt.0) 
     1                         xnval(norb) = xnval(norb) + 1
               if (i.eq.iion)  xnel(norb) = xnel(norb) - delion
               if (i.eq.iion)  xnval(norb) = xnval(norb) - delion
            endif
         endif
   20 continue
      norbco = norb

c     check that all occupation numbers are within limits
      do 50 i = 1, norb
         if ( xnel(i).lt.0 .or.  xnel(i).gt.2*abs(nk(i)) .or.
     1       xnval(i).lt.0 .or. xnval(i).gt.2*abs(nk(i)) ) then
            write (slog,55) i
   55       format(' error in getorb.f. Check occupation number for ',
     1      i3, '-th orbital. May be a problem with ionicity.')
            call wlog(slog)
            call par_stop('GETORB-99')
         endif
  50  continue
c      do 60 i=1,norb
c60    xnval(i) = 0.0d0
c60    xnval(i) = xnel(i)
            
      return
      end
      double precision function getxk (e)
      implicit double precision (a-h, o-z)

c     Make xk from energy(in Hartrees) as
c          k =  sqrt(2*e)  for e > 0  (above the edge)
c          k = -sqrt(-2*e)  for e < 0  (below the edge)
      getxk = sqrt(abs(2*e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
      subroutine sthead (ntitle, title, nph, iz, rmt, rnrm,
     1                  xion, ihole, ixc,
     2                  vr0, vi0, gamach, xmu, xf, vint, rs,
     2                  nohole, lreal,  rgrd)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.bin

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/vers.h
      character*12 vfeff
c                       123456789012  
      parameter (vfeff='Feff 8.50   ')
c= ../HEADERS/vers.h}

      dimension xion(0:nphx)
      dimension iz(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(nheadx), store
      character*16 s1, s2

      character*10 shole(0:29)
      character*8  sout(0:7)
      data shole /'no hole',   'K  shell',  'L1 shell',  'L2 shell',
     2            'L3 shell',  'M1 shell',  'M2 shell',  'M3 shell',
     3            'M4 shell',  'M5 shell',  'N1 shell',  'N2 shell',
     4            'N3 shell',  'N4 shell',  'N5 shell',  'N6 shell',
     5            'N7 shell',  'O1 shell',  'O2 shell',  'O3 shell',
     6            'O4 shell',  'O5 shell',  'O6 shell',  'O7 shell',
     7            'P1 shell',  'P2 shell',  'P3 shell',  'P4 shell',
     8            'P5 shell',  'R1 shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 'DH - HL ',
     1           'DH + HL ', 'val=s+d ', 'sigmd(r)', 'sigmd=c '/


c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      if (ntitle .ge. 1 ) then
         ii = istrln(title(1)) 
         if (ii.gt.1)  then
            write(store,100)  title(1)(1:), vfeff
         else
            write(store,102)  vfeff
         endif
      else
         write(store,102)   vfeff
      endif
  100 format( a55, t66, a12)
  102 format( t66, a12)
      title(1) = store
      nstor = 1

c     remove empty title lines
      do 120  ititle = 2, ntitle
         ii = istrln ( title (ititle) ) 
         if (ii.le.1)  goto 120
         nstor = nstor+1
         title(nstor) = title (ititle)
  120 continue
      ntitle = nstor

c     add more title lines
      if (xion(0) .ne. 0)  then
         ntitle = ntitle + 1
         write(title(ntitle),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, xion(0), shole(ihole)
      else
         ntitle = ntitle + 1
         write(title(ntitle),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2, ' Rmt=',f6.3, ' Rnm=',f6.3,
     1       ' Ion=',f5.2,  1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)
c     if (nohole.ge.0)  then
c        ntitle = ntitle + 1
c        write(title(ntitle),142)
c 142    format ('Calculations done with no core hole.')
c     endif
      if (lreal.ge.1 .or. (abs(rgrd - 0.05) .gt. 1.0e-5)) then
        ntitle = ntitle + 1
        s1 = ' '
        if (lreal.gt.1)  then
c        write(title(ntitle),144)
c 144    format ('Calculations done using only real phase shifts.')
         s1 = 'RPHASES'
        elseif (lreal.eq.1) then
c        ntitle = ntitle + 1
c        write(title(ntitle),145)
c 145    format ('Calculations done using only real self energy.')
         s1 = 'RSIGMA'
        endif
        s2 = '  '
        if (abs(rgrd - 0.05) .gt. 1.0e-5)  then
         write(s2,146)  rgrd
  146    format ('  RGRID', f7.4)
        endif
        ilen = istrln(s1)
        title(ntitle) = s1(1:ilen) // s2
      endif

      do 150  iph = 1, nph
         if (xion(iph) .ne. 0)  then
            ntitle = ntitle + 1
            write(title(ntitle),160)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, xion(iph)
         else
            ntitle = ntitle + 1
            write(title(ntitle),170)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',f5.2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         ntitle = ntitle + 1
         write(title(ntitle),180)  gamach*hart, sout(ixc), vi0*hart,
     1                           vr0*hart
      else
         ntitle = ntitle + 1
         write(title(ntitle),190)  gamach*hart, sout(ixc)
      endif
      ntitle = ntitle + 1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(title(ntitle),200)  xmu*hart, xf/bohr, vint*hart, rs

      return
      end

      subroutine wthead (io, ntitle, title)
c     Dump title lines to unit io, which must be open. 
      integer io, i, ll
      character*80 title(ntitle)

c     nice for UNIX to use with gnuplot etc.,
      do 310 i = 1, ntitle
         ll = istrln(title(i))
         write(io,300)  title(i)(1:ll)
  300    format (a)
  310 continue

      return
      end
      function itoken (word,flname)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w
      character*20 flname
      integer itoken

      w = word(1:4)
      call upper(w)
      
c     Tokens for feff.inp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
         if     (w .eq. 'ATOM')  then
            itoken = 1
         elseif (w .eq. 'HOLE')  then
            itoken = 2
         elseif (w .eq. 'OVER')  then
            itoken = 3
         elseif (w .eq. 'CONT')  then
            itoken = 4
         elseif (w .eq. 'EXCH')  then
            itoken = 5
         elseif (w .eq. 'ION ')  then
            itoken = 6
         elseif (w .eq. 'TITL')  then
            itoken = 7
         elseif (w .eq. 'FOLP')  then
            itoken = 8
         elseif (w .eq. 'RPAT' .or. w .eq. 'RMAX')  then
            itoken = 9
         elseif (w .eq. 'DEBY')  then
            itoken = 10
         elseif (w .eq. 'RMUL')  then
            itoken = 11
         elseif (w .eq. 'SS  ')  then
            itoken = 12
         elseif (w .eq. 'PRIN')  then
            itoken = 13
         elseif (w .eq. 'POTE')  then
            itoken = 14
         elseif (w .eq. 'NLEG')  then
            itoken = 15
         elseif (w .eq. 'CRIT')  then
            itoken = 16
         elseif (w .eq. 'NOGE')  then
            itoken = 17
         elseif (w .eq. 'IORD')  then
            itoken = 18
         elseif (w .eq. 'PCRI')  then
            itoken = 19
         elseif (w .eq. 'SIG2')  then
            itoken = 20
         elseif (w .eq. 'XANE')  then
            itoken = 21
         elseif (w .eq. 'CORR')  then
            itoken = 22
         elseif (w .eq. 'AFOL')  then
            itoken = 23
         elseif (w .eq. 'EXAF')  then
            itoken = 24
         elseif (w .eq. 'POLA')  then
            itoken = 25
         elseif (w .eq. 'ELLI')  then
            itoken = 26
         elseif (w .eq. 'RGRI')  then
            itoken = 27
         elseif (w .eq. 'RPHA')  then
            itoken = 28
         elseif (w .eq. 'NSTA')  then
            itoken = 29
         elseif (w .eq. 'NOHO')  then
            itoken = 30
         elseif (w .eq. 'SIG3')  then
            itoken = 31
         elseif (w .eq. 'JUMP')  then
            itoken = 32
         elseif (w .eq. 'MBCO')  then
            itoken = 33
         elseif (w .eq. 'SPIN')  then
            itoken = 34
         elseif (w .eq. 'EDGE')  then
            itoken = 35
         elseif (w .eq. 'SCF ')  then
            itoken = 36
         elseif (w .eq. 'FMS ')  then
            itoken = 37
         elseif (w .eq. 'LDOS')  then
            itoken = 38
         elseif (w .eq. 'INTE')  then
            itoken = 39
         elseif (w .eq. 'CFAV')  then
            itoken = 40
         elseif (w .eq. 'S02 ')  then
            itoken = 41
         elseif (w .eq. 'XES ')  then
            itoken = 42
         elseif (w .eq. 'DANE')  then
            itoken = 43
         elseif (w .eq. 'FPRI')  then
            itoken = 44
         elseif (w .eq. 'RSIG')  then
            itoken = 45
         elseif (w .eq. 'XNCD')  then
            itoken = 46
         elseif (w .eq. 'XMCD')  then
            itoken = 46
         elseif (w .eq. 'MULT')  then
            itoken = 47
         elseif (w .eq. 'UNFR')  then
            itoken = 48
         elseif (w .eq. 'TDLD')  then
            itoken = 49
         elseif (w .eq. 'PMBS')  then
            itoken = 50
         elseif (w .eq. 'PLAS')  then
            itoken = 51
         elseif (w .eq. 'SO2C')  then
            itoken = 52
         elseif (w .eq. 'SELF')  then
            itoken = 53
         elseif (w .eq. 'SFSE')  then
            itoken = 54
         elseif (w .eq. 'RCONV') then
            itoken = 55
         elseif (w .eq. 'ELNE') then !KJ new card for EELS 1-06
            itoken = 56
         elseif (w .eq. 'EXEL') then !KJ new card for EELS 1-06
            itoken = 57
         elseif (w .eq. 'MAGI') then !KJ new card for EELS 1-06
            itoken = 58
         elseif (w .eq. 'ABSO') then !KJ new card 3-06
            itoken = 59    
         elseif (w .eq. 'EGRI')  then !Josh Kas
            itoken = 60
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
      elseif (flname(1:10).eq.'spring.inp') then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     These tokens are for spring.inp (input for eq of motion method)
         if (w .eq. 'STRE')  then
            itoken = 1
         elseif (w .eq. 'ANGL')  then
            itoken = 2
         elseif (w .eq. 'VDOS')  then
            itoken = 3
         elseif (w .eq. 'PRDOS') then
            itoken = 4
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif
      return
      end


c====================================================================
      integer function nxtunt(iunit)

c  this function returns the value of the next unopened unit number
c  equal to or larger than iunit.  it will return neither unit numbers
c  0, 5, or 6 nor a negative unit number
c $Id: nxtunt.f,v 1.1.1.1 2006/01/12 06:37:42 hebhop Exp $
c $Log: nxtunt.f,v $
c Revision 1.1.1.1  2006/01/12 06:37:42  hebhop
c New version of feff. feff8.5 (Extension of feff8.4)
c Includes:
c 	1) All feff8.4 capabilities.
c 	2) Screened core hole (calculation of W).
c 	3) Multiple pole self energy calculation.
c 	4) Convolution with spectral function.
c New cards and options:
c 	1) NOHOLE 2      (screened hole)
c 	2) PLASMON ipl   (multiple pole self energy)
c 	3) SO2CONV       (convolve output with spectral function)
c 	4) SELF          (print on shell self energy as calculated by Luke)
c 	5) SFSE k0        (print off shell self energy Sigma(k0,e) )
c
c Revision 1.1.1.1  2000/02/11 02:23:58  alex
c Initialize feff82
c
c Revision 1.10  1999/04/02 21:32:47  newville
c cleaned up nxtunt (matt)
c
c Revision 1.9  1999/02/11 20:08:08  alex
c x39 version: dim.h + misc. small changes
c
c Revision 1.8  1998/12/29 23:59:07  alex
c feff8x35 version
c
c Revision 1.7  1998/11/19 03:23:11  alex
c feff8x32 version
c
c Revision 1.6  1998/10/26 14:11:16  ravel
c no comments beyond column 71
c
c Revision 1.5  1998/10/18 21:47:51  alex
c feff8x30 version implements Broyden algorithm for self-consistency
c
c Revision 1.4  1998/02/24 18:31:37  ravel
c I should really be more careful.  This is the last commitment done
c      cright.
c
c Revision 1.1.1.1  1997/04/27 20:18:03  ravel
c Initial import of xanes sources, version 0.37
c
c Revision 1.1  1996/06/23 16:05:02  bruce
c Initial revision
c

       integer iunit
       logical open

       nxtunt = max(1, iunit) - 1
 10    continue
       nxtunt = nxtunt + 1
       if ((nxtunt.eq.5).or.(nxtunt.eq.6)) nxtunt = 7
       inquire (unit=nxtunt, opened=open)
       if (open) go to 10
       return
c  end integer function nxtunt
       end

c====================================================================
c     Periodic table of the elements
c     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

c     atwts(iz)  single precision fn, returns atomic weight
c     atwtd(iz)  double precision fn, returns atomic weight
c     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real function atwts (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      character*2 sym
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
c     PeRiodic TaBle Block Data

c     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      common /atwtco/ weight(103)

      character*2 sym
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
      subroutine pijump (ph, old)
      implicit double precision (a-h, o-z)

c     removes jumps of 2*pi in phases

c     ph = current value of phase (may be modified on output, but
c          only by multiples of 2*pi)
c     old = previous value of phase

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
      parameter (twopi = 2 * pi)
      dimension xph(3)

      xph(1) = ph - old
      jump =  (abs(xph(1))+ pi) / twopi
      xph(2) = xph(1) - jump*twopi
      xph(3) = xph(1) + jump*twopi


      xphmin = min (abs(xph(1)), abs(xph(2)), abs(xph(3)))
      isave = 0
      do 10  i = 1, 3
         if (abs (xphmin - abs(xph(i))) .le. 0.01)  isave = i
   10 continue
      if (isave .eq. 0)  then
         call par_stop('pijump')
      endif

      ph = old + xph(isave)

      return
      end
      subroutine rdhead (io, nhead, head, lhead)
      implicit double precision (a-h, o-z)

c     Reads title line(s) from unit io.  Returns number of lines
c     read.  If more than nheadx lines, skips over them.  End-of-header
c     marker is a line of 1 blank, 71 '-'s.
c     lhead is length of each line w/o trailing blanks.
c     header lines returned will have 1st space on line blank for
c     carriage control

      character*80 head(nhead)
      dimension lhead(nhead)
      character*80  line

      n = 0
      nheadx = nhead
      nhead = 0
   10 read(io,20)  line
   20    format(a)
         if (line(4:11) .eq. '--------')  goto 100
         n = n+1
         if (n .le. nheadx)  then
            head(n) = line
            lhead(n) = istrln(head(n))
            nhead = n
         endif
      goto 10
  100 continue
      return
      end
      subroutine rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel, 
     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c  opens pot.bin file and reads following information
c  General:
c     ntitle - number of title lines
c     title  - title itself
c     emu    - edge position (x-ray energy for final state at Fermi level)
c  Muffin-tin geometry
c     rmt    - muffin-tin radii
c     imt    - index of radial grid just below muffin-tin radii
c     rnrm   - Norman radii
c     inrm   - index of radial grid just below Norman radii
c     rnrmav - average Norman radius
c     folp   - overlap parameter for rmt
c     folpx  - maximum value for folp
c     xnatph - number of atoms of each potential type
c  Atomic orbitals info (Dirac spinors)
c     dgc0   - upper component for initial orbital
c     dpc0   - lower component for initial orbital
c     dgc    - upper components for all atomic orbitals
c     dpc    - lower components for all atomic orbitals
c     adgc   - development coefficient for upper components
c     adpc   - development coefficient for lower components
c     xnval  - number of valence electrons for each atomic orbital
c              used for core-valence separation and non-local exchange
c     eorb  - atomic enrgy of each orbital for the absorber
c  Electron density information 
c     rhoint - interstitial density
c     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
c     xf     - estimate of momentum at Fermi level from rhoint
c     edens  - total electron density
c     edenvl - density from valence electrons
c     dmag   - density for spin-up minus density for spin-down
c     qtotel - total charge of a cluster
c     qnrm   - charge accumulated inside Norman sphere as result of SCF
c     xnmues - occupation numbers of valence orbitals from SCF procedure
c  Potential information
c     xmu    - Fermi level position
c     ecv    - core-valence separation energy
c     vint   - muffin-tin zero energy (interstitial potential)
c     vclap  - Coulomb potential
c     vtot   - vclap + xc potential from edens
c     vvalgs - vclap + xc potential from edenvl (EXCHANGE 5 model)
c  Specific data for convolution with excitation spectrum (see mbconv)
c     s02    - many body reduction factor S_0^2 
c     erelax - estimate of relaxation energy = efrozen - emu, where
c              efrozen is edge position estimate from Koopmans theorem
c     wp     - estimate of plasmon frequency from rhoint

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)

      character*80 title(nheadx)

      dimension dum(13)

  10  format(a)
   20 format (bn, i15)

      open (unit=3, file='pot.bin', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,
     1            jumprm, iunf
  30  format(9(1x,i4))
c     nph and npadx are not passed to calling subroutine
      do 133  i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
  133 continue
c     Misc double precision stuff from pot.bin
      call rdpadd(3, npadx, dum(1), 13)
      rnrmav = dum(1)
      xmu    = dum(2)
      vint   = dum(3)
      rhoint = dum(4)
      emu    = dum(5)
      s02    = dum(6)
      erelax = dum(7)
      wp     = dum(8)
      ecv    = dum(9)
      rs     = dum(10)
      xf     = dum(11)
      qtotel = dum(12)
      totvol = dum(13)

c     read imt
      read (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      call rdpadd(3, npadx, rmt(0), nph+1)
c     read inrm
      read (3, 40) (inrm(i),i=0,nph)
      read (3, 40) (iz(i),i=0,nph)
      read (3, 40) (kappa(i),i=1,30)
      call rdpadd(3, npadx, rnrm(0), nph+1)
      call rdpadd(3, npadx, folp(0), nph+1)
      call rdpadd(3, npadx, folpx(0), nph+1)
      call rdpadd(3, npadx, xnatph(0), nph+1)
      call rdpadd(3, npadx, xion(0), nph+1)
      call rdpadd(3, npadx, dgc0(1), 251)
      call rdpadd(3, npadx, dpc0(1), 251)
      call rdpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call rdpadd(3, npadx, eorb(1), 30)
      do 50 iph=0,nph
 50   read (3, 60) (iorb(i,iph),i=-4,3)
 60   format(8(1x,i2))
      call rdpadd(3, npadx, qnrm(0), nph+1 )
      nn = (lx+1)*(nph+1)
      call rdpadd(3, npadx, xnmues(0,0), nn )
      close (unit=3)

      return
      end
      subroutine rdxsph ( ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1               ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phase.bin' 
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
c     kinit, linit, ilinit,  - initial state kappa and ang. mom.
c     lmaxp1  -largest lmax in problem + 1

c     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

c     Local staff
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
c     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp 
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp 
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

c     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
c        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.
     3          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end
      subroutine setkap(ihole, kinit, linit)
      implicit double precision (a-h, o-z)

c     Set initial state ang mom and quantum number kappa
c     ihole  initial state from ihole    
c     1      K    1s      L=0 -> linit=0 
c     2      LI   2s      L=0 -> linit=0
c     3      LII  2p 1/2  L=1 -> linit=1
c     4      LIII 2p 3/2  L=1 -> linit=1
c     5+     etc.
      if (ihole.le. 2 .or. ihole.eq. 5 .or. ihole.eq.10 .or.
     1    ihole.eq.17 .or. ihole.eq.24 .or. ihole.eq.27)  then
c        hole in s state
         linit = 0
         kinit = -1
      elseif (ihole.eq. 3 .or. ihole.eq. 6 .or. ihole.eq.11 .or.
     1        ihole.eq.18 .or. ihole.eq.25 .or. ihole.eq.30)  then
c        hole in p 1/2 state
         linit = 1
         kinit = 1
      elseif (ihole.eq. 4 .or. ihole.eq. 7 .or. ihole.eq.12 .or.
     1        ihole.eq.19 .or. ihole.eq.26)  then
c        hole in p 3/2 state
         linit = 1
         kinit = -2
      elseif (ihole.eq. 8 .or. ihole.eq.13 .or.
     1        ihole.eq.20 .or. ihole.eq.27)  then
c        hole in d 3/2 state
         linit = 2
         kinit = 2
      elseif (ihole.eq. 9 .or. ihole.eq.14 .or.
     1        ihole.eq.21 .or. ihole.eq.28)  then
c        hole in d 5/2 state
         linit = 2
         kinit = -3
      elseif (ihole.eq.15 .or. ihole.eq.22)  then
c        hole in  f 5/2 state
         linit = 3
         kinit = 3
      elseif (ihole.eq.16 .or. ihole.eq.23)  then
c        hole in  f 7/2 state
         linit = 3
         kinit = -4
      else
c        some unknown hole
         call par_stop('invalid hole number in setkap')
      endif

      return
      end
C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
C SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
C SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks or tabs, or a comma and zero or more blanks.
C
C     ARGS        I/O      DESCRIPTION
C     ----        ---      -----------
C     S            I       CHAR*(*)  String to be broken up
C     NWORDS      I/O      Input:  Maximum number of words to get
C                          Output: Number of words found
C     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
C                          Contains words found.  WORDS(J), where J is
C                          greater then NWORDS found, are undefined on
C                          output.
C
C      Written by:  Steven Zabinsky, September 1984
C      Tab char added July 1994.
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
C     there is a tab character here               ^.

C  -- BETW    .TRUE. if between words
C     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

C  -- Maximum number of words allowed
      WORDSX = NWORDS

C  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

C  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

C  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END
      SUBROUTINE UNTAB (STRING)
C REPLACE TABS WITH BLANKS :    TAB IS ASCII DEPENDENT
      INTEGER        ITAB , I, ILEN, ISTRLN
      PARAMETER      (ITAB = 9)
      CHARACTER*(*)  STRING, TAB*1
      EXTERNAL ISTRLN
      TAB  = CHAR(ITAB)
      ILEN = MAX(1, ISTRLN(STRING))
 10   CONTINUE 
        I = INDEX(STRING(:ILEN), TAB ) 
        IF (I .NE. 0) THEN
            STRING(I:I) = ' '
            GO TO 10
        END IF
      RETURN
C END SUBROUTINE UNTAB
      END

      logical function iscomm (line)
c     returns true if line is a comment or blank line, false otherwise
c#mn{ rewritten to allow ";*%#" as comment characters
       character*(*) line
       iscomm = ((line.eq.' ').or.(index(';*%#',line(1:1)).ne.0))
c#mn}
      return
      end
      subroutine str2dp(str,dpval,ierr)
c  return dp number "dpval" from character string "str"
c  if str cannot be a number, ierr < 0 is returned.
      character*(*) str, fmt*15 
      double precision dpval
      integer  ierr , lenmax
      parameter ( lenmax = 40)
      logical  isnum
      external isnum
      ierr = -99
      if (isnum(str)) then
         ierr = 0
         write(fmt, 10) min(lenmax, len(str))
 10      format('(bn,f',i3,'.0)')
         read(str, fmt, err = 20, iostat=ierr) dpval
      end if    
      if (ierr.gt.0) ierr = -ierr
      return
 20   continue
      ierr = -98
      return
c end subroutine str2dp
      end
      subroutine str2re(str,val,ierr)
c  return real from character string "str"
      character*(*) str 
      double precision dpval
      real     val
      integer  ierr
      call str2dp(str,dpval,ierr)
      if (ierr.eq.0) val = dpval
      return
c end subroutine str2re
      end
      subroutine str2in(str,intg,ierr)
c  return integer from character string "str"
c  returns ierr = 1 if value was clearly non-integer
      character*(*) str 
      double precision val, tenth
      parameter (tenth = 1.d-1)
      integer  ierr, intg
      call str2dp(str,val,ierr)
      if (ierr.eq.0) then
         intg = int(val)
         if ((abs(intg - val) .gt. tenth))  ierr = 1
       end if
      return
c end subroutine str2in
      end
       logical function isnum (string)
c  tests whether a string can be a number. not foolproof!
c  to return true, string must contain:
c    - only characters in  'deDE.+-, 1234567890' (case is checked)
c    - no more than one 'd' or 'e' 
c    - no more than one '.'
c  matt newville
       character*(*)  string,  number*20
c note:  layout and case of *number* is important: do not change!
       parameter (number = 'deDE.,+- 1234567890')
       integer   iexp, idec, i, j, istrln
       external  istrln
       iexp  = 0
       idec  = 0
       isnum = .false. 
       do 100  i = 1, max(1, istrln(string))
          j = index(number,string(i:i))
          if (j.le.0)               go to 200
          if((j.ge.1).and.(j.le.4)) iexp = iexp + 1
          if (j.eq.5)               idec = idec + 1
 100   continue
c  every character in "string" is also in "number".  so, if there are 
c  not more than one exponential and decimal markers, it's a number
       if ((iexp.le.1).and.(idec.le.1)) isnum = .true.
 200   continue
       return
c  end logical function isnum
       end
      subroutine wlog (string)
      character*(*) string

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c     This output routine is used to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.

c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.

c     make sure not to write trailing blanks

   10 format (a)

c     Suppress output in sequential loops
      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print10
         if (par_type .ne. 3) write(11,10)
      else
         print10, string(1:il)
         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
c     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
c     r grid point at index j
      rr = exp (xx(j))
      return
      end

      function ii(r)
      implicit double precision (a-h, o-z)
c     index of grid point immediately below postion r
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
c
c PAD library:   Packed Ascii Data 
c   these routines contain code for handling packed-ascii-data  
c   (pad) arrays for writing printable character strings that 
c   represent real or complex scalars and arrays to a file.
c
c routines included in padlib are (dp==double precision):
c   wrpadd     write a dp array as pad character strings
c   wrpadx     write a dp complex array as pad character strings
c   rdpadr     read a pad character array as a real array
c   rdpadd     read a pad character array as a dp  array
c   rdpadc     read a pad character array as a complex array
c   rdpadx     read a pad character array as a dp complex array
c   pad        internal routine to convert dp number to pad string
c   unpad      internal routine to pad string to dp number
c
c routines not included, but required by padlib:
c     triml, istrln, wlog
c
c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2001 Matthew Newville, The University of Chicago
c Copyright (c) 1992--1996 Matthew Newville, University of Washington
c
c Permission to use and redistribute the source code or binary forms of
c this software and its documentation, with or without modification is
c hereby granted provided that the above notice of copyright, these
c terms of use, and the disclaimer of warranty below appear in the
c source code and documentation, and that none of the names of The
c University of Chicago, The University of Washington, or the authors
c appear in advertising or endorsement of works derived from this
c software without specific prior written permission from all parties.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
c//////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
       subroutine wrpadd(iout,npack,array,npts)
c
c write a dp array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       double precision array(*), xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = array(i)
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadx(iout,npack,array,npts)
c write complex*16 array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex*16 array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = dimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadr(iout,npack,array,npts)
c
c write a real array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       real    array(*)
       double precision xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = dble(array(i))
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadc(iout,npack,array,npts)
c write complex (*8) array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex    array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = aimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine rdpadd(iou,npack,array,npts)
c read dparray from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       double precision    array(*), unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = tmp
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadr(iou,npack,array,npts)
c read real array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       real    array(*)
       double precision unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = real(tmp)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadc(iou,npack,array,npts)
c read complex array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
       subroutine rdpadx(iou,npack,array,npts)
c read complex*16 array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex*16  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end

c --padlib--
       subroutine pad(xreal,npack,str)
c  convert dp number *xreal* to packed-ascii-data string *str*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer  iexp, itmp, isgn, i, npack, j
       double precision xreal, xwork, xsave,onem, tenth
       parameter (onem  =  0.99999999997d0)
       parameter (tenth =  0.099999999994d0)
       character str*(*)
c
       str      = ' '
       xsave    = min(huge, max(-huge, xreal))
       isgn     = 1
       if (xsave.le.0) isgn = 0
c
       xwork    = dabs( xsave )
       iexp     = 0
       if ((xwork.lt.huge).and.(xwork.gt.tiny))  then
          iexp  =   1 + int(log(xwork) / tenlog  )
       else if (xwork.ge.huge) then
          iexp  = ihuge
          xwork = one
       else if (xwork.le.tiny)  then
          xwork = zero
       end if
c force xwork between ~0.1 and ~1
c note: this causes a loss of precision, but 
c allows backward compatibility
       xwork    = xwork / (ten ** iexp)
 20    continue
       if (xwork.ge.one) then
          xwork = xwork * 0.100000000000000d0
          iexp  = iexp + 1
       else if (xwork.le.tenth) then
          xwork = xwork * ten
          iexp  = iexp - 1
       endif
       if (xwork.ge.one) go to 20

       itmp     = int ( ibas2 * xwork ) 
       str(1:1) = char(iexp  + ioff + ibas2 )
       str(2:2) = char( 2 * itmp + isgn + ioff)
       xwork    = xwork * ibas2 - itmp
       if (npack.gt.2) then
          do 100 i = 3, npack
             itmp     = int( base * xwork + 1.d-9)
             str(i:i) = char(itmp + ioff)
             xwork    = xwork * base - itmp
 100      continue
       end if
       if (xwork.ge.0.5d0) then
          i = itmp + ioff + 1
          if (i.le.126) then
             str(npack:npack)= char(i)
          else 
             j = ichar(str(npack-1:npack-1))
             if (j.lt.126) then
                str(npack-1:npack-1) = char(j+1)
                str(npack:npack)     = char(37)
             endif 
          endif
       endif
       return
       end
c --padlib--
       double precision function unpad(str,npack)
c
c  convert packed-ascii-data string *str* to dp number *unpad*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       double precision sum
       integer   iexp, itmp, isgn, i, npack
       character str*(*)
       unpad = zero
       if (npack.le.2) return
       iexp  =     (ichar(str(1:1)) - ioff   ) - ibas2
       isgn  = mod (ichar(str(2:2)) - ioff, 2) * 2 - 1
       itmp  =     (ichar(str(2:2)) - ioff   ) / 2
       sum   = dble(itmp/(base*base))
c       do 100 i = 3, npack
c          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
c 100   continue
       do 100 i = npack, 3, -1
          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
 100   continue
       unpad = 2 * isgn * ibase * sum * (ten ** iexp)
cc       print*, sum, iexp,unpad
       return
       end
c --padlib--
c end of pad library
c ----------
       integer function iread(lun,string)
c
c generalized internal read:
c    read a string the next line of an opened file 
c    unit, returning the real length of string
c 
c inputs:   
c   lun     opened file unit number
c outputs:
c   string  string read from file
c returns:
c   iread   useful length of string, as found from 
c                  sending string to 'sclean' to 
c                  remove non-printable characters
c                   and then istrln  
c           or
c              -1   on 'end-of-file'
c              -2   on 'error'
c
c copyright (c) 1999  Matthew Newville
       implicit none
       character*(*) string
       integer    lun, istrln
       external   istrln
       string = ' '
 10    format(a)
       read(lun, 10, end = 40, err = 50) string
       call sclean(string)
       iread = istrln(string)
       return
 40    continue 
       string = ' '
       iread = -1
       return
 50    continue 
       string = ' '
       iread = -2
       return
       end
       subroutine sclean(str) 
c
c  clean a string, especially for strings passed between 
c  different file systems, or from C functions:
c
c   1. characters in the range char(0), or char(10)...char(15) 
c      are interpreted as end-of-line characters, so that all
c      remaining characters are explicitly blanked.
c   2. all other characters below char(31) (including tab) are
c      replaced by a single blank
c
c  this is mostly useful when getting a string generated by a C 
c  function and for handling dos/unix/max line-endings.
c
c copyright (c) 1999  Matthew Newville
       character*(*) str, blank*1
       parameter (blank = ' ')
       integer i,j,is
       do 20 i = 1, len(str)
          is = ichar(str(i:i))
          if ((is.eq.0) .or. ((is.ge.10) .and. (is.le.15))) then
             do 10 j= i, len(str)
                str(j:j) = blank
 10          continue
             return
          elseif (is.le.31)  then
             str(i:i)  = blank
          end if
 20    continue 
       return
c end subroutine sclean
       end

      SUBROUTINE rdcmt(iUnt,Cmt)
      INTEGER iUnt, i1
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER TmpCmt(4), ch
      LOGICAL CmtLin

      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,*,END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) goto 5
      END DO
      
      BACKSPACE(iUnt)
      
 10   CONTINUE
      
      RETURN
      END
      subroutine setgam (iz, ihole, gamach)

c     Sets gamach, core hole lifetime.  Data comes from graphs in
c     K. Rahkonen and K. Krause,
c     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.
c     output gamach is in eV

      implicit double precision (a-h, o-z)

      dimension zh(8,16), gamh(8,16)

      dimension zk(8), gamkp(8)
      parameter (ryd  = 13.605 698d0)
      parameter (hart = 2*ryd)
      character*512 slog


c     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
c     trouble.
c     Gam arrays contain the gamma values.
c     We will take log10 of the gamma values so we can do linear
c     interpolation from a log plot.

      data  zh   / 0.99,  10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1,
     2              0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0,  95.1,
     3              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     4              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     5              0.99,  20.0, 28.0, 30.0, 36.0, 53.0,  80.0, 95.1,
     6              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     7              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     8              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     9              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     *              0.99,  30.0, 40.0, 47.0, 50.0, 63.0,  80.0, 95.1,
     1              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     2              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     3              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     4              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     5              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0,
     6              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0/

      data  gamh / 0.02,  0.28, 0.75,  4.8, 10.5, 21.0, 60.0, 105.0,
     2              0.07,  3.9,  3.8,  7.0,  6.0,  3.7,  8.0,  19.0,
     3              0.001, 0.12,  1.4,  0.8,  2.6,  4.1,   6.3, 10.5,
     4              0.001, 0.12, 0.55,  0.7,  2.1,  3.5,   5.4,  9.0,
     5              0.001,  1.0,  2.9,  2.2,  5.5, 10.0,  22.0, 22.0,
     6              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  15.0, 16.0,
     7              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  10.0, 10.0,
     8              0.0006,0.09, 0.07, 0.48,  1.0,  4.0,   2.7,  4.7,
     9              0.0006,0.09, 0.07, 0.48, 0.87,  2.2,   2.5,  4.3,
     *              0.001,0.001,  6.2,  7.0,  3.2, 12.0,  16.0, 13.0,
     1              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     2              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     3              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     4              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     5              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9,
     6              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9/

c     Since feff8 can be called any number of times . ALA

      if (ihole .le. 0)  then
         gamach = 0
         write(slog,'(a,1pe13.5)') ' No hole in SETGAM, gamach = ', 
     1                             gamach
         call wlog(slog)
         return
      endif
      if (ihole .gt. 16)  then
         call wlog(' This version of FEFF will set gamach = 0.1 eV ' //
     1             ' for O1 and higher hole')
         call wlog(' You can use CORRECTIONS card  to set ' //
     1   ' gamach = 0.1 + 2*vicorr ')
c        stop 'SETGAM-2'
      endif

      zz = iz
      if (ihole .le. 16)  then
         do 10  i = 1, 8
            gamkp(i) = log10 (gamh(i,ihole))
            zk(i) = zh(i,ihole)
   10    continue
         call terp (zk, gamkp, 8, 1, zz, gamach)
      else
c     include data from the tables later.
c     Now gamach=0.1eV for any O-hole for any element.
         gamach = -1.0
      endif

c     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach


      return
      end
      subroutine iniptz(ptz,iptz,modus)
        !KJ This routine rewrites the ptz-matrix.

      implicit none
c     which polarization tensor to create      
      integer iptz
c     two ways of working
      integer modus
c     the polarization tensor
      complex*16 ptz(-1:1,-1:1)

      complex*16 zero,one,coni
      parameter (zero=(0,0),one=(1,0),coni=(0,1))
      integer i,j
      real*8 unity(3,3)



      if (iptz.lt.1.or.iptz.gt.10) then
          call wlog('Inieln sees weird iptz - returns without 
     1      changing ptz - danger of calculating nonsense !!')
      endif


      do i=1,3
      do j=1,3
	unity(i,j)=dble(0)
      enddo
      unity(i,i)=dble(1)/dble(3)
      enddo
      do i=-1,1
      do j=-1,1
        ptz(i,j)=zero
      enddo
      enddo
      

      if (modus.eq.1) then
! work in spherical coordinates

         if(iptz.eq.10) then
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
         else
            i=(iptz-1)/3+1  ! row index
            j=iptz-3*(i-1)  ! column index
            i=i-2 !shift from 1..3 to -1..1
            j=j-2
            ptz(i,j)=one
         endif  


      elseif(modus.eq.2) then
! work in carthesian coordinates      

      if (iptz.eq.10) then ! orientation averaged spectrum
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
	
        elseif (iptz.eq.1) then   ! x x*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=-one/dble(2)
        ptz(1,-1)=-one/dble(2)
        elseif (iptz.eq.5) then ! y y*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=one/dble(2)
        ptz(1,-1)=one/dble(2)
         elseif (iptz.eq.9) then ! z z*
        ptz(0,0)=one
        elseif (iptz.eq.2) then ! x y*
          ptz(1,1)=one*coni/dble(2)
        ptz(-1,-1)=-one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.4) then ! x* y
          ptz(1,1)=-one*coni/dble(2)
        ptz(-1,-1)=one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.3) then ! x z*
          ptz(-1,0)=one/dsqrt(dble(2))
        ptz(1,0)=-one/dsqrt(dble(2))
        elseif (iptz.eq.7) then ! x* z
          ptz(0,-1)=one/dsqrt(dble(2))
        ptz(0,1)=-one/dsqrt(dble(2))
        elseif (iptz.eq.6) then ! y z*
          ptz(-1,0)=-one*coni/dsqrt(dble(2))
        ptz(1,0)=-one*coni/dsqrt(dble(2))
        elseif (iptz.eq.8) then ! y* z
          ptz(0,-1)=one*coni/dsqrt(dble(2))
        ptz(0,1)=one*coni/dsqrt(dble(2))
      endif
      
      
        else
          stop 'alien modus in inieln'
        endif


      return
      end


C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.


      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      DOUBLE PRECISION X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      DOUBLE PRECISION A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) THEN
         RETURN
      END IF
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
c///////////////////////////////////////////////////////////////////////
c PAR Subroutines
c Written by J. Sims, NIST, 2001

c This software was developed at the National Institute of Standards
c and Technology by employees of the Federal Government in the course
c of their official duties. Pursuant to title 17 Section 105 of
c the United States Code this software is not subject to copyright
c protection and is in the public domain. PAR is an experimental
c system. NIST assumes no responsibility whatsoever for its use by
c other parties, and makes no guarantees, expressed or implied, about 
c its quality, reliability, or any other characteristic. We would
c appreciate acknowledgement if the software is used.

c This software can be redistributed and/or modified freely provided
c that any derivative works bear some notice that they are derived from
c it, and any modified versions bear some notice that they have been
c modified.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c  **************************************************
c  Parallel feff8 routines
c  Jim Sims
c  **************************************************

      subroutine par_begin
c  **************************************************
c  Initializations for parallel version(s)
c  **************************************************

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}

c-- So cvd or dbx can attach to a running process
c     call sleep(30) 

      call MPI_INIT(ierrorflag)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierrorflag)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierrorflag)
      this_process = my_rank

      par_type = 0
      parallel_run = .true.
c-- The following variable will be used for IO that should only be
c-- done in one process.
      master = (my_rank .eq. 0)

      worker = (.not. master)
      if (worker) par_type = 1

c     write(6,*) 'this process = ',this_process, ' worker = ',worker

      if (master) write(6,*) 'Number of processors = ',numprocs

      return
      end

      subroutine par_stop (string)
c  **************************************************
c  Abnormal termination of the parallel session
c  **************************************************
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
c     For abnormal exits 
c     If open, close unit = 11
c     Go to the barrier that workers are sitting at
c     Then everyone will call par_end and stop
      logical is_open
      character*(*) string

      inquire(unit=11,opened=is_open)
      if (is_open) then
        call wlog(string)
        close(unit=11)
      else if (string .ne. ' ') then
	print *,string
	print *,'Abnormal termination on processor ',this_process
      endif
      call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierrorflag)

      stop ' '
      end

      subroutine par_end
c  **************************************************
c  Terminate the parallel session
c  **************************************************
      call MPI_FINALIZE(ierrorflag)
      return
      end

      subroutine par_barrier
c  **************************************************
c  Calls mpi_barrier
c  **************************************************
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BARRIER(MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_int(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for integer arrays
c  **************************************************
      integer count,dest,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_INTEGER,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_cmplx(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for complex arrays
c  **************************************************
      integer count,dest,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_dc(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for double_complex arrays
c  **************************************************
      integer count,dest,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_DOUBLE_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_recv_int(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for integer arrays
c  **************************************************
      integer count,source,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_INTEGER,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_cmplx(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for complex arrays
c  **************************************************
      integer count,source,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_dc(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for double complex arrays
c  **************************************************
      integer count,source,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_DOUBLE_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_bcast_int(buf,count,source)
c  **************************************************
c  Call mpi_bcast for integer arrays
c  **************************************************
      integer count,source
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_INTEGER,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_cmplx(buf,count,source)
c  **************************************************
c  Call mpi_bcast for complex arrays
c  **************************************************
      integer count,source
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_dc(buf,count,source)
c  **************************************************
c  Call mpi_bcast for double_complex arrays
c  **************************************************
      integer count,source
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_DOUBLE_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine MPE_DECOMP1D( n, num_procs, myid, s, e )
c  ******************************************************
c  A routine for producing a decomposition of a 1-d 
c  array when given a number of processors.  It may 
c  be used in "direct" product decomposition.  The 
c  values returned assume a "global" domain in [1:n]
c  ******************************************************
c  MPE_Decomp1d - Compute a balanced decomposition of
c  a 1-D array
c  ******************************************************
c  Input Parameters:
c  n  - Length of the array
c  num_procs - Number of processors in decomposition
c  myid  - Rank of this processor in the decomposition 
c  (0 <= rank < size)
c  ******************************************************
c  Output Parameters:
c  s,e - Array my_particles are s:e, with the original 
c  array considered as 1:n.  
c  ******************************************************

      integer n, num_procs, myid, s, e
      integer nloc, deficit
 
      nloc  = n / num_procs
      s       = myid * nloc + 1
      deficit = mod(n,num_procs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
        nloc = nloc + 1
      endif
      e = s + nloc - 1
      if (e .gt. n .or. myid .eq. num_procs-1) e = n

      return
      end

      SUBROUTINE SECONDS( W)
c  ***************************************************
c  SECONDS returns the wall clock times for a process
c  in seconds.
c  ***************************************************

      real*8 W, MPI_Wtime

      W = MPI_Wtime()

      RETURN
      END
c///////////////////////////////////////////////////////////////////////
c Distribution:  FEFF_MATH 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: FEFF_MATH 1.0
c      FEFF_MATH 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,
     1                 kind, lind, bmat)
c     written by alexei ankudinov; march 2000
c     calculate bmat: the energy independent sum over polarization and
c     angular momenta indices
c     bmat = \sum_{p,p', all m_j} <LS|J><J|R|J1><J1|\alpha_p exp(i kz)|I>
c                    ptz(p,p') 
c            <I|\alpha_p'^* exp(-i kz) J2><J2|R'|J'><J'|L'S'>
c     where R is rotation from spin vector to x-ray k-vector
c     and R' is rotation back
c     see Eq.10 and 11 in Ankudinov,Rehr, Phys.Rev.B (accepted),
c     Theory of solid state contribution to the x-ray elastic scattering
c     aditional rotation matrices are needed when x-ray k-vector
c     is not along the spin-axis (see rotations in rdinp)

c     more precisely it is
c     bmat(l1 l1' j l ml ms; l2 l2' j' l' ml' ms') =
c        (-)**(j-j'+l2'+1) i**(l'-l) \sum_{p,p',mi,m1,mj,m2,mj'}
c        <LS|J>   r^j_{m1,mj}(angks)   3j( j l1 i ; -m1 p mi)
c        (-p)**(l1+l1'+1) ptz(p,p') (-p')**(l2+l2'+1) 
c        3j( j' l2 i ; -m2  p' mi)   r^j'_{m2,mj'}(angks)   <J'|L'S'>
c     where l1 l1' are set by the multipole moment(E1-l1=1,l1'=0;
c     E2-l1=2,l1'=1; M1-l1=1,l1'=1; etc.;
c     j and l define quantum number kappa and for each multipole moment
c     Only few final kappa's are allowed and  it is convinient
c     to denote (l1 l1' j l) by one index 'k'
c     thus  k=1-8 to include both E1 and E2 transitions;
c     ml and ms are projections of orbital and spin moments.

c     bmat  is used to calculate absorption fine structure (chi) via
c       chi = \sum_{k ms,k' ms'}  rkk(k,ms)  rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms')  G_(l' ml' ms';l ml ms)
c     where sum over spins can be moved from first sum to second for
c     spin independent systems. The above expression is suitable for FMS
c     and for MS expansion on can use Eq.15 in RA paper to obtain
c     expression for the termination   matrix
c     T_{lam1 ms,lamN ms'} = \sum_{k k'} rkk(k,ms) rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms') gam(l,lam1,rho1,ms)
c        gamtl(l',lamN,rhoN,ms')
c     Notice that for spin-dependent systems the scattering F matrices
c     in RA paper also should have additional spin indices. In genfmt
c     we currently neglect spin-flip processes which simplifies
c     calculations with MS expansion. (T and F are diagonal in ms,ms')
       
c     This subroutine is written for general spin-dependent asymmetric
c     system and arbitrary polarization tenzor. The symmetry of the 
c     system and polarization tenzor can be used
c     to speed up FMS or MS calculations in appropriate subroutines.
c     (see comments in subroutines mpprmp, fmstot)

c     input:
c       kinit - kappa for initial orbital
c       ipol - polarization type measurement
c       ptz  - polarization tensor (needed only for ipol=1 case)
c       le2  - sets which multipole moments to include (see mkptz)
c       ltrace- .true. for xsect.f, where need to perform trace over ml
c       angks - angle between k-vector and spin-vector 

c     output
c       lind  - orb.mom.(kappa)  needed in fmstot only (for indexing)
c       bmat  - energy independent matrix to calculate absorption 
c       in many cases bmat is diagonal due to the choice of xyz frame,
c       but for general case full 16*(2*lx+1)*16*(2*lx+1) matrix is kept

      implicit double precision (a-h,o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      complex*16 coni
      parameter (coni = (0,1))

c     need only parameter lx to set max orb momentum
      complex*16 ptz, bmat, pmat, tmat
      dimension ptz(-1:1,-1:1),  bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
c       to include all possible dipole and quadrupole transitions 
c       final kp, and kpp have 8 possibilities
      logical ltrace

c     local staff
      dimension  t3j( 8, 0:1, -lx:lx+1), x3j(8, -1:1, -lx:lx+1)
c     qmat = <J2|R'|J'><J'|L'S'> - diagonal in kappa index
      dimension qmat( -lx:lx+1, -lx:lx, 0:1, 8)
c     pmat = <J1|\alpha_j exp(i kz)|I> ptz <I|\alpha_k^* exp(-i kz)|J2>
      dimension pmat( -lx:lx+1, 8, -lx:lx+1, 8)
c     tmat = pmat*qmat ; bmat = qmat^T*tmat
      dimension tmat( -lx:lx+1, 8, -lx:lx, 0:1, 8)
c     total and orbital momenta for 8 possible final kappa
      dimension jind(8), lind(8), kind(8)

      external cwig3j

      do 10 i6 = 1, 8
      do 10 i5 = 0 ,1
      do 10 i4 = -lx,lx
      do 10 i3 = 1, 8
      do 10 i2 = 0 ,1
      do 10 i1 = -lx,lx
         bmat( i1, i2, i3, i4, i5, i6) = 0
  10  continue

c     3 dipole transitions
      do 20 k=-1,1
         kap=kinit+k
         if (k.eq.0) kap=-kap
         jkap = abs(kap)
         lkap = kap
         if (kap.le.0) lkap = abs(kap) -1
c        check that orbital momentum does not exceed max allowed
         if (lkap .gt. lx) then
c          set final j and l to unphysical values
           jkap = 0
           lkap = -1 
           kap = 0
         endif
         jind(k+2) = jkap
         lind(k+2) = lkap
         kind(k+2) = kap
  20  continue

c     include 5 quadrupole or 3 mag.dipole  transitions
      do 120 k=-2,2
         jkap = abs(kinit) + k
         if (jkap.le.0) jkap = 0
         kap= jkap
         if (kinit.lt.0 .and. abs(k).ne.1) kap=-jkap
         if (kinit.gt.0 .and. abs(k).eq.1) kap=-jkap
         lkap = kap
         if(kap.le.0) lkap = - kap - 1
         if (lkap.gt.lx .or. le2.eq.0
     1                  .or. (le2.eq.1 .and. abs(k).eq.2)) then
c           set unphysical jkap and lkap to make shorter calculations
            jkap = 0
            lkap = -1
            kap = 0
         endif
         jind(k+6) = jkap
         lind(k+6) = lkap
         kind(k+6) = kap
 120  continue

      if (ipol.eq.0) then
c       polarization average case; bmat is diagonal and simple
        do 100 k = 1, 8
        do 100 ms = 0 ,1
        do 100 ml = -lind(k), lind(k)
c         i2 = (2*l1+1) , where l1 is defined by multipole moment
          i2 = 3
          if (le2.eq.2 .and. k.gt.3) i2 = 5
          bmat(ml,ms,k, ml,ms,k) = 0.5d0 / (2*lind(k)+1.d0) / i2
          if (k.le.3) bmat(ml,ms,k, ml,ms,k) = - bmat(ml,ms,k, ml,ms,k)
 100    continue
      else
c       more complicated bmat for linear(ipol=1) and circular(ipol=2)
c       polarizations
c       Put 3j factors in x3j and t3j. t3j are multiplied by
c       sqrt(2*j'+1) for  further convinience.
        do 30  mp=-lx,lx+1
        do 30  ms=0,1
        do 30  k1=1,8
  30    t3j(k1,ms,mp) = 0.0d0
        do 40  mp=-lx,lx+1
        do 40  ms=-1,1
        do 40  k1=1,8
  40      x3j(k1,ms,mp) = 0.0d0

        do 70  k1 = 1,8
        do 70  mp = -jind(k1)+1,jind(k1)
          do 50 ms=0,1
            j1 = 2 * lind(k1)
            j2 = 1
            j3 = 2 * jind(k1) - 1
            m1 = 2*(mp-ms)
            m2 = 2*ms - 1
            t3j(k1,ms,mp)=sqrt(j3+1.0d0) * cwig3j(j1,j2,j3,m1,m2,2)
            if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0) 
     1          t3j(k1,ms,mp) = - t3j(k1,ms,mp)
c           t3j(m0,i)    are Clebsch-Gordon coefficients
  50      continue
          do 60 i=-1,1
            j1 = 2 * jind(k1) - 1
            j2 = 2
            if (k1.gt.3 .and. le2.eq.2) j2 = 4
            j3 = 2 * abs(kinit) - 1
            m1 = -2*mp + 1
            m2 = 2*i
            x3j(k1,i,mp)= cwig3j(j1,j2,j3,m1,m2,2)
  60      continue
  70    continue

c       calculate qmat
        do 220 i=1,8
        do 220 ms=0,1
        do 220 ml= -lind(i), lind(i)
        do 220 mj= -jind(i)+1, jind(i)
          mp = ml+ms
          jj = 2*jind(i) - 1
          mmj = 2*mj - 1
          mmp = 2*mp - 1
          value = rotwig(angks, jj, mmj, mmp, 2)
          qmat(mj,ml,ms,i) = value * t3j(i,ms,mp)
 220    continue

c       calculate pmat
        do 240 i2 = 1,8
        do 240 m2 = -jind(i2)+1, jind(i2)
        do 240 i1 = 1,8
        do 240 m1 = -jind(i1)+1, jind(i1)
          pmat(m1,i1,m2,i2) = 0
          if (abs(m2-m1).le.2) then
            do 230 j=-1,1
            do 230 i=-1,1
c             check that initial moment is the same
              if (m1-i.eq.m2-j) then
                is = 1
c               (-p) factors for M1 transitions
                if (le2.eq.1 .and. i.gt.0 .and. i1.gt.3) is = -is
                if (le2.eq.1 .and. j.gt.0 .and. i2.gt.3) is = -is
                pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) +
     1          is * x3j(i1,i,m1) * ptz(i,j) * x3j(i2,j,m2)
              endif
 230        continue
c           multiply by (-)^(j-j'+l2'+1) i**(l'-l) factor
c           additional (-) is from Eq.10 (-2*ck)
            is = 1
            if (mod(jind(i1)-jind(i2), 2) .ne.0) is = -is
            if (i2.le.3) is = -is
            pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) * is
     1           * coni**(lind(i2)-lind(i1))
          endif
 240    continue

c       calculate tmat = pmat*qmat
        do 270 i1=1,8
        do 270 ms=0,1
        do 270 ml=-lind(i1), lind(i1)
        do 270 i2=1,8
        do 270 mj=-jind(i2)+1, jind(i2)
          tmat(mj,i2, ml,ms,i1) = 0
          do 260 mp = -jind(i1)+1, jind(i1)
            tmat(mj,i2, ml,ms,i1) = tmat(mj,i2, ml,ms,i1)+
     1           pmat(mj,i2,mp,i1) * qmat(mp,ml,ms,i1)
 260      continue
 270    continue
         
c       calculate bmat = qmat^T * tmat
        do 300 i1=1,8
        do 300 ms1=0,1
        do 300 ml1=-lind(i1), lind(i1)
        do 300 i2=1,8
        do 300 ms2=0,1
        do 300 ml2=-lind(i2), lind(i2)
          bmat(ml2,ms2,i2, ml1,ms1,i1) = 0
          do 280 mj=-jind(i2)+1, jind(i2)
            bmat(ml2,ms2,i2, ml1,ms1,i1) = bmat(ml2,ms2,i2, ml1,ms1,i1)+
     1      qmat(mj,ml2,ms2,i2) * tmat(mj,i2,ml1,ms1,i1) 
 280      continue
 300    continue
c       end of ipol=1,2 cases
      endif 

      if (ltrace) then
c       need to trace bmat over ml for xsect.f
        do 390 i1 = 1, 8
        do 390 ms1 = 0,1
        do 390 i2 = 1, 8
        do 390 ms2 = 0,1
          if (lind(i1).ne.lind(i2) .or. ms1.ne.ms2) then
               bmat(0,ms2,i2, 0,ms1,i1) = 0
          else
             do 360 ml = 1, lind(i1)
               bmat(0,ms1,i2, 0,ms1,i1) =  bmat(0,ms1,i2, 0,ms1,i1) +
     1         bmat(-ml,ms1,i2, -ml,ms1,i1) + bmat(ml,ms1,i2, ml,ms1,i1)
 360         continue
          endif
 390    continue
      endif

      if (ispin .eq. 0) then
c       G(Ls,L's') is spin diagonal; trace over spin
        do 480 i1 = 1, 8
        do 480 i2 = 1, 8
        do 480 ml1 = -lind(i1), lind(i1)
        do 480 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,0,i2, ml1,0,i1) +
     1                                  bmat(ml2,1,i2, ml1,1,i1)
 480    continue
      elseif (ispin.eq.2 .or. (ispin.eq.1 .and. nspx.eq.1)) then
c       move spin up part into the position of spin-down
        do 490 i1 = 1, 8
        do 490 i2 = 1, 8
        do 490 ml1 = -lind(i1), lind(i1)
        do 490 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,1,i2, ml1,1,i1)
 490    continue

      endif

      return
      end
      subroutine besjn (x, jl, nl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 30 (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this array nl = abramowitz yl.
c       jl and nl must be dimensioned 
c            complex*16 jl(ltot+2), nl(ltot+2), with ltot defined in 
c            dim.h.
c
c     notes:  jl and nl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 cjl(ltot+2), sjl(ltot+2), cnl(ltot+2), snl(ltot+2)

      complex*16 xjl,xnl,asx,acx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)

      if (dble(x) .le. 0)  stop 'Re(x) is .le. zero in besjn'

      lmaxp1 = ltot+2

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 il = 1,lmaxp1
            l = il-1
            ifl = 0
            call bjnser (x,l,xjl,xnl,ifl)
            jl(il) = xjl
            nl(il) = xnl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmaxp1-1,xjl,xnl,1)
         jl(lmaxp1) = xjl

         call bjnser (x,lmaxp1-2,xjl,xnl,1)
         jl(lmaxp1-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(lp1) = tlxp3 * jl(lp1+1) / x  -  jl(lp1+2)
   60    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 80 ie = 1,11
            snl(ie) = cjl(ie)
            cnl(ie) = -sjl(ie)
   80    continue
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
            snl(lp1) = tlxp1*xi*snl(lp1-1)-snl(lp1-2)
            cnl(lp1) = tlxp1*xi*cnl(lp1-1)-cnl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         do 110 lp1 = 1,lmaxp1
            jl(lp1) = asx*sjl(lp1)+acx*cjl(lp1)
            nl(lp1) = asx*snl(lp1)+acx*cnl(lp1)
  110    continue
      endif

      return
      end
      subroutine besjh (x, lbmax, jl, hl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and hl
c               for l = 0 to lbmax (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       lbmax
c       jl = jl bessel function (abramowitz conventions)
c       hl = hl^+ bessel function (messiah conventions) for Im x >=0
c       hl = hl^- bessel function (messiah conventions) for Im x < 0
c       jl and hl must be dimensioned 
c            complex*16 jl(0:lbmax), hl(0:lbmax), 
c
c     notes:  jl and hl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c     rewritten for jl and hl by a.l. ankudinov feb 2000
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(0:lbmax), nl(ltot+2)
      complex*16 hl(0:lbmax)
      complex*16 cjl(ltot+2), sjl(ltot+2)

      complex*16 xjl,xnl,asx,acx, epx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)
      complex*16 coni
      parameter (coni=(0,1))

      if (dble(x) .lt. 0)  stop 'Re(x) is .lt. zero in besjh'

      lmax = min(lbmax, ltot+1)
      lmaxp1 = lmax + 1

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 ll = 0,lmax
            ifl = 0
            call bjnser (x,ll,xjl,xnl,ifl)
            jl(ll) = xjl
            hl(ll) = -xnl + coni*xjl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmax,xjl,xnl,1)
         jl(lmax) = xjl

         call bjnser (x,lmax-1,xjl,xnl,1)
         jl(lmax-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(l) = tlxp3 * jl(l+1) / x  -  jl(l+2)
   60    continue

         do 65 il = 1, lmaxp1
            l = il - 1
            hl(l) = -nl(il) + coni*jl(l)
   65    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         if (dimag(x).ge. 0.d0) then
           epx = exp(coni*x)
         else 
           epx = exp(-coni*x)
         endif
         do 110 ll = 0,lmax
            lp1 = ll + 1
            jl(ll) = asx*sjl(lp1)+acx*cjl(lp1)
            if (dimag(x).ge. 0.d0) then
              hl(ll) = (sjl(lp1)+coni*cjl(lp1)) * epx
            else
              hl(ll) = (sjl(lp1)-coni*cjl(lp1)) * epx
            endif
  110    continue
      endif

      return
      end
      subroutine bjnser (x, l, jl, nl, ifl)

c-----------------------------------------------------------------------
c
c     subroutine: bjnser (x,l,jl,nl,ifl)
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c
c     arguments:
c       x = argument of jl and nl
c       l = l value calculated (no offset)
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c       ifl = 0 return both jl and nl
c             1 return jl only
c             2 return nl only
c
c     notes:  jl and nl are calculated by a series
c             expansion according to 10.1.2 and 10.1.3
c             in abramowitz and stegun (ninth printing),
c             page 437
c
c             error msgs written with PRINT statements.
c
c     first coded by r. c. albers on 26 jan 83
c
c     version 2
c
c     last modified: 27 jan 83 by r. c. albers
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      character*512 slog

      parameter (niter = 160, tol = 1.e-15)

      if (l .lt. 0) then
         call wlog(' l .lt. 0 in bjnser')
         stop 'bjnser 1'
      endif
      if (dble(x).lt. 0.) then
         write(slog,30) x
         call wlog(slog)
   30    format (' x = ', 1p, 2e14.6, ' is .le. 0 in bjnser')
         stop 'bjnser 2'
      endif

      lp1 = l+1
      u = x**2 / 2

c     make djl = 1 * 3 * 5 * ... * (2*l+1),
c          dnl = 1 * 3 * 5 * ... * (2*l-1)
      djl = 1
      fac = -1
      do 50 il = 1, lp1
         fac = fac + 2
         djl = fac * djl
   50 continue
      dnl = djl / (2*l+1)


      if (ifl .eq. 2)   goto 90
c     make jl
c     pj is term in { } in 10.1.2, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pj = 1
      nf = 1
      nfac = 2*l + 3
      den = nfac
      sgn = -1
      ux = u
      do 60 il = 1, niter
         del = sgn*ux / den
         pj = pj + del
         trel = abs (del / pj)
         if (trel .le. tol)  goto 80
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
   60 continue
      stop  'jl does not converge in bjnser'
   80 jl = pj * (x**l) / djl

   90 if (ifl.eq.1) return
c     make nl
c     pn is term in { } in 10.1.3, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pn = 1
      nf = 1
      nfac = 1 - 2*l
      den = nfac
      sgn = -1
      ux = u
      do 100  il = 1, niter
         del = sgn * ux / den
         pn = pn + del
         trel = abs (del / pn)
         if (trel .le. tol) goto 120
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
  100 continue
      stop  'nl does not converge in bjnser'
  120 nl = -pn * dnl / (x**lp1)

      return
      end
      subroutine conv(omega,xsec,ne1,vicorr)
c     multiply xsec by theta(omega-efermi) and
c     convolute xsec(omega) with  xloss/((omega-omega0)**2+xloss**2)/pi
c     the result is xsec0(omega0)

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      dimension  omega(nex)
      complex*16 xsec(nex), xsec0(nex), xsecdx

      complex*16 conv1
      external conv1

      do 100 ie = 1,ne1
         xsec0(ie) = 0.0d0
         omega0 = omega(ie)
c        Add one more point to correct for the finite grid
c        at large energies. Use linear interpolation.
         dx = max( omega(ne1) - omega(ne1-1), 50*vicorr)
         xlast = omega(ne1)+dx
         dx = dx / ( omega(ne1) - omega(ne1-1))
         xsecdx = xsec(ne1)+ (xsec(ne1)-xsec(ne1-1)) * dx

c        first interval
         do 50  i = 1, ne1-1
            xsec0(ie) = xsec0(ie) +
     1      conv1(omega(i),omega(i+1),xsec(i),xsec(i+1),omega0,vicorr)
  50     continue
c        last interval
         xsec0(ie) = xsec0(ie) +
     1   conv1(omega(ne1),xlast,xsec(ne1),xsecdx,omega0,vicorr)
         xsec0(ie) = xsec0(ie) /real(pi)
  100 continue
      do 200 ie = 1, ne1
  200 xsec(ie) = xsec0(ie)

      return
      end

      complex*16 function conv1(x1,x2,y1,y2,x0,xloss)
c     convolution of function 1/(omega-omega0-i*xloss)/pi
c     makes linear interpolation for function between x1,x2 and
c     takes advantage that the integral can be taken analytically.
      implicit double precision (a-h, o-z)
      complex*16  y1, y2, t, coni,dum, a, b
      parameter (coni = (0.0,1.0))

      d = (x2-x1) / 2.0
      a = dble(y2-y1) / 2.0
      b = dble(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = dimag (dum)

      d = (x2-x1) / 2.0
      a = dimag(y2-y1) / 2.0
      b = dimag(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = conv1 + coni* dimag( dum)

      return
      end
      subroutine cpl0 (x, pl0, lmaxp1)
      implicit double precision (a-h, o-z)

c-----------------------------------------------------------------------
c
c     cpl0:  Calculate associated legendre polynomials p_l0(x)
c            by recursion.
c            Adapted from aslgndr.
c
c     first written: (25 june 86) by j. j. rehr
c
c     version 1 (25 june 86) (aslgndr)
c     version 2 (March, 1992) siz
c
c-----------------------------------------------------------------------

      dimension pl0 (lmaxp1)

      lmax = lmaxp1-1

c     calculate legendre polynomials p_l0(x) up to l=lmax
      pl0(1) = 1
      pl0(2) = x
      do 10  il = 2, lmax
         l = il-1
         pl0(il+1) = ( (2*l+1)*x*pl0(il) - l*pl0(l) ) / il
   10 continue

      return
      end
      subroutine csomm (dr,dp,dq,dpas,da,m,np)
c Modified to use complex p and q.  SIZ 4/91
c integration by the method of simpson of (dp+dq)*dr**m from 
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),dq(*),da,dc
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      da=da+dc
      dc=dq(i)*dl
      da=da+dc
   70 continue
      da=dpas*da/3
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine csomm2 (dr,dp,dpas,da,rnrm,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),da,dc

      d1=dble(da)+1
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)/db
      dd=(dr(1))*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)
c     wigner 3j coefficient for integers  (ient=1)
c                         or semiintegers (ient=2)
c     other arguments should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
      character*512 slog
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      m3=-m1-m2
      if (ini) 1,21,1
c        initialisation of the log's of the factorials
 1    ini=0
      al(1)=0.0d 00
      do 11 i=1,idim
         b=i
 11      al(i+1)=al(i)+ log(b)
 21   cwig3j=0.0d 00
      if (((ient-1)*(ient-2)).ne.0) go to 101
      ii=ient+ient
c        test triangular inequalities, parity and maximum values of m
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99
      m(1)=j1+j2-j3
      m(2)=j2+j3-j1
      m(3)=j3+j1-j2
      m(4)=j1+m1
      m(5)=j1-m1
      m(6)=j2+m2
      m(7)=j2-m2
      m(8)=j3+m3
      m(9)=j3-m3
      m(10)=j1+j2+j3+ient
      m(11)=j2-j3-m1
      m(12)=j1-j3+m2
      do 41 i=1,12
         if (i.gt.10) go to 31
         if (m(i).lt.0) go to 99
 31      if (mod(m(i),ient).ne.0) go to 101
         m(i)=m(i)/ient
         if (m(i).gt.idim) go to 101
 41   continue

c        calculate 3j coefficient
      max0= max(m(11),m(12),0)+1
      min0= min(m(1),m(5),m(6))+1
      isig=1
      if (mod(max0-1,2).ne.0) isig=-isig
      c=-al(m(10)+1)
      do 61 i=1,9
 61   c=c+al(m(i)+1)
      c=c/2.0d 00
      do 71 i=max0,min0
      j=2-i
      b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))
      cwig3j=cwig3j+isig* exp(c-b)
 71   isig=-isig
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j
 99   return
 101     write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient
         call wlog(slog)
      stop
      end
      double precision function determ(array,nord,nrows)
c
c     calculate determinate of a square matrix
c        (from bevington "data reduction and error analysis
c         for the physical sciences" pg 294)
c     array: matrix to be analyzed
c     nord: order of matrix
c     nrows:  first dimension of matrix in calling routine
c
      double precision array(nrows,nrows)
      determ = 1.
      do 150 k=1,nord
c
c
        if (array(k,k).ne.0) go to 130
        do 100 j=k,nord
          if (array(k,j).ne.0) go to 110
  100   continue
        determ = 0.
        go to 160
c
  110   do 120 i=k,nord
          saved = array(i,j)
          array(i,j) = array(i,k)
  120   array(i,k) = saved
        determ = -determ
c
  130   determ = determ*array(k,k)
        if (k.ge.nord) go to 150
        k1 = k+1
        do 140 i=k1,nord
          do 140 j=k1,nord
  140   array(i,j) = array(i,j)-array(i,k)*array(k,j)/array(k,k)
  150 continue
  160 return
c end double precision function determ
      end
      double precision function dist (r0, r1)
c     find distance between cartesian points r0 and r1
      implicit double precision (a-h, o-z)
      dimension r0(3), r1(3)
      dist = 0
      do 10  i = 1, 3
         dist = dist + (r0(i) - r1(i))**2
   10 continue
      dist = sqrt (dist)
      return
      end
      double precision function rotwig (beta, jj, m1, m2, ient)
c     uses Wigner formula (Messiah eq.C.72) to calculate rotation matrix
c     for integers  (ient=1)  or semiintegers (ient=2)
c     other arguments (except beta) should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      if (((ient-1)*(ient-2)).ne.0) stop ' Illegal ient in rotwig.'

      if (ini.eq.1) then
c       initialisation of the log's of the factorials
        ini=0
        al(1)=0.0d 00
        do 11 i=1,idim
           b=i
 11        al(i+1)=al(i)+ log(b)
      endif
      rotwig = 0.d0

      if ( m1.ge.0 .and. abs(m1).ge.abs(m2)) then
         m1p = m1 
         m2p = m2
         betap = beta
         isign = 1
      elseif (m2.ge.0 .and. abs(m2).ge.abs(m1)) then
         m1p = m2
         m2p = m1
         betap = - beta
         isign = 1
      elseif (m1.le.0 .and. abs(m1).ge.abs(m2)) then
         m1p = - m1
         m2p = - m2
         betap = beta
         isign = (-1)**( (m1-m2)/ient ) 
      else
         m1p = - m2
         m2p = - m1
         betap = - beta
         isign = (-1)**( (m2-m1)/ient ) 
      endif

      temp = 0.d0
      zeta = cos ( betap / 2.d0 )
      eta  = sin ( betap / 2.d0 )
      do 100 it = m1p - m2p, jj - m2p, ient
        m(1) = 1 + (jj+m1p) / ient
        m(2) = 1 + (jj-m1p) / ient
        m(3) = 1 + (jj+m2p) / ient
        m(4) = 1 + (jj-m2p) / ient
        m(5) = 1 + (jj+m1p-it) / ient
        m(6) = 1 + (jj-m2p-it) / ient
        m(7) = 1 + it / ient
        m(8) = 1 + (m2p-m1p+it) / ient
        m(9)  = (2*jj+m1p-m2p-2*it) / ient 
        m(10) = (2*it-m1p+m2p) / ient 
        factor = 0.d0
        do 110 i = 1,4
  110     factor = factor + al(m(i))/2.d0 - al(m(i+4))
c       special cases to resolve 0.d0**0 problem (a.ankudinov, may 2001)
        if (m(10).eq.0 .and. m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*exp(factor)
        elseif (m(10).eq.0) then
          temp = temp + (-1)**(it/ient)*zeta**m(9)*exp(factor)
        elseif (m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*eta**m(10)*exp(factor)
        else
c         general expression
          temp = temp+ (-1)**(it/ient)*zeta**m(9)*eta**m(10)*exp(factor)
        endif
  100 continue

      rotwig = isign * temp
     
      return
      end
      subroutine phamp (rmt, pu, qu, ck, jl, nl, jlp, nlp, ikap,
     1                  ph, amp)
c     calculate phase shift at mt radius
c     needs to calculate atan of complex variable (coded below)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      external besjn, atan2c

      complex*16 pu, qu, ck,  jl, nl, jlp, nlp, ph, amp
      complex*16 xkr, a, b, factor

c     initialize staff
      xkr = ck*rmt
      isign=1
      if (ikap.lt.0) isign = -1
      a = ck*alphfs
      factor = isign*a/(1+sqrt(1+a**2))

c     find a and b that pu = rmt*(a*jl+b*nl), qu=factor*rmt*(a*jlp+b*nlp)
      a = isign*ck*xkr* (pu*nlp - qu*nl/factor)
      b = isign*ck*xkr* (qu*jl/factor - pu*jlp)

c     pu =  amp * rmt * (jl*cos(ph) - nl*sin(ph))
c     qu =  amp * rmt * (jlp*cos(ph) - nlp*sin(ph)) * factor
c     tan(ph) = - b/a
      b = -b
      call atan2c ( a, b, amp, ph)

      return
      end
      subroutine atancc(temp, phx)
c     phx=atan(temp), for complex numbers
      implicit double precision (a-h, o-z)
      complex*16 temp, phx

      xx = dble (temp)
      yy = dimag(temp)
      if (xx .ne. 0)  then
         alph = (1 - xx**2 - yy**2)
         alph = sqrt(alph**2 + 4*xx**2) - alph
         alph = alph / (2 * xx)
         alph = atan (alph)
      else
         alph = 0
      endif
      beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
      beta = log(beta) / 4
      phx = dcmplx (alph, beta)

      return
      end

      subroutine atan2c(a, b, ampl, phx)
c     for complex a, b find complex ampl, phx such that:
c     a= ampl*cos(phx)  and  b= ampl*sin(phx)
c     phx=atan(b/a)
      implicit double precision (a-h, o-z)
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      complex*16 a, b, ampl, phx, temp

      aa = abs(a)
      bb = abs(b)
      if (aa+bb.eq. 0) then
         ampl=0.d0
         phx =0.d0
      elseif ( aa.gt.bb) then
         temp = b/a
         call atancc ( temp, phx)
         ampl = a / cos(phx)
      else
         temp = a/b
         call atancc ( temp, phx)
         phx = pi / 2 - phx
         ampl = b/sin(phx)
      endif

      if (dble(ampl).lt. 0.d0) then
         ampl = -ampl
         phx = phx + pi
      endif

      return
      end
      subroutine exjlnl (z, l, jl, nl)

c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 6  using exact analytic expression
c
c     arguments:
c       z = argument of jl and nl
c       l = integer order of spherical bessel function
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this nl = abramowitz yl.
c
c       analytic expressions from abramowitz 10.1.11 and 10.1.12
c       recurrence relation to get analytic j4,n4  eqns 10.1.19-22 ala

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

c     Exact formulae unstable for very small z, so use series
c     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3)  then
         call bjnser (z, l, jl, nl, 0)
      else
c        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3/z**3 - 1/z)*sinz - 3*cosz/z**2
            nl = (-3/z**3 + 1/z)*cosz - 3*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15/z**4 - 6/z**2)*sinz + (-15/z**3 + 1/z)*cosz
            nl = (-15/z**4 + 6/z**2)*cosz + (-15/z**3 + 1/z)*sinz

         elseif (l .eq. 4)  then
            jl = ( 105/z**5 - 45/z**3 + 1/z )*sinz + 
     1                ( -105/z**4 + 10/z**2 )*cosz
            nl = (-105/z**5 + 45/z**3 - 1/z )*cosz + 
     1                ( -105/z**4 + 10/z**2 )*sinz

         elseif (l .eq. 5)  then
            jl = ( 945/z**6 - 420/z**4 + 15/z**2 )*sinz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*cosz
            nl = (-945/z**6 + 420/z**4 - 15/z**2 )*cosz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*sinz

         elseif (l .eq. 6)  then
            jl = ( 10395/z**7 - 4725/z**5 + 210/z**3 - 1/z )*sinz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*cosz
            nl = (-10395/z**7 + 4725/z**5 - 210/z**3 + 1/z )*cosz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*sinz

         else
            stop 'exjlnl, l out of range'
         endif
      endif

      return
      end
      subroutine polint( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) pause 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      function sdist (r0, r1)
c     find distance squared between cartesian points r0 and r1
c     single precision
      dimension r0(3), r1(3)
      sdist = 0
      do 10  i = 1, 3
         sdist = sdist + (r0(i) - r1(i))**2
   10 continue
      sdist = sqrt(sdist)
      return
      end
      subroutine somm (dr,dp,dq,dpas,da,m,np)
c
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(np), dp(np), dq(np)
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      if (dc) 20,40,30
   20 db=db+dc
      go to 40
   30 da=da+dc
   40 dc=dq(i)*dl
      if (dc) 50,70,60
   50 db=db+dc
      go to 70
   60 da=da+dc
   70 continue
      da = dpas * (da + db) / 3.0
      dc=exp(dpas)-1.0
      db=d1*(d1+1.0)*dc*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.0+1.0/(dc*(d1+1.0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine somm2 (dr,dp,dpas,da,rnrm,m,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      dimension  dp(*)

      mm = m + 1
      d1=dble(da)+mm
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)**mm*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)**mm*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)**mm*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)**mm*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)**mm*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)**mm*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)**mm*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)**mm
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      subroutine strap (x, y, n, sum)

c     Trapeziodal integration of y(x), result in sum
c     SINGLE PRECISION
c     modified by ala to handle cases for E<Efermi
c     sum only positive numbers

      dimension x(n), y(n)

      sum = y(1) * abs(x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * abs(x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * abs(x(n) - x(n-1))
      sum = sum/2

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terp (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polint( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end


c     These routines, terp1 and locat1, are special versions to
c     be used with ff2chi, which uses some single and some double
c     precision.  They are the same as the routines in terp.f.

      subroutine terp1 (x, y, n, x0, y0)
      implicit double precision (a-h, o-z)

      real x(n), y(n)

c     Find out between which x points x0 lies
      i = locat1 (x0, n, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      i = max (i, 1)
      i = min (i, n-1)

      if (x(i+1) - x(i) .eq. 0)  stop 'TERP-1'

      y0 = y(i) +  (x0 - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))

      return
      end

      function locat1 (x, n, xx)
      integer  u, m, n
      double precision x
      real xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat1 = 0
      u = n+1

   10 if (u-locat1 .gt. 1)  then
         m = (u + locat1) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat1 = m
         endif
         goto 10
      endif

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terpc (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      complex*16 y, y0, dy
      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polinc( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      subroutine polinc( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      complex*16 ya,y,dy,c,d,w,den
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) stop 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      subroutine trap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2

      return
      end
      SUBROUTINE CQdrtc(Coef,Sol,NSol)
c     Combutes the zeros of a quadratic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(3)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a and b are zero
      COMPLEX*16 Sol(2)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 q, Sqrt
      DOUBLE PRECISION Sgn

      IF(Coef(1).eq.0.d0) THEN
         IF(Coef(2).eq.0.d0) THEN
            NSol = -1
            RETURN
         ELSE
            NSol = 1
            Sol(1) = -Coef(3)/Coef(2)
         END IF
      ELSE
         NSol = 2
         Root = Sqrt(Coef(2)**2-4.d0*Coef(1)*Coef(3))
         Sgn  = SIGN(DBLE(CONJG(Coef(2))*Root),1.d0)
         q    = -0.5d0*(Coef(2) + Sgn*Root)
         
         Sol(1) = q/Coef(1)
         Sol(2) = Coef(3)/q
      END IF

      RETURN
      END


      SUBROUTINE CCubic(Coef,Sol,NSol)
c     Combutes the zeros of a cubic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(4)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a, b, and c are zero
      COMPLEX*16 Sol(4)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 P1, P2, Q, R, Coef2(3), a, b, c
      DOUBLE PRECISION Sgn, Theta
c     PARAMETERS
      COMPLEX*16 I
      PARAMETER(I = (0.d0, 1.d0))
      DOUBLE PRECISION Pi
      PARAMETER(Pi = 3.141592653589793238462643d0)

      IF(Coef(1).eq.0.d0) THEN
         Coef2(1) = Coef(2)
         Coef2(2) = Coef(3)
         Coef2(3) = Coef(4)         
         CALL CQdrtc(Coef2,Sol,NSol)
      ELSE
         a = Coef(2)/Coef(1)
         b = Coef(3)/Coef(1)
         c = Coef(4)/Coef(1)
         NSol = 3
         Q = (a**2 - 3.d0*b)/9.d0
         R = (2.d0*a**3 - 9.d0*a*b + 27.d0*c)/54.d0

         IF(((DIMAG(Q).eq.0.d0).and.(DIMAG(R).eq.0.d0)).and.
     &        (DIMAG(R**2).lt.DIMAG(Q**3))) THEN
            Theta = ACOS (DBLE(R/SQRT(Q**3)))
            Sol(1) = -2*SQRT(Q)*Cos(Theta/3.d0) - a/3.d0
            Sol(2) = -2*SQRT(Q)*Cos((Theta+2.d0*Pi)/3.d0) - a/3.d0
            Sol(3) = -2*SQRT(Q)*Cos((Theta-2.d0*Pi)/3.d0) - a/3.d0
         ELSE
            Sgn = SIGN(1.d0, DBLE(CONJG(R)*SQRT(R**2-Q**3)))
            P1 = -(R + Sgn*SQRT(R**2-Q**3))**(1.d0/3.d0)
            IF(P1.eq.0.d0) THEN
               P2 = 0.d0
            ELSE
               P2 = Q/P1
            END IF
            Sol(1) = (P1 + P2) - a/3.d0
            Sol(2) = -0.5d0*(P1 + P2) - a/3.d0 +
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
            Sol(3) = -0.5d0*(P1 + P2) - a/3.d0 -
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
         END IF
      END IF

      RETURN
      END
      
