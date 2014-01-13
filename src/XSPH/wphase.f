      subroutine wphase (nph, em, eref, lmax, ne, ph, ntitle, title)

c     Writes phase data to file PHASExx.DAT for each shell

      implicit double precision (a-h, o-z)

      include '../HEADERS/dim.h'

      complex*16 eref(nex, nspx)
      complex*16 ph( nex, -ltot:ltot, nspx, 0:nphx)
      complex*16  em(nex)
      dimension lmax(0:nphx)
      character*30  fname
      character*80  title(ntitle)
      character*2 coment
      parameter (coment='# ')

c     Dump phase data, eref and complex phase for each shell
      do 200  iph = 0, nph
         linit = 0
         if (linit .ge. lmax(iph)-1) linit = lmax(iph)-2
         if (linit .lt. 0) linit = 0

c        prepare files for shell's phase data

         write(fname,20)  iph
  20     format('phase', i2.2, '.dat')
         open (unit=1, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'wphase')

         write(fname,30)  iph
  30     format('phmin', i2.2, '.dat')
         open (unit=2, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'wphase')

         do 50 i = 1, ntitle
            ll = istrln(title(i))
            write(1,40)  coment, title(i)(1:ll)
            write(2,40)  coment, title(i)(1:ll)
  40        format (a,a)
  50     continue
c        write out unique pot and lmax
         write(1,60)   coment, iph, lmax(iph), ne
         write(2,60)   coment, iph, lmax(iph), ne
  60     format (a, 1x, 3i4, '   unique pot,  lmax, ne')
         write(2,70) coment, linit,linit+1,linit+2
  70     format (a,'      energy      re(eref)     re(p)    phase( ',i2,
     1         ')  phase(',i2,') phase(',i2,')' ) 

c        for each energy
c        ie, em, eref, p=sqrt(2*(em-eref))
c        ph array from 0 to ltot, 5 values per line
         do 150  ie = 1, ne
           write(1,110) coment, ie, dble(em(ie)), eref(ie,1),
     1                  sqrt(2*(em(ie)-eref(ie,1)))
  110      format (a, '   ie        energy      re(eref)',
     1             '      im(eref)',
     2             '         re(p)         im(p)', /,
     3             1x, i4, 1p, 5e14.6)

           write(1,120)  (ph(ie,ll,1,iph), ll=0,lmax(iph))
  120      format (1x, 1p, 4e14.6)

           write(2,130) dble(em(ie)), dble(eref(ie,1)),
     1     dble(sqrt(2*(em(ie)-eref(ie,1)))),
     2     (dble(ph(ie,ll,1,iph)), ll=linit,linit+2)
  130       format (1p, 6e13.5)
  150    continue
         close(unit=1)
         close(unit=2)
  200 continue

      return
      end
