      subroutine json_read_pot(mpot, nph, ntitle, ihole, ipr1, iafolp,
     1       ixc, ispec, nmix, nohole, jumprm, inters, nscmt, icoul,
     2       lfms1, iunf, gamach, rgrd, ca1, ecv, totvol, rfms1,
     3       title, iz, lmaxsc, xnatph, xion, folp, novr,
     4       iphovr, nnovr, rovr )


      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      character*7 vname
      type(json_file) :: json   !the JSON structure read from the file:
      double precision toss
      integer,dimension(:),allocatable :: intgs
      character*80,dimension(:),allocatable :: strings
      double precision,dimension(:),allocatable :: dbpcs


      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    mod1.inp
        character*80 title(nheadx)
c        character*80 head(nheadx)
c        integer lhead(nheadx)
        integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec,
     1     iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
        integer iz(0:nphx), lmaxsc(0:nphx)
        real rfms1
        double precision gamach, rgrd, ca1, ecv, totvol
        double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c       for OVERLAP option
        integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
        double precision  rovr(novrx,0:nphx)
      

      call json%load_file('pot.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read pot.json"
         stop
      else
         call json%get('mpot',   mpot, found)
                   if (.not. found) call bailout('mpot', 'pot.json')
         call json%get('nph',    nph, found)
                   if (.not. found) call bailout('nph', 'pot.json')
         call json%get('ntitle', ntitle, found)
                   if (.not. found) call bailout('ntitle', 'pot.json')
         call json%get('ihole',  ihole, found)
                   if (.not. found) call bailout('ihole', 'pot.json')
         call json%get('ipr1',   ipr1, found)
                   if (.not. found) call bailout('ipr1', 'pot.json')
         call json%get('iafolp', iafolp, found)
                   if (.not. found) call bailout('iafolp', 'pot.json')
         call json%get('ixc',    ixc, found)
                   if (.not. found) call bailout('ixc', 'pot.json')
         call json%get('ispec',  ispec, found)
                   if (.not. found) call bailout('ispec', 'pot.json')

         call json%get('nmix',   nmix, found)
                   if (.not. found) call bailout('nmix', 'pot.json')
         call json%get('nohole', nohole, found)
                   if (.not. found) call bailout('nohole', 'pot.json')
         call json%get('jumprm', jumprm, found)
                   if (.not. found) call bailout('jumprm', 'pot.json')
         call json%get('inters', inters, found)
                   if (.not. found) call bailout('inters', 'pot.json')
         call json%get('nscmt',  nscmt, found)
                   if (.not. found) call bailout('nscmt', 'pot.json')
         call json%get('icoul',  icoul, found)
                   if (.not. found) call bailout('icoul', 'pot.json')
         call json%get('lfms1',  lfms1, found)
                   if (.not. found) call bailout('lfms1', 'pot.json')
         call json%get('iunf',   iunf, found)
                   if (.not. found) call bailout('iunf', 'pot.json')

         call json%get('gamach', gamach, found)
                   if (.not. found) call bailout('gamach', 'pot.json')
         call json%get('rgrd',   rgrd, found)
                   if (.not. found) call bailout('rgrd', 'pot.json')
         call json%get('ca1',    ca1, found)
                   if (.not. found) call bailout('ca1', 'pot.json')
         call json%get('ecv',    ecv, found)
                   if (.not. found) call bailout('ecv', 'pot.json')
         call json%get('totvol', totvol, found)
                   if (.not. found) call bailout('totvol', 'pot.json')
         call json%get('rfms1',  toss, found)   
                   if (.not. found) call bailout('rfms1', 'pot.json')
         rfms1 = real(toss)

         call json%get('titles', strings, found)
                   if (.not. found) call bailout('titles', 'pot.json')
         do 1000 itit = 1, nheadx
            title(itit) = strings(itit)            
 1000    continue
         call json%get('iz', intgs, found)
                   if (.not. found) call bailout('iz', 'pot.json')
         do 1010 iph = 0, nphx
            iz(iph) = intgs(iph+1)            
 1010    continue
         call json%get('lmaxsc', intgs, found)
                   if (.not. found) call bailout('lmaxsc', 'pot.json')
         do 1020 iph = 0, nphx
            lmaxsc(iph) = intgs(iph+1)            
 1020    continue
         call json%get('xnatph', dbpcs, found)
                   if (.not. found) call bailout('xnatph', 'pot.json')
         do 1030 iph = 0, nphx
            xnatph(iph) = dbpcs(iph+1)            
 1030    continue
         call json%get('xion', dbpcs, found)
                   if (.not. found) call bailout('xion', 'pot.json')
         do 1040 iph = 0, nphx
            xion(iph) = dbpcs(iph+1)            
 1040    continue
         call json%get('folp', dbpcs, found)
                   if (.not. found) call bailout('folp', 'pot.json')
         do 1050 iph = 0, nphx
            folp(iph) = dbpcs(iph+1)            
 1050    continue
         call json%get('novr', intgs, found)
                   if (.not. found) call bailout('novr', 'pot.json')
         do 1060 iph = 0, nphx
            novr(iph) = intgs(iph+1)            
 1060    continue

c        the following reconstructed all the overlap 2D arrays, see
c        RDINP/wrtjsn.f line 131 and following
         do 1200 iph = 0, nph
            write (vname, "(A6,I1)") "iphovr", iph
            call json%get(vname, intgs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1220 iovr = 1, novr(iph)
               iphovr(iovr, iph) = intgs(iovr)
 1220       continue

            write (vname, "(A5,I1)") "nnovr", iph
            call json%get(vname, intgs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1230 iovr = 1, novr(iph)
               nnovr(iovr, iph) = intgs(iovr)
 1230       continue

            write (vname, "(A4,I1)") "rovr", iph
            call json%get(vname, dbpcs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1240 iovr = 1, novr(iph)
               rovr(iovr, iph) = dbpcs(iovr)
 1240       continue
 1200    continue

         call json%destroy()
      end if

      return
      end
