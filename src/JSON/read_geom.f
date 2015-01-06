      subroutine json_read_geom(nat, nph, iatph, rat, iphat, ibounc)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      integer,dimension(:),allocatable :: intgs, intgb
      double precision,dimension(:),allocatable :: dbpcs, dbpcy, dbpcz

      integer nat, nph, iatph(0:nphx), iphat(natx), ibounc(natx)
      double precision  rat(3,natx)

      nph = 0
      do 10 iph = 0, nphx
         iatph(iph) = 0
 10   continue

      call json_initialize()
      call json%load_file(filename='geom.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read geom.json"
         stop
      else
         call json%get('natt',   nat, found)
                   if (.not. found) call bailout('natt', 'geom.json')
c         call json%get('nph',   nph, found)
c                   if (.not. found) call bailout('nph', 'geom.json')
c         call json%get('iatph', intgs, found)
c                   if (.not. found) call bailout('iatph', 'geom.json')
c         do 2000 iph = 0, nphx
c            iatph(iph) = intgs(iph+1)            
c 2000    continue
         call json%get('x', dbpcs, found)
                   if (.not. found) call bailout('x', 'geom.json')
         call json%get('y', dbpcy, found)
                   if (.not. found) call bailout('y', 'geom.json')
         call json%get('z', dbpcz, found)
                   if (.not. found) call bailout('z', 'geom.json')
         call json%get('iph', intgs, found)
                   if (.not. found) call bailout('iph', 'geom.json')
         call json%get('ibo', intgb, found)
                   if (.not. found) call bailout('ibo', 'geom.json')
         do 20 i=1,nat
            iphat(i)  = intgs(i)
            rat(1,i)  = dbpcs(i)
            rat(2,i)  = dbpcy(i)
            rat(3,i)  = dbpcz(i)
            ibounc(i) = intgb(i)
c     iatph is NOT set from the similar line written to geom.json
c     this line follows how iatph was set as data was read from geom.dat
            if (iphat(i).gt.nph) nph = iphat(i)
            if ( iatph(iphat(i)).eq.0) iatph(iphat(i)) = i
 20      continue
      end if
      call json%destroy()
         
      return
      end
