      subroutine json_read_onepath(index, iorder, ipol, nleg, deg, rat,
     &        ipot, elpty, evec, xivec, nnnn_out, json_out)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision,dimension(:),allocatable :: dbpcs, dbpcx

      integer index, nleg, ipot(0:legtot)
      double precision rat(3,0:legtot+1)
      logical nnnn_out, json_out
      double precision evec(3), xivec(3)

      character*5 vname

      elpty = 0
      do 5 ix=1,3
         evec(ix)  = 0
         xivec(ix) = 0
 5    continue


      call json%load_file('feffpath.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read feffpath.json"
         stop
      else
         call json%get('index', index, found)
             if (.not. found) call bailout('index', 'feffpath.json')
         call json%get('nleg',  nleg,  found)
             if (.not. found) call bailout('nleg',  'feffpath.json')
         call json%get('deg',   deg,   found)
             if (.not. found) call bailout('deg',  'feffpath.json')

         call json%get('iorder', iorder, found)
             if (.not. found) call bailout('iorder', 'feffpath.json')
         call json%get('ipol',   ipol,   found)
             if (.not. found) call bailout('ipol',   'feffpath.json')
         call json%get('ellipticity', elpty,  found)
           if (.not. found) call bailout('ellipticity', 'feffpath.json')


         call json%get('nnnn_out', nnnn_out, found)
             if (.not. found) call bailout('nnnn_out', 'feffpath.json')
         call json%get('json_out', json_out, found)
             if (.not. found) call bailout('json_out', 'feffpath.json')

         do 10 iat=1,nleg
            write (vname, "(A4,I1)") "atom", iat
            call json%get(vname, dbpcs, found)
                 if (.not. found) call bailout(vname, 'feffpath.json')
            
c           convert distances to code units
            rat(1,iat)  = dbpcs(1)/bohr
            rat(2,iat)  = dbpcs(2)/bohr
            rat(3,iat)  = dbpcs(3)/bohr
            ipot(iat)   = int(dbpcs(4)+0.5)
 10      continue

         call json%get('evec',  dbpcs, found)
              if (.not. found) call bailout('evec',  'feffpath.json')
         call json%get('xivec', dbpcx, found)
              if (.not. found) call bailout('xivec', 'feffpath.json')
         do 20 ix=1,3
            evec(ix)  = dbpcs(ix)
            xivec(ix) = dbpcx(ix)
 20      continue


         call json%destroy()
      end if


      return
      end


      subroutine json_read_geometry(nleg, rat, ipot)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      integer nleg, ipot(0:legtot)
      double precision rat(3,0:legtot+1)
      character*5 vname


      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
c      double precision,dimension(:),allocatable :: dbpcs

      call json%load_file('feffpath.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read feffpath.json"
         stop
      else
         call json%get('nleg',  nleg,  found)
             if (.not. found) call bailout('nleg',  'feffpath.json')
         do 10 iat=1,nleg
            write (vname, "(A4,I1)") "atom", iat
            call json%get(vname, dbpcs, found)
                 if (.not. found) call bailout(vname, 'feffpath.json')
            
c           convert distances to code units
            rat(1,iat)  = dbpcs(1)/bohr
            rat(2,iat)  = dbpcs(2)/bohr
            rat(3,iat)  = dbpcs(3)/bohr
            ipot(iat)   = int(dbpcs(4)+0.5)
 10      continue

         call json%destroy()
      end if


      return
      end
