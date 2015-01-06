      subroutine json_read_global(nabs, iphabs, rclabs, ipol, ispin,
     1                            le2, elpty, angks, evec, xivec,
     2                            spvec, ptz)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision,dimension(:),allocatable :: dbpcs

      integer nabs, iphabs, ipol, ispin, le2
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)

      call json%load_file('global.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read global.json"
         stop
      else
         call json%get('nabs', nabs, found)
                 if (.not. found) call bailout('nabs',   'global.json')
         call json%get('iphabs', iphabs, found)
                 if (.not. found) call bailout('iphabs', 'global.json')
         call json%get('rclabs', rclabs, found)
                 if (.not. found) call bailout('rclabs', 'global.json')
         call json%get('ipol', ipol, found)
                 if (.not. found) call bailout('ipol',   'global.json')
         call json%get('ispin', ispin, found)
                 if (.not. found) call bailout('ispin',  'global.json')
         call json%get('le2', le2, found)
                 if (.not. found) call bailout('le2',    'global.json')
         call json%get('elpty', elpty, found)
                 if (.not. found) call bailout('elpty',  'global.json')
         call json%get('angks', angks, found)
                 if (.not. found) call bailout('angks',  'global.json')

         call json%get('evec',  dbpcs, found)
                 if (.not. found) call bailout('evec',   'global.json')
         do 10 i=1,3
            evec(i) = dbpcs(i)
 10      continue
         call json%get('xivec',  dbpcs, found)
                 if (.not. found) call bailout('xivec',  'global.json')
         do 20 i=1,3
            xivec(i) = dbpcs(i)
 20      continue
         call json%get('spvec',  dbpcs, found)
                 if (.not. found) call bailout('spvec',  'global.json')
         do 30 i=1,3
            spvec(i) = dbpcs(i)
 30      continue

         call json%get('ptz0',  dbpcs, found)
                 if (.not. found) call bailout('ptz0',  'global.json')
         ptz(-1,-1) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0,-1) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1,-1) = dcmplx(dbpcs(5), dbpcs(6))
         call json%get('ptz1',  dbpcs, found)
                 if (.not. found) call bailout('ptz1',  'global.json')
         ptz(-1, 0) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0, 0) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1, 0) = dcmplx(dbpcs(5), dbpcs(6))
         call json%get('ptz2',  dbpcs, found)
                 if (.not. found) call bailout('ptz2',  'global.json')
         ptz(-1, 1) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0, 1) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1, 1) = dcmplx(dbpcs(5), dbpcs(6))
         call json%destroy()
      end if

      return
      end
