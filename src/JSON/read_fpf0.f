      subroutine json_read_fpf0(nosc, oscstr, enosc)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision,dimension(:),allocatable :: dbpcs, dbpce

      integer nosc
      double precision oscstr(14), enosc(14)

      call json%load_file('fpf0.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read fpf0.json"
         stop
      else
         call json%get('nosc', nat, found)
                   if (.not. found) call bailout('nosc', 'fpf0.json')
         call json%get('oscstr', dbpcs, found)
                   if (.not. found) call bailout('x', 'atoms.json')
         call json%get('enosc',  dbpce, found)
                   if (.not. found) call bailout('y', 'atoms.json')
         do 20 i=1,nosc
            oscstr(i) = dbpcs(i)
            enosc(i)  = dbpce(i)
 20      continue
         call json%destroy()
      end if

      return
      end
