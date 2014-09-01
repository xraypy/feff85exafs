      subroutine json_read_atoms(nat, rat, iphat)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      integer,dimension(:),allocatable :: intgs
      double precision,dimension(:),allocatable :: dbpcs, dbpcy, dbpcz

      integer nat, iphat(natx)
      double precision rat(3,natx)

      call json%load_file('atoms.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read atoms.json"
         stop
      else
         call json%get('natt', nat, found)
                   if (.not. found) call bailout('natt', 'atoms.json')
         call json%get('x',    dbpcs, found)
                   if (.not. found) call bailout('x', 'atoms.json')
         call json%get('y',    dbpcy, found)
                   if (.not. found) call bailout('y', 'atoms.json')
         call json%get('z',    dbpcz, found)
                   if (.not. found) call bailout('z', 'atoms.json')
         call json%get('iphatx', intgs, found)
                   if (.not. found) call bailout('iphatx', 'atoms.json')
         do 20 i=1,nat
            iphat(i)  = intgs(i)
            rat(1,i)  = dbpcs(i)
            rat(2,i)  = dbpcy(i)
            rat(3,i)  = dbpcz(i)
 20      continue
         call json%destroy()
      end if

      return
      end
