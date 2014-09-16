      subroutine read_xsect(ntit, titles, s02x, erelax, wpx, edge,
     1                      emu, gamma, ne, ne1, ik0,
     2                      er, ei, xsn, col4, col5)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision,dimension(:),allocatable :: dbpc1, dbpc2, dbpc3
      double precision,dimension(:),allocatable :: dbpc4, dbpc5
      character*80,dimension(:),allocatable :: strings


      character*80 titles(nheadx)
      double precision er(nex), ei(nex), xsn(nex)
      double precision col4(nex), col5(nex)

      call json%load_file('xsect.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read xsect.json"
         stop
      else
         call json%get('ntitle', ntit, found)
                   if (.not. found) call bailout('s02',   'xsect.json')
         call json%get('title',   strings, found)
                   if (.not. found) call bailout('title', 'xsect.json')
         do 10 i=1,ntit
            titles(i) = strings(i)
 10      continue
         

         call json%get('s02',    s02x, found)
                   if (.not. found) call bailout('s02',    'xsect.json')
         call json%get('erelax', erelax, found)
                   if (.not. found) call bailout('erelax', 'xsect.json')
         call json%get('wp',     wpx, found)
                   if (.not. found) call bailout('wp',     'xsect.json')
         call json%get('edge',   edge, found)
                   if (.not. found) call bailout('edge',   'xsect.json')
         call json%get('emu',    emu, found)
                   if (.not. found) call bailout('emu',    'xsect.json')
         call json%get('gamach', gamma, found)
                   if (.not. found) call bailout('gamach', 'xsect.json')
         call json%get('ne',     ne, found)
                   if (.not. found) call bailout('ne',     'xsect.json')
         call json%get('ne1',    ne1, found)
                   if (.not. found) call bailout('ne1',    'xsect.json')
         call json%get('ik0',    ik0, found)
                   if (.not. found) call bailout('ik0',    'xsect.json')

         call json%get('ereal',   dbpc1, found)
                   if (.not. found) call bailout('er',     'xsect.json')
         call json%get('eimag',   dbpc2, found)
                   if (.not. found) call bailout('ei',     'xsect.json')
         call json%get('xsnorm',  dbpc3, found)
                   if (.not. found) call bailout('xsn',    'xsect.json')
         call json%get('dum1',   dbpc4, found)
                   if (.not. found) call bailout('col4',   'xsect.json')
         call json%get('dum2',   dbpc5, found)
                   if (.not. found) call bailout('col5',   'xsect.json')

         do 20 i=1,ne
            er(i)   = dbpc1(i)
            ei(i) = dbpc2(i)
            xsn(i) = dbpc3(i)
            col4(i) = dbpc4(i)
            col5(i) = dbpc5(i)
 20      continue
         call json%destroy()
      end if

      return
      end

      subroutine read_titles(ntit, titles)
      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      character*80,dimension(:),allocatable :: strings

      character*80 titles(nheadx)

      call json%load_file('xsect.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read xsect.json"
         stop
      else
         call json%get('ntitle', ntit, found)
                   if (.not. found) call bailout('s02',   'xsect.json')
         call json%get('title',   strings, found)
                   if (.not. found) call bailout('title', 'xsect.json')
         do 10 i=1,ntit
            titles(i) = strings(i)
 10      continue
         call json%destroy()
      end if

      return
      end
