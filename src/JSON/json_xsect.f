      subroutine json_xsect(ntit, titles, s02x, erelax, wpx, edge,
     1                      emu, gamma, ne, ne1, ik0,
     2                      er, ei, xsn, col4, col5)

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      character*80 titles(nheadx)
      double precision er(nex), ei(nex), xsn(nex)
      double precision col4(nex), col5(nex)

      integer  iunit
      type(json_value),pointer :: xsectj
      call json_value_create(xsectj)
      call to_object(xsectj,'xsect.json')

      call json_value_add(xsectj, 'vfeff',  vfeff)
      call json_value_add(xsectj, 'vf85e',  vf85e)

      call json_value_add(xsectj, 'ntitle', ntit)
      call json_value_add(xsectj, 'title',  titles(1:ntit))

      call json_value_add(xsectj, 's02',    s02x)
      call json_value_add(xsectj, 'erelax', erelax)
      call json_value_add(xsectj, 'wp',     wpx)
      call json_value_add(xsectj, 'edge',   edge)
      call json_value_add(xsectj, 'emu',    emu)

      call json_value_add(xsectj, 'gamach', gamma)
      call json_value_add(xsectj, 'ne',     ne)
      call json_value_add(xsectj, 'ne1',    ne1)
      call json_value_add(xsectj, 'ik0',    ik0)

      call json_value_add(xsectj, 'ereal',  er(1:ne))
      call json_value_add(xsectj, 'eimag',  ei(1:ne))
      call json_value_add(xsectj, 'xsnorm', xsn(1:ne))
      call json_value_add(xsectj, 'dum1',   col4(1:ne))
      call json_value_add(xsectj, 'dum2',   col5(1:ne))

      open(newunit=iunit, file='xsect.json', status='REPLACE')
      call json_print(xsectj,iunit)
      close(iunit)
      call json_destroy(xsectj)

      return
      end
