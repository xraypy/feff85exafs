      subroutine json_nnnn(fjson, ntit, titles, rat, ipot, ri, beta,eta,
     &          index, iorder, nleg, deg, reff, rnrmav, edge,
     &          ne, col1, col2, col3, col4, col5, col6, col7)

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      character*13 fjson
      character*5  vname
      character*80 titles(nheadx)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      double precision rat(3,0:legtot+1), atom(4)
      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      double precision ria(legtot), betad(legtot+1), etad(0:legtot+1)
      integer ipot(0:legtot)

      integer  iunit
      type(json_value),pointer :: nnnnj
      call json_value_create(nnnnj)
      call to_object(nnnnj,fjson)

      call json_value_add(nnnnj, 'vfeff',   vfeff)
      call json_value_add(nnnnj, 'vf85e',   vf85e)

      call json_value_add(nnnnj, 'titles',  titles(1:ntit))

      call json_value_add(nnnnj, 'index',   index)
      call json_value_add(nnnnj, 'iorder',  iorder)
      call json_value_add(nnnnj, 'nleg',    nleg)
      call json_value_add(nnnnj, 'degen',   deg)
      call json_value_add(nnnnj, 'reff',    reff*bohr)
      call json_value_add(nnnnj, 'rnorman', rnrmav)
      call json_value_add(nnnnj, 'edge',    edge*hart)
c     call json_value_add(nnnnj, 'ne',      ne)

      do 10 i=1,nleg
         ria(i)   = ri(i)*bohr
         betad(i) = beta(i)*180/pi 
         etad(i)  = eta(i)*180/pi 
 10   continue
      call json_value_add(nnnnj, 'ri',       ria(1:nleg))
      call json_value_add(nnnnj, 'beta',     betad(1:nleg))
      call json_value_add(nnnnj, 'eta',      etad(1:nleg))
c     call json_value_add(nnnnj, 'ipot',     ipot(0:nleg))

      do 20 iat=1,nleg
         write (vname, "(A4,I1)") "atom", iat
         do 30 ix=1,3
            atom(ix) = rat(ix, iat)*bohr
 30      continue
         atom(4) = dble(ipot(iat))
         call json_value_add(nnnnj, vname, atom)
 20   continue

      call json_value_add(nnnnj, 'k',        col1(1:ne))
      call json_value_add(nnnnj, 'real_phc', col2(1:ne))
      call json_value_add(nnnnj, 'mag_feff', col3(1:ne))
      call json_value_add(nnnnj, 'pha_feff', col4(1:ne))
      call json_value_add(nnnnj, 'red_fact', col5(1:ne))
      call json_value_add(nnnnj, 'lam',      col6(1:ne))
      call json_value_add(nnnnj, 'rep',      col7(1:ne))

      open(newunit=iunit, file=fjson, status='REPLACE')
      call json_print(nnnnj,iunit)
      close(iunit)
      call json_destroy(nnnnj)

      return
      end
