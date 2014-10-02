      subroutine json_nnnn(fjson, ntit, titles, rat, ipot, ri, beta,eta,
     &          index, iorder, nleg, deg, reff, rnrmav, edge,
     &          ne, col1, col2, col3, col4, col5, col6, col7)
c+----------------------------------------------------------------------
c  Write a JSON file conatining all the information of a traditional
c  feffNNNN.dat file.
c
c  The keys for the entries in the JSON file are chosen to be the same 
c  as the corresponding attributes in a Larch FeffPath or _feffdat group
c  See http://xraypy.github.io/xraylarch/xafs/feffpaths.html#feffpath-and-feffpath-groups
c
c  Much of the header of a traditional feffNNNN.dat file is stored as
c  an array of text strings.  Degeneracy, nleg, and Reff have their
c  own keys in the JSON.
c
c  The leg lengths and beta & eta angles are also included in the JSON
c
c  Also included are the atom positions of the atoms in the path.  The
c  keys are "atomN" where N is an integer from 1 to nleg.  Each is a 
c  4-element array containing x,y,z,ipot of the atom.  Note that the
c  ipot is stored as a double and will have to be converted back to an
c  integer when imported.
c
c  fjson contains the path/name of the output JSON file
c
c  the colN inputs are the columns of the traditional feffNNNN.dat file.
c  they are computed in GENFMT/fdtarr.f
c+----------------------------------------------------------------------

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      parameter (eps = 1.0D-8)

      character*13 fjson
      character*5  vname
      character*80 titles(nheadx)
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      double precision rat(3,0:legtot+1), atom(4)
      double precision ri(legtot),  beta(legtot+1),  eta(0:legtot+1)
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
         if (abs(etad(i)-360) .lt. eps) then
            etad(i) = 0   
         end if
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
