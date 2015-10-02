!=======================================================================
!     GLOBAL
!=======================================================================

module global_inp
  
  implicit none

  integer,private ::  ihuge
  double precision,private :: ten, huge, tiny, one
  parameter(ihuge = 10)
  parameter(one=1.d0, ten = 10.d0)
  parameter(huge = ten**ihuge, tiny = one/huge)

  !	the variables evnorm, xivnorm, spvnorm and l2lp are exclusive to nrixs (feffq) calculations
  !	le2 has different meaning for feffq calculations
  !	xivec serves many different functions depending on spectroscopy : xas/eels/nrixs
  integer do_nrixs,lj,ldecmx !KJ 7-09 for feff8q
  !	configuration average data :
  integer nabs, iphabs
  real*8 rclabs
  !	global polarization data :
  integer ipol, ispin, le2, l2lp
  real*8 elpty, angks
  real*8 evec(3), xivec(3), spvec(3)
  complex*16 ptz(-1:1,-1:1)
  double precision evnorm, xivnorm, spvnorm
  !   moved here because I think it belongs here !KJ 7-09
  integer ispec
  character(*),parameter,private :: filename='global.inp'  !KJ used to be global.dat !!!
  !       How many q-vectors:  (impulse transfer)
  integer nq
  !       Are we doing direction averaged impulse transfer?  (Note: this means q || e_z, which is not really averaging!)		
  logical qaverage
  !       The list of q-vectors and their norm:
  real*8,allocatable :: qs(:,:),qn(:)
  !       Weights of q-vectors in the cross-section (probably calculated by another code):
  complex*16,allocatable :: qw(:)
  !       Are we doing q,q' crossterms in the NRIXS code?
  logical mixdff
  !       If entering q, q' as length(q), length(q'), angle(q,q'), this is cosine(angle(q,q'))
  real*8,allocatable :: cosmdff(:,:)
  !       and this is norm(q')
  real*8 qqmdff
  !       A rotation matrix for each q-vector  (containing cos(theta),sin(theta),cos(fi),sin(fi) for each q)
  real*8,allocatable :: qtrig(:,:)  !compare to Adam's code qtrig(iq,1)=qcst(iq); 2)=qsnt; 3)=qcsf; 4)=qsnf
  !       Should the mdff program run? !11-2010
  integer imdff

contains

  subroutine init_feffq
    !	called to calculate some variables for nrixs
    integer i
    evnorm=0.0d0
    xivnorm=0.0d0
    spvnorm=0.0d0
    do i=1,3
       evnorm=evnorm+evec(i)*evec(i)
       xivnorm=xivnorm+xivec(i)*xivec(i)
       spvnorm=spvnorm+spvec(i)*spvec(i)
    end do
    spvnorm=sqrt(spvnorm)
    xivnorm=sqrt(xivnorm)
    evnorm=sqrt(evnorm)
  end subroutine init_feffq

  subroutine global_write(iniq)
    integer i
    logical,intent(in) :: iniq
    if(iniq) call init_feffq
    open (file=filename, unit=3, status='unknown')
    write (3, 10) ' nabs, iphabs - CFAVERAGE data'
    write (3, 45) nabs, iphabs, rclabs
45  format ( 2i8, f13.5)
    write (3,10) ' ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj' !KJ last 4 added for feff8q.  Note le2 new meaning in feff8q
    write (3, 50)  ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj !KJ
50  format ( 3i5, 2f12.4, 10i5)  !KJ
    write (3, 10) 'evec                   xivec            spvec'
    do 60 i = 1,3
       write (3,30) evec(i), xivec(i), spvec(i)
60     continue
       write (3, 10) ' polarization tensor '
       do 70 i = -1, 1
          write(3,30) dble(ptz(-1,i)), dimag(ptz(-1,i)), dble(ptz(0,i)), dimag(ptz(0,i)),  dble(ptz(1,i)), dimag(ptz(1,i))
70        continue
          !KJ for feff8q - was in different place in file in feff8q:
          write(3,10) 'evnorm, xivnorm, spvnorm - only used for nrixs'
          write (3,30) evnorm, xivnorm, spvnorm !KJ
          !KJ for a list of q-vectors (NRIXS) and MDFF calculation (NRIXS) - only relevant for NRIXS calculations  12-2010
          write(3,10) "nq,    imdff,   qaverage,   mixdff"
          write(3,*) nq,imdff,qaverage,mixdff
          write(3,*) 'q-vectors : qx, qy, qz, q(norm), weight, qcosth, qsinth, qcosfi, qsinfi'
          if(nq.gt.0) then    !note that this is redundant with xivec if nq=1, but ah well.
             do i=1,nq
                write(3,30) qs(i,:),qn(i),qw(i),qtrig(i,1:4)
             enddo
          endif
          if(mixdff) then
             write(3,*) "   qqmdff,   cos<q,q'>"
             write(3,*) qqmdff,cosmdff
          endif
          close(3)
          ! standard formats for string, integers and real numbers
10        format(a)
!20        format (20i4)
30        format (20f13.5)

        end subroutine global_write

        subroutine global_read
          real*8 aa1,bb1,aa2,bb2,aa3,bb3
          integer i
          open (file=filename, unit=3, status='old')
          read  (3,*)
          read  (3,*) nabs, iphabs, rclabs
          read  (3,*)
          read  (3,*)  ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj
          read  (3,*)
          do i = 1,3
             read  (3,*) evec(i), xivec(i), spvec(i)
          enddo
          read  (3,*)
          do i = -1, 1
             read (3,*) aa1, bb1, aa2, bb2, aa3, bb3 !KJ changed names of dummies to avoid confusion with my WELL DEFINED arrays (f*cking "implicit" people !#$&%)
             ptz(-1,i)= dcmplx(aa1,bb1)  !KJ changed cmplx to dcmplx to satisfy thorough compilers
             ptz(0,i) = dcmplx(aa2,bb2)
             ptz(1,i) = dcmplx(aa3,bb3)
          enddo
          read (3,*)
          read (3,*) evnorm, xivnorm, spvnorm !KJ
          if(do_nrixs .ne. 0) then !compatibility with (most) old files
             read(3,*)
             read(3,*) nq,imdff,qaverage,mixdff
             read(3,*)
             call make_qlist(nq)
             if(nq.gt.0) then    !note that this is redundant with xivec if nq=1, but ah well.
                do i=1,nq
                   read(3,30) qs(i,:),qn(i),qw(i),qtrig(i,1:4)
                enddo
             endif
             if(mixdff) then
                read(3,*)
                read(3,*) qqmdff,cosmdff
             endif
30           format (20f13.5)
          endif
          close(3)
        end subroutine global_read

        subroutine make_qlist(n)
          implicit none
          integer,intent(in) :: n
          allocate(qs(n,3),qn(n),qw(n),qtrig(n,4),cosmdff(n,n))
          qs(:,:)=0.d0
          qn(:)=0.d0
          qw(:)=1.d0
          qtrig(:,:)=0.d0
          qtrig(:,1)=1.d0
          qtrig(:,3)=1.d0 !corresponding to not rotating at all
          cosmdff(:,:)=1.d0
          return
        end subroutine make_qlist

        subroutine make_qtrig
          !          simple routine to get rotation angles; copied from mkptz      
          implicit none
          double precision rr,rsp
          integer iq
          do iq=1,nq
             if (qn(iq).gt.0.0d0) then
                rsp = qn(iq)
                rr = qs(iq,1)**2 + qs(iq,2)**2
                if (rr.lt. tiny) then
                   qtrig(iq,1) = - 1.d0
                   qtrig(iq,2) = 0.d0
                   qtrig(iq,3) = 1.d0
                   qtrig(iq,4) = 0.d0
                elseif (qs(iq,3).lt.0) then !meaning forward scattering ??
                   !                  rotation is defined by angles theta and fi
                   rr = sqrt(rr)
                   qtrig(iq,1) = qs(iq,3) / rsp
                   qtrig(iq,2) = rr / rsp
                   qtrig(iq,3) = qs(iq,1) / rr
                   qtrig(iq,4) = qs(iq,2) / rr
                else
                   qtrig(iq,1)=1.0d0
                   qtrig(iq,2)=0.0d0
                   qtrig(iq,3)=1.0d0 !surely this is a bug??  Shouldn't this be 1? !KJ 12-2011 changed 0->1 because produces NaN in genfmt otherwise
                   qtrig(iq,4)=0.0d0
                end if
             else
                call wlog(' FATAL error: one of the q-vectors is zero')
                call par_stop(' ')
             endif
          enddo !iq
          return
        end subroutine make_qtrig

        subroutine global_init
          ispec = 0
          ldecmx=-1 ! initialize the number of decomposition channels - KJ 7-09 for feff8q
          nabs = 1
          iphabs = 0
          rclabs = 0.d0
          ipol = 0
          ispin = 0
          le2 = 0
          l2lp = 0
          elpty = 0.d0
          angks = 0.d0
          evec(:) = 0.d0
          xivec(:) = 0.d0
          spvec(:) = 0.d0
          ptz(:,:) = dcmplx(0.d0,0.d0)
          evnorm=0.0d0
          xivnorm=0.0d0
          spvnorm=0.0d0
          do_nrixs=0 ! no nrixs calculation
          lj = -1
          nq=0
          qaverage=.true.
          imdff=0 !no mdff
          mixdff=.false. !no mdff in NRIXS
          !	cosmdff=1.d0  ! q || q'  => cos(0)=1     !KJ 11-2011 this is now an allocatable array.  Instruction fails on gfortran.
          qqmdff=-1.d0 ! leads to q=q' (norm only)
        end subroutine global_init

      end module global_inp
