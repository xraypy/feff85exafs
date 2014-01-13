c====================================================================
      integer function nxtunt(iunit)

c  this function returns the value of the next unopened unit number
c  equal to or larger than iunit.  it will return neither unit numbers
c  0, 5, or 6 nor a negative unit number
c $Id: nxtunt.f,v 1.1.1.1 2006/01/12 06:37:42 hebhop Exp $
c $Log: nxtunt.f,v $
c Revision 1.1.1.1  2006/01/12 06:37:42  hebhop
c New version of feff. feff8.5 (Extension of feff8.4)
c Includes:
c 	1) All feff8.4 capabilities.
c 	2) Screened core hole (calculation of W).
c 	3) Multiple pole self energy calculation.
c 	4) Convolution with spectral function.
c New cards and options:
c 	1) NOHOLE 2      (screened hole)
c 	2) PLASMON ipl   (multiple pole self energy)
c 	3) SO2CONV       (convolve output with spectral function)
c 	4) SELF          (print on shell self energy as calculated by Luke)
c 	5) SFSE k0        (print off shell self energy Sigma(k0,e) )
c
c Revision 1.1.1.1  2000/02/11 02:23:58  alex
c Initialize feff82
c
c Revision 1.10  1999/04/02 21:32:47  newville
c cleaned up nxtunt (matt)
c
c Revision 1.9  1999/02/11 20:08:08  alex
c x39 version: dim.h + misc. small changes
c
c Revision 1.8  1998/12/29 23:59:07  alex
c feff8x35 version
c
c Revision 1.7  1998/11/19 03:23:11  alex
c feff8x32 version
c
c Revision 1.6  1998/10/26 14:11:16  ravel
c no comments beyond column 71
c
c Revision 1.5  1998/10/18 21:47:51  alex
c feff8x30 version implements Broyden algorithm for self-consistency
c
c Revision 1.4  1998/02/24 18:31:37  ravel
c I should really be more careful.  This is the last commitment done
c      cright.
c
c Revision 1.1.1.1  1997/04/27 20:18:03  ravel
c Initial import of xanes sources, version 0.37
c
c Revision 1.1  1996/06/23 16:05:02  bruce
c Initial revision
c

       integer iunit
       logical open

       nxtunt = max(1, iunit) - 1
 10    continue
       nxtunt = nxtunt + 1
       if ((nxtunt.eq.5).or.(nxtunt.eq.6)) nxtunt = 7
       inquire (unit=nxtunt, opened=open)
       if (open) go to 10
       return
c  end integer function nxtunt
       end

c====================================================================
