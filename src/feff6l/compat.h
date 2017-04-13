c{compat.h  -*-fortran-*-
c  used so setsca/getsca works with libxafs for autobk,diffkk,feff6, etc.
c  values for scalars
       integer          mxsca
       parameter       (mxsca=512)
       double precision sscaval(mxsca)
       character*96     sscanam(mxsca), snamtmp
       common /sca_dp/  sscaval
       common /sca_ch/  sscanam, snamtmp
c}
