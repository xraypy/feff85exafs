c{echo.h -*-fortran-*-
c i_echo controls echo outputs:
c   0   save to echo buffer  in echo_str
c   1   write to screen 
c   2   write to echo file
c   3   write to both screen and echo file
       integer  mxecho, n_echo, i_echo, lun_echo
       parameter(mxecho =  512)
       character*512    echo_str(mxecho), echo_file
       common /echo_s/  echo_str, echo_file
       common /echo_i/  n_echo, i_echo, lun_echo
c}
