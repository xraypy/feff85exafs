c -*-fortran-*-
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
