c  **************************************************
c  Parallel feff8 routines
c  Jim Sims
c  **************************************************

      subroutine par_begin
c  **************************************************
c  Initializations for parallel version(s)
c  **************************************************

      include '../HEADERS/parallel.h'
      include 'mpif.h'

c-- So cvd or dbx can attach to a running process
c     call sleep(30) 

      call MPI_INIT(ierrorflag)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierrorflag)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierrorflag)
      this_process = my_rank

      par_type = 0
      parallel_run = .true.
c-- The following variable will be used for IO that should only be
c-- done in one process.
      master = (my_rank .eq. 0)

      worker = (.not. master)
      if (worker) par_type = 1

c     write(6,*) 'this process = ',this_process, ' worker = ',worker

      if (master) write(6,*) 'Number of processors = ',numprocs

      return
      end

      subroutine par_stop (string)
c  **************************************************
c  Abnormal termination of the parallel session
c  **************************************************
      include '../HEADERS/parallel.h'
      include 'mpif.h'
c     For abnormal exits 
c     If open, close unit = 11
c     Go to the barrier that workers are sitting at
c     Then everyone will call par_end and stop
      logical is_open
      character*(*) string

      inquire(unit=11,opened=is_open)
      if (is_open) then
        call wlog(string)
        close(unit=11)
      else if (string .ne. ' ') then
	print *,string
	print *,'Abnormal termination on processor ',this_process
      endif
      call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierrorflag)

      stop ' '
      end

      subroutine par_end
c  **************************************************
c  Terminate the parallel session
c  **************************************************
      call MPI_FINALIZE(ierrorflag)
      return
      end

      subroutine par_barrier
c  **************************************************
c  Calls mpi_barrier
c  **************************************************
      include 'mpif.h'
      call MPI_BARRIER(MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_int(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for integer arrays
c  **************************************************
      integer count,dest,tag
      integer buf(*)
      include 'mpif.h'
      call MPI_SEND(buf,count,MPI_INTEGER,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_cmplx(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for complex arrays
c  **************************************************
      integer count,dest,tag
      complex buf(*)
      include 'mpif.h'
      call MPI_SEND(buf,count,MPI_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_dc(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for double_complex arrays
c  **************************************************
      integer count,dest,tag
      complex*16 buf(*)
      include 'mpif.h'
      call MPI_SEND(buf,count,MPI_DOUBLE_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_recv_int(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for integer arrays
c  **************************************************
      integer count,source,tag
      integer buf(*)
      include 'mpif.h'
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_INTEGER,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_cmplx(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for complex arrays
c  **************************************************
      integer count,source,tag
      complex buf(*)
      include 'mpif.h'
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_dc(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for double complex arrays
c  **************************************************
      integer count,source,tag
      complex*16 buf(*)
      include 'mpif.h'
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_DOUBLE_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_bcast_int(buf,count,source)
c  **************************************************
c  Call mpi_bcast for integer arrays
c  **************************************************
      integer count,source
      integer buf(*)
      include 'mpif.h'
      call MPI_BCAST(buf,count,MPI_INTEGER,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_cmplx(buf,count,source)
c  **************************************************
c  Call mpi_bcast for complex arrays
c  **************************************************
      integer count,source
      complex buf(*)
      include 'mpif.h'
      call MPI_BCAST(buf,count,MPI_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_dc(buf,count,source)
c  **************************************************
c  Call mpi_bcast for double_complex arrays
c  **************************************************
      integer count,source
      complex*16 buf(*)
      include 'mpif.h'
      call MPI_BCAST(buf,count,MPI_DOUBLE_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine MPE_DECOMP1D( n, num_procs, myid, s, e )
c  ******************************************************
c  A routine for producing a decomposition of a 1-d 
c  array when given a number of processors.  It may 
c  be used in "direct" product decomposition.  The 
c  values returned assume a "global" domain in [1:n]
c  ******************************************************
c  MPE_Decomp1d - Compute a balanced decomposition of
c  a 1-D array
c  ******************************************************
c  Input Parameters:
c  n  - Length of the array
c  num_procs - Number of processors in decomposition
c  myid  - Rank of this processor in the decomposition 
c  (0 <= rank < size)
c  ******************************************************
c  Output Parameters:
c  s,e - Array my_particles are s:e, with the original 
c  array considered as 1:n.  
c  ******************************************************

      integer n, num_procs, myid, s, e
      integer nloc, deficit
 
      nloc  = n / num_procs
      s       = myid * nloc + 1
      deficit = mod(n,num_procs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
        nloc = nloc + 1
      endif
      e = s + nloc - 1
      if (e .gt. n .or. myid .eq. num_procs-1) e = n

      return
      end

      SUBROUTINE SECONDS( W)
c  ***************************************************
c  SECONDS returns the wall clock times for a process
c  in seconds.
c  ***************************************************

      real*8 W, MPI_Wtime

      W = MPI_Wtime()

      RETURN
      END
