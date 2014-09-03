      subroutine bailout(msg, file)
      character*(*) msg, file
      print *, "Error reading "//file//", failed to find >"//msg//"<"
      stop
      end
