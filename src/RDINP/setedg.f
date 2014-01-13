      subroutine setedg (a2, ihole)
      integer i, ihole
      character*2 a2, edglbl, edglbp
      dimension edglbl(0:29), edglbp(0:29)

      data edglbl / 'NO', 'K ', 'L1', 'L2', 'L3',
     3            'M1','M2','M3','M4','M5',
     4            'N1','N2','N3','N4','N5','N6','N7',
     5            'O1','O2','O3','O4','O5','O6','O7',
     6            'P1','P2','P3','P4','P5','R1' /
      data edglbp / '0', '1 ', '2', '3', '4',
     3            '5','6','7','8','9',
     4            '10','11','12','13','14','15','16',
     5            '17','18','19','20','21','22','23',
     6            '24','25','26','27','28','29' /

      ihole  = -1
      do 10 i = 0,29
  10     if (a2 .eq. edglbl(i) .or. a2 .eq. edglbp(i) ) ihole  = i
      if (ihole  .lt. 0) call par_stop('unknown EDGE')

      return
      end
