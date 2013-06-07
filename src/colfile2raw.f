      subroutine f8px( num, str)
      real*8 num
      integer isigma
      character *8 str
      character *9 arse
! this still gets the wrong decimal places on sigma
! obviously it misses some trick
!     8 places to use 
!               12345678
      if( abs(num).lt.999.999) then
         write(str, '(F8.3)') num
         return
      endif
      if (abs(num).lt.9999.99) then
        write(str, '(F8.2)') num
        return
      endif
!               12345678
      if (num.lt.99999.9) then
        write(str, '(F8.1)') num
        return
      endif
!               12345678
      if (num.lt.999999.) then
        write(str, '(F8.0)') num
        return
      endif
      write(arse, '(ES9.3E1)') num
      str(1:1) = ' '
      str(2:6) = arse(2:6)
      str(7:7) = 'e'
      str(8:8) = arse(9:9)
      return
      end subroutine f8px
!
!
!
      program colfile2raw
! 38 values on a line, all are numbers
      real*8 d(38)
      integer*8 i(38),j,ikey
      character *8 iobs, sigi
! titles
      read(*,*)
 1    read(*,*,end=100, err=99) d
! rounding, whatever
      i = d
      call f8px( d(4), iobs, 0 )
      call f8px( d(5), sigi, 1 )
      if (iobs(4:4).eq.'*') goto 99
      ikey = 1048576*(i(1)+511) + 1024*(i(2)+511) + i(3) + 511
      write(*,1000)i(1:3),iobs,sigi,i(6),d(7:12),i(13),d(14:24),
     &i(25:28),d(29),i(30),d(31),i(32),i(33),d(34:38),1
 1000 format( 3I4,A8,A8,I4,6F8.5,I3,F7.2,F7.2,F8.2,F7.2,F7.2,F8.2,
     &  F6.3,F5.2,F7.2,F7.2,F7.2,I2,I5,I9,I7,F7.2,I4,F6.3,I11,I3,F6.3,
     & F8.2,F8.2,F8.3,F8.3,I4)
      goto 1
 99   write(*,*)'got an error'
      write(*,*)d(4),d(5)
 100  stop
      end program colfile2raw
      
