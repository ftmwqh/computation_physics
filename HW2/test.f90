program test
   real, allocatable :: u(:)

   u = random_array(3)
   print *, (u)
   read (*, *)

contains
  function random_array(n) result(array)
      implicit none
      integer :: n

      integer  z(3)
      integer i, m, a, r, q
      real :: array(n)

      m = 2**31 - 1
      a = 7**5
      r = mod(m, a)
      q = m/a
      z(1) = seed()                          !初始值由seed产生

      do i = 1, n
         z(i + 1) = a*mod(z(i), q) - r*(z(i)/q)
         if (z(i + 1) < 0) z(i + 1) = z(i + 1) + m
      end do
      array = real(z)/real(m)
   end function random_array

   integer function seed()                            !产生随机数种子值
      integer, dimension(8) :: time
      integer I_0
      call date_and_time(VALUES=time)             !获取系统时间
      I_0 = time(1) + 70*(time(2) + 12*(time(3) + 31*(time(5) + 23*(time(6) + 59*time(7)))))
      seed = I_0
      return
   end function

   
end program test