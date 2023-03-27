program half_sphere
   !use random,only:random_array !此处试用其他文件module,为提交方便直接在文件内定义
   implicit none
   real pi
   real, allocatable :: x(:), y(:), z(:), r(:), u(:), v(:), t(:)
   integer n, i, count
   pi = 4*atan(1.0)!定义pi值
   n = 10**4
   t = 2*random_array(2*n) - 1
   allocate (u(n))
   allocate (v(n))
   allocate (x(n))
   allocate (y(n))
   allocate (z(n))
   do i = 1, n
      u(i) = t(2*i - 1)
      v(i) = t(2*i)
   end do

   r = (u**2 + v**2)**(0.5)!Marsaglia法抽样
   x = 2*u*sqrt(1 - r**2)
   y = 2*v*sqrt(1 - r**2)
   z = 1 - 2*r**2

   open (unit=1, file='xy.txt')!写入文件
   do i = 1, n
      valid3: if (r(i) <= 1) then

         write (1, "(F8.6XF8.6XF8.6)") x(i), y(i), z(i)!此为输出8位(7小数位)，中间空格X
      end if valid3
   end do

   read (*, *)

contains
   integer function seed()                            !产生随机数种子值
      integer, dimension(8) :: time
      integer I_0
      call date_and_time(VALUES=time)             !获取系统时间
      I_0 = time(1) + 70*(time(2) + 12*(time(3) + 31*(time(5) + 23*(time(6) + 59*time(7)))))
      seed = I_0
      return
   end function

   function random_array(n) result(array)!产生长度为n的[0,1]的随机数
      implicit none
      integer :: n

      integer, allocatable ::  z(:)
      integer i, m, a, r, q
      real, allocatable :: array(:)
      allocate (z(n))
      allocate (array(n))
      m = 2**31 - 1
      a = 7**5
      r = mod(m, a)
      q = m/a
      z(1) = seed()                          !初始值由seed产生

      do i = 1, n - 1                          !数组超上限会出现SIGSEGV或者自己断点的情况
         z(i + 1) = a*mod(z(i), q) - r*(z(i)/q)
         if (z(i + 1) < 0) z(i + 1) = z(i + 1) + m

      end do
      array = real(z)/real(m)

   end function

end program half_sphere
