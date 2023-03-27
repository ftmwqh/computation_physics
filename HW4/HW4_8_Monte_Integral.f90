program Monte_Carlo_Integral
   real*8, allocatable ::x(:), y(:), z(:), u(:), v(:), array(:)
   integer n, i, j, count
   real*8 epsilon, a, b, c, d, e, M, sum
   i = 1
   n = 10**7
   b = 5.
   a = 0.
   M = f(b)
   allocate (x(n))
   allocate (y(n))
   allocate (z(n))
   allocate (u(n))
   allocate (v(n))
   allocate (array(5*n))
   array = random_array(5*n)

   !掷石法求解第一问
   x = array(1:n)*(b - a) + a
   y = array(n + 1:2*n)*M

   count = 0
   do i = 1, n
      if (y(i) < f(x(i))) count = count + 1
   end do
   print*,"1."
   print *, "Monte Carlo (n=", (n), ")", ((b - a)*M*real(count)/n)
   print *, "standard", "     15.439010735567484"

   !平均值法求解第二问
   a = 0.7
   b = 4./7
   c = 0.9
   d = 2.
   e = 13./11.
   x = array(1:n)*a
   y = array(n + 1:2*n)*b
   z = array(2*n + 1:3*n)*c
   u = array(3*n + 1:4*n)*d
   v = array(4*n + 1:5*n)*e
   sum = 0
   do i = 1, n
      sum = sum + g(x(i), y(i), z(i), u(i), v(i))
   end do
   sum = sum/n*a*b*c*d*e
   print*,"2."
   print *, "Monte Carlo (n=", (n), ")", (sum)
   print *, "standard", "   5.67712092"
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
      real*8, allocatable :: array(:)
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

   real*8 function f(x)!第一问被积函数
      real*8 x
      f = sqrt(x**2 + 2.*sqrt(x))
      return
   end function

   real*8 function g(x, y, z, u, v)!第二问被积函数
      real*8 x, y, z, u, v
      g = 5 + x**2 - y**2 + 3*x*y - z**2 + u**3 - v**3
      return
   end function

end program Monte_Carlo_Integral
