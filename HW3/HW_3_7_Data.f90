program Data
!取分段函数作为舍选法的F(x)
   real*8, allocatable ::xi_x(:), xi_y(:), xi_1(:), xi_2(:), x(:), e(:), y(:), array(:)
   integer n, i, j
   real*8 a, b, sum
   j = 1
   n = 10**6
   allocate (xi_1(n))
   allocate (xi_2(n))
   allocate (xi_x(n))
   allocate (xi_y(n))
   allocate (e(114))
   allocate (y(114))
   allocate (x(n))
   allocate (array(2*n))
   array = random_array(2*n)!array前n用于存储xi_1的随机数，后n存储xi_2的随机数
   !舍选法
   xi_1 = array(1:n)
   xi_2 = array(n + 1:2*n)

   open (unit=1, file='data.TXT')!从data.TXT文件读取
   read (1, *)!跳过第一行的中文
   do i = 1, 114
      read (1, *) e(i), y(i)
   end do
   y = y/sum_array(y, 114)

   do i = 1, n!计算0-1的均匀随机数xi_1对应的xi_x
      if (xi_1(i) >= 0 .and. xi_1(i) <= 93./181) then
         xi_x(i) = 2.715*xi_1(i)/0.015 + 2900.
      else if (xi_1(i) > 93./181 .and. xi_1(i) < 173./181) then
         xi_x(i) = (2.715*xi_1(i) - 1.395)/0.1 + 2993.
      else if (xi_1(i) >= 173./181 .and. xi_1(i) <= 1.) then
         xi_x(i) = (2.715*xi_1(i) - 2.595)/0.015 + 3005.
      end if

      xi_y(i) = xi_2(i)*F(xi_x(i))!xi_2对应的xi_y

      if (xi_y(i) < y(int(xi_x(i) + 0.5) - 2899)) then!选取满足舍选法要求的计入数组x
         x(j) = int(xi_x(i) + 0.5)
         j = j + 1
      end if
   end do

   open (unit=2, file='x_1.txt')!x写入文件
   do i = 1, j - 1
      write (2, *) x(i)
   end do
   print *, "rate:", (real(j)/n) !计算抽样率

   !直接抽样法
   open (unit=3, file='x_2.txt')
   do i = 1, n
      sum = 0
      j = 1
      do while (xi_1(i) > sum)!求取均匀随机变量xi_1对应的积分到哪个点
         sum = sum + y(j)
         j = j + 1
      end do
      x(i) = e(j - 1)
      write (3, *) x(i)
   end do
   print *, "end"

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

   real*8 function F(x)!舍选法权重函数
      real*8 x, a, b
      a = 2993.
      b = 3005.
      if (x > 2993 .and. x < 3005) then
         F = 0.1
      else
         F = 0.015
      end if
      return
   end function

   real*8 function sum_array(x, n)!对一个数列求和
      real*8::x(n)
      integer n
      sum_array = 0
      do i = 1, n
         sum_array = sum_array + x(i)
      end do
      return
   end function
end program Data
