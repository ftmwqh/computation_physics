program Random_walk
   real*8, allocatable ::x(:, :, :), ave(:, :), ave_2(:), C(:), vx(:, :), vy(:, :), array(:), vx_0(:), vy_0(:)
   integer n, i, j, count, t
   real*8 A, w, b, d, e, M, sum, vx0, vy0, u, v, pi
   pi = 3.1415926536
   i = 1
   t = 100
   n = 10**5
   b = 5.
   w = 0.5
   A = 1.0!令粒子在沿电场运动的概率~1+A*sin(wt)

   allocate (x(n, t, 2))!记录每个（n个）粒子在每个（t个）时间的坐标(x,y)对应(1,2)
   allocate (ave(t, 2))
   allocate (ave_2(t))
   allocate (C(t - 1))
   allocate (vx(n, t - 1))
   allocate (vy(n, t - 1))
   allocate (vx_0(n))
   allocate (vy_0(n))
   allocate (array(t*n))
   array = random_array(t*n)
   particles: do i = 1, n!进行MonteCarlo模拟
      !初始化从初始时刻在0位置走一次
      if (array((i - 1)*100 + 1) > 0 .and. array((i - 1)*100 + 1) < 0.25*(1 + A*sin(w*1))) then
         x(i, 1, 1) = +1.0!x+方向
         x(i, 1, 2) = 0.0
      else if (array((i - 1)*100 + 1) > 0.25*(1 + A*sin(w*1)) .and. array((i - 1)*100 + 1) < 0.5) then
         x(i, 1, 1) = -1.0!x-方向
         x(i, 1, 2) = 0.0
      else if (array((i - 1)*100 + 1) > 0.5 .and. array((i - 1)*100 + 1) < 0.75) then
         x(i, 1, 1) = 0.0
         x(i, 1, 2) = 1.0!y+方向
      else if (array((i - 1)*100 + 1) > 0.75 .and. array((i - 1)*100 + 1) < 1.0) then
         x(i, 1, 1) = 0.0
         x(i, 1, 2) = -1.0!y-方向
      end if
      time: do j = 2, t
         !根据[0,1]均匀分布的随机数所在区间判断随机游走方向
         if (array((i - 1)*100 + j) > 0 .and. array((i - 1)*100 + j) < 0.25*(1 + A*sin(w*j))) then
            x(i, j, 1) = x(i, j - 1, 1) + 1.0!x+方向
            x(i, j, 2) = x(i, j - 1, 2)
            vx(i, j - 1) = 1.0!速度x+方向
            vy(i, j - 1) = 0.0
         else if (array((i - 1)*100 + j) > 0.25*(1 + A*sin(w*j)) .and. array((i - 1)*100 + j) < 0.5) then
            x(i, j, 1) = x(i, j - 1, 1) - 1.0!x-方向
            x(i, j, 2) = x(i, j - 1, 2)
            vx(i, j - 1) = -1.0!速度x-方向
            vy(i, j - 1) = 0.0
         else if (array((i - 1)*100 + j) > 0.5 .and. array((i - 1)*100 + j) < 0.75) then
            x(i, j, 1) = x(i, j - 1, 1)
            x(i, j, 2) = x(i, j - 1, 2) + 1!y+方向
            vx(i, j - 1) = 0.0
            vy(i, j - 1) = 1.0
         else if (array((i - 1)*100 + j) > 0.75 .and. array((i - 1)*100 + j) < 1.0) then
            x(i, j, 1) = x(i, j - 1, 1)
            x(i, j, 2) = x(i, j - 1, 2) - 1!y-方向
            vx(i, j - 1) = 0.0
            vy(i, j - 1) = -1.0
         end if
      end do time
   end do particles

   do j = 1, t!求n次模拟的平均值<x(n)>
      ave(j, 1) = 0
      ave(j, 2) = 0
      do i = 1, n
         ave(j, 1) = ave(j, 1) + x(i, j, 1)/real(n)
         ave(j, 2) = ave(j, 2) + x(i, j, 2)/real(n)
      end do
   end do
   open (unit=1, file='A1w05x.txt')!x写入文件,导入python画图坐标随时间变化
   open (unit=2, file='A1w05y.txt')
   do i = 1, t
      write (1, "(F8.4)") ave(i, 1)!此为输出t时间段内x,y坐标
      write (2, "(F8.4)") ave(i, 2)
   end do

   do j = 1, t!求n次模拟的平均值<x^2(n)>
      ave_2(j) = 0
      do i = 1, n
         ave_2(j) = ave_2(j) + (x(i, j, 1)**2 + x(i, j, 2)**2)/real(n)
      end do
   end do
   open (unit=3, file='A1w05x2.txt')!x^2写入文件,导入python画图坐标随时间变化
   do i = 1, t
      write (3, "(F8.4)") ave_2(i)!此为输出t时间段内x^2
   end do

   vx0 = 1.0!设定初始时刻速度
   vy0 = 1.0
   do j = 1, t !求n次模拟的C(t)=<v(t)*v(0)>/2
      C(j) = 0
      if (j == 1) then
         do i = 1, n
            C(j) = C(j) + vx0*vx0/(2.0*real(n)) + vy0*vy0/(2.0*real(n))
         end do
      else
         do i = 1, n
            C(j) = C(j) + vx(i, j - 1)*vx0/(2.0*real(n)) + vy(i, j - 1)*vy0/(2.0*real(n))
         end do
      end if
   end do
   open (unit=4, file='C_A1w05.txt')
   do i = 1, t
      write (4, "(F8.4)") C(i)
   end do

   u = rand_num(seed())
   v = rand_num(int(array(1)*(2**31 - 1)))
   print *, (u)
   print *, (v)
   do i = 1, n!初始速度按照玻尔兹曼分布即高斯分布抽样
      vx_0(i) = sqrt(-2*log(u))*cos(2*pi*v)
      vy_0(i) = sqrt(-2*log(u))*sin(2*pi*v)
      u = rand_num(int(u*(2**31 - 1)))
      v = rand_num(int(v*(2**31 - 1)))
   end do
   do j = 1, t!粒子初始速度为Boltzmann分布下的速度自相关函数
      C(j) = 0
      if (j == 1) then
         do i = 1, n
            C(j) = C(j) + vx_0(i)*vx_0(i)/(2.0*real(n)) + vy_0(i)*vy_0(i)/(2.0*real(n))
         end do
      else
         do i = 1, n
            C(j) = C(j) + vx(i, j - 1)*vx_0(i)/(2.0*real(n)) + vy(i, j - 1)*vy_0(i)/(2.0*real(n))
         end do
      end if
   end do
   open (unit=9, file='C_A1w05_1.txt')
   do i = 1, t
      write (9, "(F8.4)") C(i)
   end do

   !输出第一次随机游走坐标的变化
   open (unit=5, file='walkx.txt')
   open (unit=6, file='walky.txt')
   do j = 1, t
      write (5, "(F8.4)") x(1, j, 1)
      write (6, "(F8.4)") x(1, j, 2)
   end do

   !输出t=100时粒子的分布情况
   open (unit=7, file='distribution_x.txt')
   open (unit=8, file='distribution_y.txt')
   do i = 1, n
      write (7, "(F8.4)") x(i, t, 1)
      write (8, "(F8.4)") x(i, t, 2)
   end do

   do i = 1, 9
      close (unit=i)
   end do

   print *, "finish"
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

   real*8 function rand_num(seed)!根据上一个使用的随机数生成下一个随机数
      implicit none
      integer :: n, seed, z
      integer i, m, a, r, q
      m = 2**31 - 1
      a = 7**5
      r = mod(m, a)
      q = m/a

      z = a*mod(seed, q) - r*(seed/q)
      if (z < 0) z = z + m
      rand_num = real(z)/real(m)
      return
   end function

end program Random_walk
