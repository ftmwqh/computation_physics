program Markov
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer n, i, j
   real*8 x0, y0, beta, gamma, rand, r, Integral, efficiency
   real*8, allocatable:: xy(:, :), ave(:),standard(:)

   n = 10**6!粒子数
   m = 2**31 - 10!随机数列最大值
   x0 = 5.0
   y0 = 5.0!设定初始坐标

   allocate (xy(2, n))
   allocate (ave(3))
   allocate(standard(3))
!beta=0.2
   standard(1)=1.561886003600077
   standard(2)=1.561886003600077
   standard(3)=2*1.561886003600077
   do i=1,6!结果绝对误差随模拟步数变化
      n=10**i
      beta=0.2
      xy = Metropolis(n, beta, x0, y0)
      ave = average(n, xy)
      print *,i,(abs(ave-standard))
   end do
!beta=1
   standard(1)=1.682470341271205
   standard(2)=1.682470341271205
   standard(3)=2*1.682470341271205
   do i=1,6!结果绝对误差随模拟步数变化
      n=10**i
      beta=1
      xy = Metropolis(n, beta, x0, y0)
      ave = average(n, xy)
      print *,i,(abs(ave-standard))
   end do
!beta=5
   standard(1)=1.947830003176043
   standard(2)=1.947830003176043
   standard(3)=2*1.947830003176043
   do i=1,6!结果绝对误差随模拟步数变化
      n=10**i
      beta=5
      xy = Metropolis(n, beta, x0, y0)
      ave = average(n, xy)
      print *,i,(abs(ave-standard))
   end do

   print *, "beta=0.2:"!beta=0.2下结果
   beta = 0.2
   xy = Metropolis(n, beta, x0, y0)
   open (unit=1, file="02_x.txt")
   open (unit=2, file="02_y.txt")
   do i = 1, n
      write (1, *) xy(1, i)
      write (2, *) xy(2, i)
   end do
   ave = average(n, xy)
   print *, (ave)
   close(1)
   close(2)   

   print *, "beta=1:"!beta=1下结果
   beta = 1.0
   xy = Metropolis(n, beta, x0, y0)
   open (unit=3, file="1_x.txt")
   open (unit=4, file="1_y.txt")
   do i = 1, n
      write (3, *) xy(1, i)
      write (4, *) xy(2, i)
   end do
   ave = average(n, xy)
   print *, (ave)
   close(4)
   close(3)   

   print *, "beta=5:"!beta=5下结果
   beta = 5.0
   xy = Metropolis(n, beta, x0, y0)
   ave = average(n, xy)
   print *, (ave)
   open (unit=5, file="5_x.txt")
   open (unit=6, file="5_y.txt")
   do i = 1, n
      write (5, *) xy(1, i)
      write (6, *) xy(2, i)
   end do
   close(6)
   close(5)   

   print*, "finish"!问题：打开文件后无法print到屏幕上，程序直接闪退，需要加断点
   read (*, *)
   stop
contains
   integer function seed()                            !产生随机数种子值
      integer, dimension(8) :: time
      integer I_0
      call date_and_time(VALUES=time)             !获取系统时间
      I_0 = time(1) + 70*(time(2) + 12*(time(3) + 31*(time(5) + 23*(time(6) + 59*time(7)))))
      seed = I_0
      return
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

 function random_array(n, seed) result(array)!产生长度为n的[0,1]的随机数
      implicit none
      integer :: n, seed
      integer, allocatable ::  z(:)
      integer i, m, a, r, q
      real*8, allocatable :: array(:)
      allocate (z(n))
      allocate (array(n))
      m = 2**31 - 1
      a = 7**5
      r = mod(m, a)
      q = m/a
      z(1) = seed

      do i = 1, n - 1                          !数组超上限会出现SIGSEGV或者自己断点的情况
         z(i + 1) = a*mod(z(i), q) - r*(z(i)/q)
         if (z(i + 1) < 0) z(i + 1) = z(i + 1) + m

      end do
      array = real(z)/real(m)

   end function

   real*8 function Hamilton(x, y)!哈密顿量值
      implicit none
      real*8 x, y
      Hamilton = -2.0*(x**2 + y**2) + (x**4 + y**4)/2 + (x - y)**4/2.0
      return
   end function

   real*8 function min_2(a, b)!求两个数之间的较小数
      real*8 a, b
      if (a > b) min_2 = b
      if (a <= b) min_2 = a
      return
   end function

   function Metropolis(n, beta, x0, y0) result(xy)!Metropolis方法生成x，y的Markov链
      implicit none
      integer i, j, k, m, n
      real*8 x0, y0, beta, rand, x, y, r
      real*8 xy(2, n),rands(3*n)
      m = 2**31 - 10!随机数列最大值
      xy(1, 1) = x0!初始坐标赋值，其中xy第一行为x坐标，第二行为y坐标
      xy(2, 1) = y0
      rand = rand_num(seed())!这里相连两个随机数关联性导致y^2总比x^2大，
      rands=random_array(3*n,seed())!选用生成随机数列，减小关联性
      do i = 2, n!进行抽样得到x0
         x = xy(1, i - 1) + 2*(rands(i) - 0.5)!行走一步在[-1,1]之间分布
         rand = rand_num(int(rand*m))
         y = xy(2, i - 1) + 2*(rands(i+n) - 0.5)
         rand = rand_num(int(rand*m))
         !判定是否进行位移
         if (rands(2*n+i) > min_2(num1, exp(-beta*(Hamilton(x, y) - Hamilton(xy(1, i - 1), xy(2, i - 1)))))) then
            xy(1, i) = xy(1, i - 1)
            xy(2, i) = xy(2, i - 1)
         else
            xy(1, i) = x
            xy(2, i) = y
         end if
         rand = rand_num(int(rand*m))
      end do
   end function

   function average(n, xy) result(ave)!求x^2,y^2,x^2+y^2的系综平均
      integer i, j, k, n
      real*8 xy(2, n), ave(3)
      ave(1) = 0
      ave(2) = 0
      do i = n/10 + 1, n
         ave(1) = ave(1) + xy(1, i)**2/(real(n)*9.0/10.0)
         ave(2) = ave(2) + xy(2, i)**2/(real(n)*9.0/10.0)
      end do
      ave(3) = ave(1) + ave(2)
   end function

end program Markov
