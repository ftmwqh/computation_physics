program Metropolis_Hasting
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer n, i, j
   real*8 alpha, beta, gamma, rand, r, Integral, efficiency
   real*8, allocatable:: x0(:), x(:), rand_array(:)

   n = 10**6!粒子数
   m = 2**31 - 1!随机数列最大值
   alpha = 2.
   beta = 1.

   allocate (x0(n))
   allocate (x(n))
   allocate (rand_array(2*n))
   rand_array=random_array(2*n,seed())

   open (unit=1, file="error_1.txt")
   open (unit=2, file="efficiency_1.txt")
   gamma = 0.1
   do while (gamma < 50.0)
      rand = rand_num(seed())
      do i = 1, n!进行抽样得到x0
         x0(i) = -gamma*log(rand_array(i))
      end do

      do i = 1, n!对于所抽取的序列进行metropolis hasting方法
         if (i == 1) then!初始位置为x=1
            r = (x0(i)/1.0)**(alpha - 1.0)*exp(-(x0(i) - 1.0)/beta)*exp((x0(i) - 1.0)/gamma)
            if (rand_array(n+i) < 1 .and. rand_array(n+i) < r) x(i) = x0(i)!要注意这里的判据对象
            if (rand_array(n+i) > 1 .or. rand_array(n+i) > r) x(i) = 1.0
         else
            r = (x0(i)/x(i - 1))**(alpha - 1.0)*exp(-(x0(i) - x(i - 1))/beta)*exp((x0(i) - x(i - 1))/gamma)
            if (rand_array(n+i) < 1 .and. rand_array(n+i) < r) x(i) = x0(i)
            if (rand_array(n+i) > 1 .or. rand_array(n+i) > r) x(i) = x(i - 1)
         end if
      end do

      Integral = 0.0!进行积分
      do i = n/10 + 1, n
         Integral = Integral + (x(i) - alpha*beta)**2/(real(n)*9.0/10.0)
      end do
      write (1, *) abs(Integral - alpha*beta**2)!将绝对误差（误差绝对值）输出到txt文件

      efficiency = 0.0
      do i = 2, n
         if (x(i) /= x(i - 1)) efficiency = efficiency + 1.0/real(n)
      end do
      write (2, *) efficiency!将抽样效率输出到txt文件
      gamma = gamma + 0.1
   end do

   !取p(x)=(x-alpha*beta)^2f(x)时的模拟
   open (unit=3, file="error_p2.txt")
   open (unit=4, file="efficiency_p2.txt")
   gamma = 0.1
   do while (gamma < 100.0)
      do i = 1, n!进行抽样得到x0
         x0(i) = -gamma*log(rand_array(i))
      end do

      do i = 1, n!对于所抽取的序列进行metropolis hasting方法
         if (i == 1) then!初始位置为x=1
            r = ((x0(i) - alpha*beta)/(1.0 - alpha*beta))**2 &!这行太长了，会产生line-truncation报错，加'&'换行即可
                *(x0(i)/1.0)**(alpha - 1.0)*exp(-(x0(i) - 1.0)/beta)*exp((x0(i) - 1.0)/gamma)
            if (rand_array(n+i) < 1 .and. rand_array(n+i) < r) x(i) = x0(i)!要注意这里的判据对象
            if (rand_array(n+i) > 1 .or. rand_array(n+i) > r) x(i) = 1.0
         else
            r = ((x0(i) - alpha*beta)/(x(i - 1) - alpha*beta))**2 &
                *(x0(i)/x(i - 1))**(alpha - 1.0)*exp(-(x0(i) - x(i - 1))/beta)*exp((x0(i) - x(i - 1))/gamma)
            if (rand_array(n+i) < 1 .and. rand_array(n+i) < r) x(i) = x0(i)
            if (rand_array(n+i) > 1 .or. rand_array(n+i) > r) x(i) = x(i - 1)
         end if
      end do

      Integral = 0.0
      do i = n/10 + 1, n
         Integral = Integral + alpha*beta**2/(real(n)*9.0/10.0)
      end do
      write (3, *) abs(Integral - alpha*beta**2)!将绝对误差（误差绝对值）输出到txt文件

      efficiency = 0.0
      do i = 2, n
         if (x(i) /= x(i - 1)) efficiency = efficiency + 1.0/real(n)
      end do
      write (4, *) efficiency!将抽样效率输出到txt文件
      gamma = gamma + 0.1
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

   real*8 function T(x, gamma)!试探步函数T
      implicit none
      real*8 x, gamma
      T = 0.5*exp(-x/gamma)
      return
   end function

   real*8 function f(x, alpha, beta, Gamma_func)!Gamma分布函数
      implicit none
      real*8 x, alpha, beta, Gamma_func
      f = 1/(beta*Gamma_func)*(x/beta)**(alpha - 1)*exp(-x/beta)
      return
   end function

end program Metropolis_Hasting
