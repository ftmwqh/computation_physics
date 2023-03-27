program Central_Limit
   real*8, allocatable ::x(:), y(:), z(:), u(:), p(:), array(:)
   integer k, N, i, j, count
   real*8 a, b, c, d, e, M, sum, lambda, mu
   i = 1

   k = 10**5

   allocate (x(k))
   allocate (y(k))
   allocate (z(k))
   allocate (u(k))
   allocate (p(k*N))
   allocate (array(k*N))
   N = 2
   !计算泊松分布
   lambda = 5.
   p = poisson(k, lambda)
   print *, "Poisson average  (n=", (k), ", lambda=", (lambda), ")", (sum(p, k)/k)
   y = Poisson_central(k, N, lambda, "poisson_2.txt")

   !计算指数分布
   mu = 5.
   u = random_array(k)
   x = -mu*log(u)
   print *, "Exp average  (n=", (k), ", mu=", (mu), ")", (sum(x, k)/k)
   y = Exp_central(k, N, mu, "exp_2.txt")

   a=0.
   b=1.
   y=Uniform(k,N,a,b,"uniform_2.txt")
   y=Gauss(k,N,"gauss_2.txt")

   N = 5!N=5的情况
   y = Poisson_central(k, N, lambda, "poisson_5.txt")
   y = Exp_central(k, N, mu, "exp_5.txt")
   y=Uniform(k,N,a,b,"uniform_5.txt")
   y=Gauss(k,N,"gauss_5.txt")

   N = 10!N=10的情况
   y = Poisson_central(k, N, lambda, "poisson_10.txt")
   y = Exp_central(k, N, mu, "exp_10.txt")
   y=Uniform(k,N,a,b,"uniform_10.txt")
   y=Gauss(k,N,"gauss_10.txt")
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

   function poisson(n, lambda) result(poisson_array)!产生泊松分布的随机数列
      integer n, i, j, count
      real*8 lambda, u(int(5*lambda)*n), poisson_array(n), p
      u = random_array(int(5*lambda)*n)!取一个较大的均匀分布随机数组，防止产生泊松分布时超过数列上限
      j = 1

      do i = 1, n
         p = 1.0
         count = 0!记录乘以几次均匀分布的随机变量
         do while (p >= exp(-lambda))
            p = p*u(j)
            j = j + 1
            count = count + 1
         end do
         poisson_array(i) = count - 1
      end do
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

   real*8 function square_sum_array(x, n)!对一个数列求平方和
      real*8::x(n)
      integer n
      square_sum_array = 0
      do i = 1, n
         square_sum_array = square_sum_array + x(i)**2
      end do
      return
   end function

   function Poisson_central(k, N, lambda, filename) result(t)!参数lambda的泊松分布
      integer k, N, i, j
      real*8 p(k*N), x(N), lambda, t(k)
      character(len=*) filename
      p = poisson(k*N, lambda)!p为生成所有误差变量的总泊松分布随机数列
      open (unit=1, file=filename)
      do i = 1, k
         x = p((i - 1)*N + 1:i*N)
         t(i) = (sum_array(x, N)/N - lambda)/(sqrt((square_sum_array(x, N)/N - (sum_array(x, N)/N)**2)/N))
         if (sqrt((square_sum_array(x, N)/N - (sum_array(x, N)/N)**2)/N) /= 0) then!防止分母为0，输出无穷大
            write (1, *) t(i)!t即为中心极限定理的误差随机变量
         end if
      end do
   end function

   function Exp_central(k, N, mu, filename) result(t)!参数为mu的指数分布
      integer k, N, i, j
      real*8 u(k*N), x(N), mu, t(k)
      character(len=*) filename
      u = random_array(k*N)!u为生成所有误差变量的总均匀分布随机数列
      open (unit=2, file=filename)

      do i = 1, k
         x = -mu*log(u((i - 1)*N + 1:i*N))
         t(i) = (sum_array(x, N)/N - mu)/(sqrt((square_sum_array(x, N)/N - (sum_array(x, N)/N)**2)/N))
         write (2, *) t(i)!t即为中心极限定理的误差随机变量
      end do
   end function

   function Uniform(k, N,a,b, filename) result(t)![a,b]均匀分布的随机数
      integer k, N, i, j
      real*8 u(k*N), x(N), mu, t(k),a,b
      character(len=*) filename
      u = random_array(k*N)!u为生成所有误差变量的总均匀分布随机数列
      mu=(a+b)/2
      open (unit=3, file=filename)

      do i = 1, k
         x = u((i - 1)*N + 1:i*N)
         t(i) = (sum_array(x, N)/N - mu)/(sqrt((square_sum_array(x, N)/N - (sum_array(x, N)/N)**2)/N))
         write (3, *) t(i)!t即为中心极限定理的误差随机变量
      end do
   end function

   function Gauss(k, N, filename) result(t)!高斯分布
      integer k, N, i, j
      real*8 u(k*N),v(k*N),array(2*k*N) , x(N), mu, t(k),pi
      character(len=*) filename
      array=random_array(2*k*N)
      u = array(1:k*N)!u,v为生成所有误差变量的总均匀分布随机数列
      v = array(k*N+1:2*k*N)
      pi=4.*atan(1.)
      open (unit=4, file=filename)

      do i = 1, k
         x = sqrt(-2*log(u((i - 1)*N + 1:i*N)))*cos(2*pi*v((i - 1)*N + 1:i*N))
         t(i) = (sum_array(x, N)/N )/(sqrt((square_sum_array(x, N)/N - (sum_array(x, N)/N)**2)/N))
         write (4, *) t(i)!t即为中心极限定理的误差随机变量
      end do
   end function
end program Central_Limit
