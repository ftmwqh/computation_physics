program DBM
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer, allocatable ::mesh(:, :)
   integer n, i, j, k, t, L, center, x, y, m, next(2)
   real*8 edge, start, rand, eta
   real*8, allocatable :: array(:)
   t = 10**3!生长次数
   n = 10**2!求解拉普拉斯方程的montecarlo模拟粒子数
   m = 2**31 - 1!随机数列最大值
   L = 501!格点阵的边长
   eta = 1.0!生长速率参数
   center = (L + 1)/2!格点阵中心坐标
   allocate (mesh(L, L))
   allocate (array(t))
   array = random_array(t, seed())

   mesh = 0
   mesh(center, center) = 1!初始化网格，在中心设1为结晶核
   rand = rand_num(seed())

   do k = 1, t!进行t次生长!这里循环用i会出现跳变的问题，可能是grow函数里也用到了i，需要注意变量与函数中变量混用可能的问题
      next = grow(mesh, L, k, n, eta, array(k))
      mesh(next(1), next(2)) = 1
      if (k == 100) then!将不同生长次数的网格导出到txt
         open (unit=10, file="DBM_1_100.txt")
         write (10, *) mesh
         close (unit=10)
      else if (k == 300) then
         open (unit=11, file="DBM_1_300.txt")
         write (11, *) mesh
         close (unit=11)
      else if (k == 700) then
         open (unit=12, file="DBM_1_700.txt")
         write (12, *) mesh
         close (unit=12)
      else if (k == 1000) then
         open (unit=13, file="DBM_1_1000.txt")
         write (13, *) mesh
         close (unit=13)
      else if (k == 1500) then
         open (unit=14, file="DBM_1500.txt")
         write (14, *) mesh
         close (unit=14)
      else if (k == 2000) then
         open (unit=15, file="DBM_2000.txt")
         write (15, *) mesh
         close (unit=14)
      end if
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

   real*8 function sum_array(x, n)!对一个数列求和
      real*8::x(n)
      integer n
      sum_array = 0
      do i = 1, n
         sum_array = sum_array + x(i)
      end do
      return
   end function

   integer function walkx(x, rand)!x方向随机游走，根据随机数值判断游走方向，返回之后的坐标
      integer x
      real*8 rand
      if (rand < 0.25) then
         x = x + 1
      else if (rand > 0.25 .and. rand < 0.5) then
         x = x - 1
      end if
      walkx = x
      return
   end function

   integer function walky(y, rand)!y方向随机游走，根据随机数值判断游走方向，返回之后的坐标
      integer y
      real*8 rand
      if (rand > 0.5 .and. rand < 0.75) then
         y = y + 1
      else if (rand > 0.75) then
         y = y - 1
      end if
      walky = y
      return
   end function

   integer function judge(mesh, L, x, y, edge)!判断当前游走位置的状况
      integer mesh(L, L)
      integer L, x, y
      real*8 edge

      !达到边界返回2
      if (sqrt(real((x - (L + 1)/2)**2 + (y - (L + 1)/2)**2)) > edge) then
         judge = 2
         !周围有已经生长的(值为1)则返回1
      else if (mesh(x + 1, y) == 1 .or. mesh(x - 1, y) == 1 .or. mesh(x, y + 1) == 1 .or. mesh(x, y - 1) == 1) then
         judge = 1
      else!继续游走
         judge = 0
      end if
      return
   end function

   function grow(mesh, L, t, n, eta, rand_grow) result(coordinate)!求解生长位置
      integer mesh(L, L)
      integer L, i, j, k, n, m, t, x, y, start(4*t, 2), coordinate(2)
      real*8 rand_grow, edge, sum, rand, phi(n), eta, a, b
      real*8, allocatable:: p(:)
      m = 2**31 - 1
      rand = rand_num(seed())
      if (3*sqrt(real(t)) < (L - 1)/2 - 1) then!取遍历寻找闪电界限的变界
         edge = 3*sqrt(real(t))
      else
         edge = (L - 1)/2 - 1
      end if

      !遍历网格寻找闪电的边缘处，即01交界的点
      i = 0
      do x = (L + 1)/2 - int(edge), (L + 1)/2 + int(edge)
         do y = (L + 1)/2 - int(edge), (L + 1)/2 + int(edge)
            if (mesh(x, y) == 0 .and. judge(mesh, L, x, y, edge) == 1) then
               i = i + 1
               start(i, 1) = x
               start(i, 2) = y!交界总点数为i
            end if
         end do
      end do

      allocate (p(i))!每个边界点的扩散概率
      do j = 1, i
         do k = 1, n
            x = start(j, 1)!对于n次模拟求平均，都从该边界初始点开始
            y = start(j, 2)
            do while (.TRUE.)
               x = walkx(x, rand)
               y = walky(y, rand)

               rand = rand_num(int(rand*m))
               if (judge(mesh, L, x, y, edge) == 2) then!先判断是否离开边界，否则可能指标溢出
                  phi(k) = 0
                  exit
               else if (mesh(x, y) == 1) then
                  phi(k) = 1.0
                  exit

               end if
            end do
         end do
         sum = sum_array(phi, n)/real(n)
         p(j) = real(t + 1)*(1.0 - sum)**eta
      end do
      sum = sum_array(p, i)
      p = p/sum!归一化得到i个边界点的概率
      sum = 0
      do j = 1, i!对概率进行求和，随机数rand在该小区间范围内则返回该概率区域对应的边界点作为生长点
         sum = sum + p(j)
         if (rand_grow < sum) then
            coordinate(1) = start(j, 1)
            coordinate(2) = start(j, 2)
            exit
         end if
      end do
   end function

end program DBM
