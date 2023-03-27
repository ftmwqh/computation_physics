program Fractal
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer, allocatable ::mesh(:, :)
   integer n, i, j, k, L, center, x, y, m, r(9), count(9)
   real*8 edge, start, rand

   n = 10**5!粒子数
   m = 2**31 - 1!随机数列最大值
   L = 1001!格点阵的边长
   center = (L + 1)/2!格点阵中心坐标
   allocate (mesh(L, L))

   mesh = DLA(n, L)
   count = 0
   open (unit=7, file="r.txt")
   open (unit=8, file="count.txt")
   do i = 1, 9!sandbox法
      r(i) = 2**i
      do j = center - r(i)/2 + 1, center + r(i)/2
         do k = center - r(i)/2 + 1, center + r(i)/2
            if (mesh(j, k) == 1) count(i) = count(i) + 1
         end do
      end do
      write (7, *) log(real(r(i)))
      write (8, *) log(real(count(i)))
   end do

   close (5)
   close (6)
   close (7)
   close (8)
   close (4)
   close (3)
   close (2)
   close (1)
   print *, (log(real(r)))
   print *, (log(real(count)))
   print *, (log(real(count))/log(real(r)))

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
      integer L, x, y, n_eff
      real*8 edge

      !达到边界返回2
      if (sqrt(real((x - (L + 1)/2)**2 + (y - (L + 1)/2)**2)) > edge) then
         judge = 2
         !周围有已经生长的(值为1)则返回1
      else if (mesh(x + 1, y) == 1 .or. mesh(x - 1, y) == 1 .or. mesh(x, y + 1) == 1 .or. mesh(x, y - 1) == 1) then
         judge = 1
      else
         judge = 0
      end if
      return
   end function

   function DLA(n, L) result(mesh)
      integer, allocatable ::mesh(:, :)
      integer n, i, j, L, center, x, y, m, count, output
      real*8 edge, start, rand, Rg, n_eff, rands(10**7), r_max

      m = 2**31 - 1!随机数列最大值
      center = (L + 1)/2!格点阵中心坐标
      allocate (mesh(L, L))

      mesh = 0
      mesh(center, center) = 1!初始化网格，在中心设1为结晶核
      rands = random_array(10**7, seed())
      n_eff = 0
      count = 1
      output = 1
      r_max = 0
      open (unit=5, file="Rg.txt")
      open (unit=6, file="N_eff.txt")
      !对n个粒子进行随机游走
      particles: do i = 1, n
         if (2*(r_max + 5.) < (L - 1)/2 - 1) then
            edge = 2*(5.+r_max)!粒子消失的边缘
            start = r_max + 5.!作为初始粒子生成的圆半径
         else 
            edge = (L - 1)/2 - 1
            start = r_max + 5.!作为初始粒子生成的圆半径
         end if
         x = center + int(start*cos(2*pi*rands(count)))
         y = center + int(start*sin(2*pi*rands(count)))
         count = count + 1
         if (count == 10**7) then
            rands = random_array(10**7, seed())
            count = 1
         end if
         do while (.TRUE.)!进行随机游走，直到出边界或碰撞已经生长点
            x = walkx(x, rands(count))
            y = walky(y, rands(count))
            count = count + 1
            if (count == 10**7) then
               rands = random_array(10**7, seed())
               count = 1
            end if
            if (judge(mesh, L, x, y, edge) == 1) then
               n_eff = n_eff + 1.0!计算实际生长的粒子数
               mesh(x, y) = 1
               Rg = Rg + (x - center)**2 + (y - center)**2!加入新生长粒子的半径平方
               if (n_eff == 2**output) then
                  write (5, *) log(sqrt(Rg/n_eff))
                  write (6, *) log(n_eff)
                  output = output + 1
               end if
               if ((x - center)**2 + (y - center)**2 > r_max**2) r_max = sqrt(real((x - center)**2 + (y - center)**2))!更新最大半径
               exit
            else if (judge(mesh, L, x, y, edge) == 2) then
               exit
            end if
         end do
         ! print*,(edge),' ',(start)
         ! print*, (n_eff)
         !输出从100到10^5个粒子生长的过程
           if(i==100)then
              open(unit=1,file="mesh_2.txt")
              write(1,*)mesh
           else if(i==1000)then
              open(unit=2,file="mesh_3.txt")
              write(2,*)mesh
           else if(i==10**4)then
              open(unit=3,file="mesh_4.txt")
              write(3,*)mesh
           else if(i==10**5)then
              open(unit=4,file="mesh_5.txt")
              write(4,*)mesh
           end if
      end do particles
   end function
end program Fractal
