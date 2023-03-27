program DLA
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer, allocatable ::mesh(:, :)
   real*8,allocatable::rands(:)
   integer n, i, j, L, center, x, y,m ,n_eff,count
   real*8 edge, start,rand,r_max

   n = 10**5!粒子数
   m = 2**31 - 1!随机数列最大值
   L = 1001!格点阵的边长
   center = (L + 1)/2!格点阵中心坐标
   allocate (mesh(L, L))
   allocate(rands(10**7))
   rands=random_array(10**7,seed())
   count=1

   mesh = 0
   r_max=0!当前最大半径
   mesh(center, center) = 1!初始化网格，在中心设1为结晶核
   rand=rand_num(seed())
   n_eff=0
   !对n个粒子进行随机游走
   particles: do i = 1, n
      ! if (2*real(i)**0.67 < (L - 1)/2-1) then
      !    edge = 2*real(i)**0.67!粒子消失的边缘
      !    start = 0.8*real(i)**0.67!作为初始粒子生成的圆半径
      ! else
      !    edge = (L - 1)/2-1
      !    start = 0.8*real(i)**0.5!作为初始粒子生成的圆半径
      ! end if
      ! if(i<100)then
      !    start=10
      !    edge=2*start
      ! else if(100<=i.and.i<300)then
      !    start=20
      !    edge=2*start
      ! else if(300<=i.and.i<1000)then
      !    start=60
      !    edge=2*start
      ! else if(1000<=i.and.i<10**4)then
      !    start=100
      !    edge=2*start
      ! else if(10**4<=i.and.i<5*10**4)then
      !    start=300
      !    edge=499
      ! else
      !    start=400
      !    edge=499
      ! end if
      if (2*(r_max+5.) < (L - 1)/2-1) then
         edge = 2*(5.+r_max)!粒子消失的边缘
         start = r_max+5.!作为初始粒子生成的圆半径
      else
         edge = (L - 1)/2-1
         start = r_max+5.!作为初始粒子生成的圆半径
      end if
      !print*,(r_max)
      x = center + int(start*cos(2*pi*rands(count)))
      y = center + int(start*sin(2*pi*rands(count)))
      !print*,(x),",",(y)
      !print*,"start==================================================",(rand)
      !!!!
      !!!!发现奇怪的事情，用单个生成的随机数rand_num会导致最终陷于循环，每次重生点位都是同一个，
      !!!!可能是由于丧失了一定的精度，导致周期变短，慎之慎之
      !!!!
      rand=rand_num(int(rand*m))
      count = count + 1
      if(count==10**7)then
         rands=random_array(10**7,seed())
         count=1
      end if
      do while (.TRUE.)!进行随机游走，直到出边界或碰撞已经生长点
         x = walkx(x, rands(count))
         y = walky(y, rands(count))
         rand=rand_num(int(rand*m))
         count = count + 1
         if(count==10**7)then
            rands=random_array(10**7,seed())
            count=1
         end if
         if (judge(mesh, L, x, y, edge) == 1) then
            mesh(x, y) = 1
            n_eff=n_eff+1!记录生长粒子数
            !print*,(n_eff),sqrt(real((x-center)**2+(y-center)**2))
            !print*,(x),",",(y)
            !print*,"r_now:",(sqrt(real((x-center)**2+(y-center)**2)))
            if((x-center)**2+(y-center)**2>r_max**2)r_max=sqrt(real((x-center)**2+(y-center)**2))!更新最大半径
            exit
         else if (judge(mesh, L, x, y, edge) == 2) then
            !print*,(x),",",(y)
            exit
         end if
      end do
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
      else if (rand >= 0.25 .and. rand <= 0.5) then
         x = x - 1
      end if
      walkx = x
      return
   end function

   integer function walky(y, rand)!y方向随机游走，根据随机数值判断游走方向，返回之后的坐标
      integer y
      real*8 rand
      if (rand > 0.5 .and. rand <= 0.75) then
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
      else
         judge = 0
      end if
      return
   end function
end program DLA
