program half_sphere
   !use random,only:random_array !此处试用其他文件module,为提交方便直接在文件内定义
   implicit none
   real pi
   real, allocatable :: theta(:), phi(:), r(:), u(:), v(:),t(:)
   integer n, i, count
   pi=4*atan(1.0)!定义pi值
   n = 10**4
   t=2*random_array(2*n)-1
   allocate (u(n))
   allocate (v(n))
   do i=1,n
   u(i) = t(2*i-1)
   v(i) = t(2*i)
   end do

   r = (u**2 + v**2)**(0.5)
   theta = asin(2*r*(1 - r**2)**(0.5))!半球面，theta只取0到pi/2
   allocate (phi(n))
   count=0
   do i = 1, n
      valid: if (r(i) <= 1) then!检验取r<1的u,v值
         count=count+1

         if ((u(i) > 0) .and. (v(i) > 0)) then!根据u,v的正负判断求取phi角的象限
            phi(i) = acos(u(i)/r(i))
         end if
         if ((u(i) > 0) .and. (v(i) < 0)) then
            phi(i) = asin(v(i)/r(i))+2*pi
         end if
         if ((u(i) < 0) .and. (v(i) > 0)) then
            phi(i) = acos(u(i)/r(i))
         end if
         if ((u(i) < 0) .and. (v(i) < 0)) then
            phi(i) = -acos(u(i)/r(i))+2*pi
         end if
      end if valid
   end do

   open(unit=10,file='halfsphere.txt')!x写入文件
   do i=1,n
      valid3: if (r(i) <= 1) then
         write(10,"(F8.6XF8.6)")theta(i),phi(i)!此为输出8位(7小数位)，中间空格X
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
