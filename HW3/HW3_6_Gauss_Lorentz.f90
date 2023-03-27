program Gauss_Lorentz
!由计算分析得取a=3,b=2,c=1.5可以使c/(1+bx^4)>=exp(-ax^2)在实轴上全部成立
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi为parameter形式
   real*8, allocatable ::xi_x(:), xi_y(:), xi_1(:), xi_2(:),x(:),array(:)
   integer n, i, j
   real*8 epsilon, a, b
   j=1
   n = 10**6
   epsilon = 10.**(-3)
   allocate(xi_1(n))
   allocate(xi_2(n))
   allocate(xi_x(n))
   allocate(xi_y(n))
   allocate(x(n))
   allocate(array(2*n))
   array = random_array(2*n)!array前n用于存储xi_1的随机数，后n存储xi_2的随机数
   xi_1=array(1:n)
   xi_2=array(n+1:2*n)
   do i = 1, n
      if (xi_1(i) > 0.5) then!由函数在实数域上单调递增，判断零点所在的区间是在正半轴还是负半轴
         a = 0.0
         b = 5.0
         do while (equation(xi_1(i), b) < 0)
            a = a + 5.0
            b = b + 5.0
         end do
         do while (max(abs(b - a),abs(equation(xi_1(i), b)-equation(xi_1(i), a))) > epsilon)
            if (equation(xi_1(i), (a + b)/2) > 0) then
               b = (a + b)/2
            else
               a = (a + b)/2
            end if
         end do
         xi_x(i) = (a + b)/2
      else
         a=-5.0
         b=0.0
         do while (equation(xi_1(i), a) > 0)
            a = a - 5.0
            b = b - 5.0
         end do
         do while (max(abs(b - a),abs(equation(xi_1(i), b)-equation(xi_1(i), a))) > epsilon)
            if (equation(xi_1(i), (a + b)/2) > 0) then
               b = (a + b)/2
            else
               a = (a + b)/2
            end if
         end do
         xi_x(i) = (a + b)/2
         
      end if
      xi_y(i)=xi_2(i)*F(xi_x(i))
         if(xi_y(i)<p(xi_x(i))) then!判断满足条件时将xi_x列入x中，最后x的有效个数为j-1
            x(j)=xi_x(i)
            j=j+1
         end if

   end do

   open(unit=10,file='x.txt')!x写入文件
   do i=1,j-1
      write(10,*)x(i)
   end do 
   print*,"rate:",(real(j)/n) !计算抽样率
   print*,"end"
      


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

   real*8 function F(x)
      real*8 x,b,c
      b=2.
      c=1.5
      F = c/(1.+b*x**4)
      return
   end function

   real*8 function p(x)
      real*8 x,a
      a=3.
      p = sqrt(a/pi)*exp(-a* x**2)!给p乘以归一化因子
      return
   end function

   real*8 function equation(x1, xx)
      real*8 x1, xx
      equation = (log((sqrt(2.)*xx**2 + 2.**0.75*xx + 1.)/(sqrt(2.)*xx**2 - 2.**0.75*xx + 1)) &
                  + 2*atan(2.**0.75*xx + 1) - 2*atan(1 - 2.**0.75*xx))/(4*pi) + 0.5 - x1
      return
   end function
end program Gauss_Lorentz
