program Chaos
   real*8, parameter :: num4 = 4., num1 = 1.
   real*8, parameter :: pi = num4*atan(num1)!定义pi
   integer n, i, j, count, step, start, judge, k0
   real*8 x0, y0, lambda, delta, epsilon, min
   real*8, allocatable:: x(:), diverge(:), F(:), d(:), alpha(:)
   step = 5000
   n = 100
   allocate (d(6))
   allocate (diverge(6))
   allocate (alpha(5))
   allocate (F(4))
   allocate (x(n + step))!前step步迭代认为达到稳定，输出后n个值作为迭代稳定结果

    open (unit=1, file="result1.txt")
    lambda=-10.0
   delta=0.01
   do while(lambda<=10)!首先输出在lambda[-10,10]范围内
    x0=1.0
    x(1)=iteration(lambda,x0)
    do i=2,step+n
        x(i)=iteration(lambda,x(i-1))
        if(i>step)then
            write(1,*)x(i)
        end if
    end do
    write(1,*)","
    lambda=lambda+delta
   end do

        open (unit=2, file="result2.txt")
    lambda=0.0
   delta=0.0001
   do while(lambda<=2.0)!输出在lambda[0,2]范围内
    x0=1.0
    x(1)=iteration(lambda,x0)
    do i=2,step+n
        x(i)=iteration(lambda,x(i-1))
        if(i>step)then
            write(2,*)x(i)
        end if
    end do
    write(2,*)","
    lambda=lambda+delta
   end do

    open (unit=3, file="result3.txt")
    lambda=-2.0
   delta=0.0001
   do while(lambda<=0.0)!输出在lambda[-2,0]范围内
    x0=1.0
    x(1)=iteration(lambda,x0)
    do i=2,step+n
        x(i)=iteration(lambda,x(i-1))
        if(i>step)then
            write(3,*)x(i)
        end if
    end do
    write(3,*)","
    lambda=lambda+delta
   end do

!计算Feigenbaum常数
   lambda = 0.6
   delta = 0.0000001
   epsilon = 0.000001
   k = 1
   k0 = 1
   y0 = 0.5
   do while (lambda <= 0.9)!周期分叉在lambda[0.6,0.9]范围内
      x0 = 1.0
      x(1) = iteration(lambda, x0)
      count = 0
      do i = 2, step + n
         x(i) = iteration(lambda, x(i - 1))
         if (i > step) then
            if (count == 0) then
               count = count + 1
               start = i
            else
               judge = 1
               do j = start, i - 1
                  judge = judge*unique(x(i), x(j), epsilon)!通过judge0或1判断新的值是否与之前已经计入的一样还是新的不一样的值
                  if (judge == 0) exit
               end do
               if (judge == 1) count = count + 1
            end if

         end if
      end do
      if (count == 2**k) then
         diverge(k) = lambda
         k = k + 1
      end if

      do i = step + 1, step + 2**k0
         if (unique(x(i), y0, epsilon) == 0) then!当该点x值等于0.5时，记录该点
            if (count == 2**k0) then
               min = 1.0
               do j = i + 1, i + 2**k0 - 1
                  if (abs(x(i) - x(j)) < min) min = abs(x(i) - x(j))
               end do!这里取距离x=0.5点在同意lambda上最接近的周期点，计算另一分型与0.5的距离
               d(k0) = min
               k0 = k0 + 1
            end if
         end if
      end do

      lambda = lambda + delta
   end do
   !计算feigrnbaum常数
   do i = 1, 4
      F(i) = (diverge(i + 1) - diverge(i))/(diverge(i + 2) - diverge(i + 1))
   end do
   do i = 1, 5
      alpha(i) = d(i)/d(i + 1)
   end do

!打表输出计算的Feigenbaum常数
   write (*, "(A10,A15,A25)") "diverge", "lambda", "Feigenbaum Constant"
   write (*, "(A10,1F15.7)") "1->2", diverge(1)
   write (*, "(A10,1F15.7,1F25.7)") "2->4", diverge(2), F(1)
   write (*, "(A10,1F15.7,1F25.7)") "4->8", diverge(3), F(2)
   write (*, "(A10,1F15.7,1F25.7)") "8->16", diverge(4), F(3)
   write (*, "(A10,1F15.7,1F25.7)") "16->32", diverge(5), F(4)
   write (*, "(A10,1F15.7)") "32->64", diverge(6)

   write (*, "(A10,A15,A25)") "diverge", "d", "Feigenbaum Constant"
   write (*, "(A10,1F15.7,1F25.7)") "2", d(1), alpha(1)
   write (*, "(A10,1F15.7,1F25.7)") "4", d(2), alpha(2)
   write (*, "(A10,1F15.7,1F25.7)") "8", d(3), alpha(3)
   write (*, "(A10,1F15.7,1F25.7)") "16", d(4), alpha(4)
   write (*, "(A10,1F15.7,1F25.7)") "32", d(5), alpha(5)
   write (*, "(A10,1F15.7)") "64", d(6)

   print *, "finish"
   read (*, *)
   stop
contains

   real*8 function iteration(lambda, x)!迭代的方程
      real*8 lambda, x
      iteration = lambda*sin(pi*x)
      return
   end function

   integer function unique(x, standard, epsilon)!判断两个相近的数是否为同一个数
      real*8 x, standard, epsilon
      if (abs(x - standard) > epsilon) then!相差大于判断值时，认为不是同一个数
         unique = 1
      else
         unique = 0
      end if
      return
   end function
end program Chaos
