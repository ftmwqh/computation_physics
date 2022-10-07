PROGRAM Compare_16807_Fibonacci
   integer :: m, a, r, q, l, k, delay_1, delay_2
   integer :: z(100000000), seed
   real :: x(100000000), f(100000000), C_l
!16807产生器
   m = 2**31 - 1
   a = 7**5
   r = mod(m, a)
   q = m/a
   z(1) = seed()                          !初始值由seed产生
   print *, "seed:", (z(1))
   do i = 1, 99999999
      z(i + 1) = a*mod(z(i), q) - r*(z(i)/q)
      if (z(i + 1) < 0) z(i + 1) = z(i + 1) + m
   end do                                 !Schrage方法产生随机数序列z

   x = real(z)/real(m)                      !随机数序列归一化为x

!FIbonacci延迟器
   delay_1 = 43
   delay_2 = 29                        !取Fibonacci的延迟分别为delay_1,delay_2
   Fibonacci: do i = delay_1 + 1, 99999999 !利用前面16807产生的随机数前delay_1位，产生后面的随机数
      z(i) = mod(z(i - delay_1) - z(i - delay_2), m)
   end do Fibonacci
   f = real(z)/m
   do i = 1, 10**8
      if (f(i) < 0) then
         f(i) = f(i) + 1
      end if !随机数列中，由于Fibonacci相减产生，可能有负值，则对其+1，可得等效分布[0,1]正值
   end do

   print *, "16807 rate of Xn-1>Xn+1>Xn:", (rate(x, 10000000))   !打印出符合Xn-1>Xn+1>Xn的占比
   print *, "Fibonacci rate of Xn-1>Xn+1>Xn:", (rate(f, 10000000))!FibonacciXn-1>Xn+1>Xn的占比
   read (*, *)

END PROGRAM Compare_16807_Fibonacci

real function rate(x, N)                      !求Xn-1>Xn+1>Xn的占比
   integer N
   real x(N), count_x
   count_x = 0
   do i = 2, N - 1                             !符合逻辑关系则count_x+1,记录符合的Xn的数量
      if (x(i - 1) > x(i + 1) .and. x(i + 1) > x(i)) count_x = count_x + 1
   end do
   rate = real(count_x)/real(N)
   return
end function

integer function seed()                            !产生随机数种子值
   integer, dimension(8) :: time
   integer I_0
   call date_and_time(VALUES=time)             !获取系统时间
   I_0 = time(1) + 70*(time(2) + 12*(time(3) + 31*(time(5) + 23*(time(6) + 59*time(7)))))
   seed = I_0
   return
end function
