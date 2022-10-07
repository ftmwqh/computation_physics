PROGRAM Schrage
   integer :: m, a, r, q, l, k, delay_1, delay_2
   integer :: z(100000000)
   real*8 :: C,x(100000000)!产生10**8个随机数点,应使用 real*8以免溢出

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
   l = 2                                    !每隔多少取一组随机数作为坐标
!    open (unit=1, file='x_coordinate.txt')   !将产生的X坐标放入文件x_coordinate.txt

!    do i = 1, 99999999
!       if (mod(i, l + 2) == 1) then
!          write (1, *) x(i)
!       end if
!    end do

!    open (unit=2, file='y_coordinate.txt')   !将产生的Y坐标放入文件y_coordinate.txt
!    do i = 1, 99999999
!       if (mod(i, l + 2) == 1) then
!          write (2, *) x(i + 1)
!       end if
!    end do                                 !已生成文件，由于数据多速度过慢，注释掉



   !打印出C_l检测结果

   l = 1
   C=C_l(x,l,50000000)
   print *, "16807C_l examine(l=1,N=5*10^7):", (C_l(x, l, 50000000))
   l = 2
   print *, "16807C_l examine(l=2,N=5*10^7):", (C_l(x, l, 50000000))
   l = 3
   print *, "16807C_l examine(l=3,N=5*10^7):", (C_l(x, l, 50000000))
   !打印出<x_k>检测结果
   k = 1
   print *, "16807<x_1> examine(k=1,N=10^2):", (x_k(x, k, 100))
   print *, "16807<x_1> examine(k=1,N=10^3):", (x_k(x, k, 1000))
   print *, "16807<x_1> examine(k=1,N=10^4):", (x_k(x, k, 10000))
   print *, "16807<x_1> examine(k=1,N=10^5):", (x_k(x, k, 100000))
   print *, "16807<x_1> examine(k=1,N=10^6):", (x_k(x, k, 1000000))

   k = 2
   print *, "16807<x_1> examine(k=2,N=10^2):", (x_k(x, k, 100))
   print *, "16807<x_1> examine(k=2,N=10^3):", (x_k(x, k, 1000))
   print *, "16807<x_1> examine(k=2,N=10^4):", (x_k(x, k, 10000))
   print *, "16807<x_1> examine(k=2,N=10^5):", (x_k(x, k, 100000))
   print *, "16807<x_1> examine(k=2,N=10^6):", (x_k(x, k, 1000000))

   read (*, *)

contains

real*8 function C_l(x, l, N)                !C_l检测函数
   real*8 :: pa_1, pa_2, pa_3, x(N+l)
   integer :: i, l, N
   pa_1 = 0
   loop1: do i = 1, N                      !产生的pa_1即为<x_k*x_k+1>
      pa_1 = pa_1 + x(i)*x(i + l)/N

   end do loop1


   pa_2 = 0
   loop2: do i = 1, N                      !pa_2为<x_k>
      pa_2 = pa_2 + x(i)/N
   end do loop2

   pa_3 = 0
   loop3: do i = 1, N                       !pa_3为<x_k^2>
      pa_3 = pa_3 + x(i)**2/N
   end do loop3

   C_l = (pa_1 - pa_2**2)/(pa_3 - pa_2**2)    !计算C_l结果，并返回
   return
end function

real*8 function x_k(x, k, N)                     !<x^k>检测函数
   real*8 sum, x(N)
   integer N, k
   sum = 0
   do i = 1, N                                 !求<x^k>，并返回
      sum = sum + x(i)**k           !real可能会溢出，先进行除法
      
   end do
   x_k=sum/N
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
END PROGRAM Schrage!加contains就不需要再program中声明函数名称，在外部定义则需要声明