program test_allocatable_functionReturn
      implicit none

      integer, allocatable :: fun(:),sub(:),noa(:)
      integer, parameter :: low = -4
      integer, parameter :: hi  = 3

      call testsub(sub,low,hi)!没有result的函数要用call
      fun = testfun(low,hi)
      noa = testfun_noalloc(low,hi)

      print '(4(a,i3),a)', 'testsub:  lbound=',lbound(sub),'(expected = ',low,'), ubound=',ubound(sub),'(expected = ',hi,')'
       print '(4(a,i3),a)', 'testfun:  lbound=',lbound(fun),'(expected = ',low,'), ubound=',ubound(fun),'(expected = ',hi,')'
       print '(4(a,i3),a)', 'no alloc: lbound=',lbound(noa),'(expected = ',low,'), ubound=',ubound(noa),'(expected = ',hi,')'
    print*,(noa)
        read(*,*)
      contains
             pure function testfun_noalloc(low,hi) result(array)
                integer, intent(in) :: low,hi
                integer :: array(low:hi)
                integer :: i
                 forall(i=low:hi) array(i) = i
             end function testfun_noalloc


             pure function testfun(low,hi) result(array)
                integer, allocatable :: array(:)
                integer, intent(in) :: low,hi
                integer :: i
                 allocate(array(low:hi))
                 forall(i=low:hi) array(i) = i
             end function testfun


             pure subroutine testsub(array,low,hi)
                integer, intent(out), allocatable :: array(:)
                integer, intent(in) :: low,hi
                integer :: i
                 allocate(array(low:hi))
                 forall(i=low:hi) array(i) = i

             end subroutine testsub

end program