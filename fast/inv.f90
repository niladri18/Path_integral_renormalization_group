module inv
 use variable
 implicit none

 private
 public :: minv

 contains 
 
 subroutine minv(a,n,det,ainv)
! invert a matrix using LAPACK
  real(kind=8), dimension(:,:), intent(inout) :: a
  real(kind=8), dimension(:,:), intent(out) :: ainv
  real(kind=8), intent(out) :: det
  integer, intent(in) :: n
  integer, allocatable :: indx(:)

  allocate(indx(n))


  call dgetrf(n,n,a,n,indx,info)
  if (info/=0) then
     write(*,*) 'Info1=',info,' in minv'
     stop
  endif
  det=1.0d0
  do i=1,n
     if (indx(i)/=i) then
        det=-1.0d0*det*a(i,i)
     else
        det=det*a(i,i)
     endif
  enddo

! print*,"entered minv",lwork

  call dgetri(n,a,n,indx,WORK,LWORK,info)
  if (info/=0) then
     write(*,*) 'Info2=',info,' in minv'
     stop
  endif
  ainv=a

  deallocate(indx)

 end subroutine minv

end module inv
