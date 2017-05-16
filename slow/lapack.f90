module lapack

interface

function ddot(n,x,incx,y,incy)
  real(kind=8) :: ddot
  integer,intent(in) :: n,incx,incy
  real(kind=8),intent(in)::x(*),y(*)
end function ddot

subroutine daxpy(n,alpha,x,incx,y,incy)
  integer,intent(in) :: n,incx,incy
  real(kind=8),intent(in)::alpha,x(*)
  real(kind=8),intent(inout)::y(*)
end subroutine daxpy

subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  character,intent(in) :: jobz,uplo
  integer,intent(in):: n,lda,lwork
  integer,intent(out):: info
  real(kind=8),intent(inout) :: a(*)
  real(kind=8),intent(out) :: w(*)
  real(kind=8),intent(out) :: work(*)
end subroutine dsyev

      
subroutine dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  integer, intent(in) :: itype,n,lda,ldb,lwork
  integer, intent(out) :: info
  character, intent(in) :: jobz,uplo
  real(kind=8), intent(out) :: w(*),work(*)
  real(kind=8), intent(inout) :: a(*),b(*)
end subroutine dsygv

subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  character,intent(in) :: transa,transb
  integer,intent(in) :: m,n,k,lda,ldb,ldc
  real(kind=8), intent(in) :: a(*),b(*)
  real(kind=8), intent(out) :: c(*)
  real(kind=8), intent(in) :: alpha,beta
end subroutine dgemm

subroutine dgetrf(m,n,a,lda,ipiv,info)
  integer, intent(in) :: m,n,lda
  integer, intent(out) :: info,ipiv(*)
  real(kind=8), intent(inout) :: a(*)
end subroutine dgetrf

subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
  integer, intent(in) :: lda,lwork,n
  integer, intent(out) :: info
  integer, intent(in) :: ipiv(*)
  real(kind=8), intent(inout) :: a(*)
  real(kind=8), intent(out) :: work(*)
end subroutine dgetri

subroutine dscal(n,da,dx,incx)
  integer, intent(in) :: n,incx
  real(kind=8), intent(in) :: da
  real(kind=8), intent(inout) :: dx(*)
end subroutine

end interface

end module lapack
