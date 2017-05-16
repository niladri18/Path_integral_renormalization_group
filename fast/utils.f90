module util

 use variable
 use inv
 use lapack
 implicit none
 private
 public::greens,hamiltonian,modgs

contains

subroutine greens(basis1,basis2,green,green2,N,Ne,ndet)

  double precision,dimension(:,:),intent(in)::basis1,basis2
  double precision,dimension(:,:),intent(out)::green,green2
  integer,intent(in)::N,ne
  double precision,dimension(:,:),allocatable::prod,overlap,green1,overlapinv
  integer::i,j
  double precision,intent(out)::ndet
  double precision::det
  !print*,basis2
  allocate(prod(ne,ne))
  allocate(overlap(ne,ne))
  allocate(overlapinv(ne,ne))
  allocate(green1(N,ne))
  call dgemm('T','N',ne,ne,N,1.0d0,basis1,N,basis2,N,0.0d0,prod,ne)
  overlap=prod
  !print*,prod
  call minv(overlap,ne,det,overlapinv)
  !print*,"det",det
  ndet=det
  !print*,"inv",overlapinv
  call dgemm('N','N',N,ne,ne,1.0d0,basis2,N,overlapinv,ne,0.d0,green1,N)
  call dgemm('N','T',N,N,ne,1.0d0,green1,N,basis1,N,0.0d0,green,N)
  
  do i=1,n
     do j=1,n
        if(i==j)then
           green2(i,j)=1.0d0-green(i,j)
        else
           green2(i,j)=-green(i,j)
        end if
     end do
  end do
  
end subroutine greens


subroutine hamiltonian(greenup,greendn,etot)
  double precision,dimension(:,:),intent(in)::greenup,greendn
  double precision,intent(out)::etot
  double precision::kinetic,potential
  double precision,dimension(:),allocatable::d_occ
  !double precision::kinetic,potential
  integer::i,i1,i2,j
! print*,"Calculating the expectation of the Hamiltonian"

 kinetic=0.0d0
 do i=1,NBND
 i1=ilbond(1,i)
 i2=ilbond(2,i)
 kinetic = kinetic + ARR(i1,i2)*(greenup(i1,i2)+greendn(i1,i2)+greenup(i2,i1)+greendn(i2,i1))
 end do
 

 potential=0.0d0

 do i=1,N
 potential = potential + U(i)*greenup(i,i)*greendn(i,i)
 !print*,potential,greenup(i,i),greendn(i,i)
 end do

 !print*,"KE:",kinetic,"PE:",potential
 etot=kinetic+potential

 !print*,'etot',etot
 end subroutine hamiltonian

 subroutine modgs(phi,ne,a)
 !modified Gramm-Schmidt orthogonalization
  integer,intent(in) :: ne
  real(kind=8),intent(inout) :: phi(N,*)
  real(kind=8),intent(out) :: a
  integer :: i,k
  real(kind=8) :: temp

  a=1.0d0
  do i=1,ne
     temp=0.0d0
     temp=ddot(N,phi(1,i),1,phi(1,i),1)
     temp=1.0d0/dsqrt(temp)
     a=a*temp
     call dscal(N,temp,phi(1,i),1)
     do k=i+1,ne
        temp=0.0d0
        temp=ddot(N,phi(1,i),1,phi(1,k),1)
        call daxpy(N,-temp,phi(1,i),1,phi(1,k),1)
     enddo
  enddo
 end subroutine modgs


end module util
