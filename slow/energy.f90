module energy
 
 use variable
 use util
 use inv
 
 implicit none
 private
 public::eground

contains

subroutine eground(basup,basdn,eham,evec,L,Nup,Ndn,N)
  double precision,dimension(:,:,:),intent(inout)::basup,basdn
  double precision,dimension(:,:),intent(out)::evec
  double precision,intent(out)::eham
  integer,intent(in)::L,Nup,Ndn,N
  integer::bit1,bit2
  double precision,dimension(:,:,:),allocatable::tbasup,tbasdn
  double precision::en,detup,detdn,kinetic,potential,aa
  double precision,dimension(:,:),allocatable::hamilt,htmp,o
  double precision,dimension(:),allocatable::eigval
  allocate(hamilt(L,L))
  allocate(htmp(L,L))
  allocate(o(L,L))
  allocate(eigval(L))
  allocate(tbasup(N,Nup,L))
  allocate(tbasdn(N,Ndn,L))

  do bit1=1,L

     do bit2=1,L

        call greens(basup(:,:,bit1),basup(:,:,bit2),greenup,green2up,N,nup,detup)
        call greens(basdn(:,:,bit1),basdn(:,:,bit2),greendn,green2dn,N,ndn,detdn)
        call hamiltonian(greenup,greendn,en,kinetic,potential,N)
        hamilt(bit1,bit2)=en*detup*detdn
        o(bit1,bit2)=detup*detdn
     end do
  end do
  htmp=hamilt
  call DSYGV( 1, 'V' , 'U' , L , htmp , L , o , L , eigval , WORK1, LWORK1, INFO )
  if (info/=0) then
     write(*,*) 'info=',info,' in dsygv'
  endif
  eham=eigval(1)
  evec=htmp
! print*,"eham",eham  
end subroutine eground

end module energy
