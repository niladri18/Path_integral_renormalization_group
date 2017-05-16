module energy
 
 use variable
 use util
 use inv
 
 implicit none
 private
 public::egr

contains



subroutine egr(eham,evec)

! Calculates the matirx elements of the Hamiltonian with 
! current wavefunction and finds the lowest energy

  double precision,dimension(:,:),intent(out)::evec
  double precision,intent(out)::eham
!  integer,intent(in)::bit1!L,Nup,Ndn,N
  integer::bit1,bit2
  double precision,dimension(:,:,:),allocatable::tbasup,tbasdn
  double precision::en,detup,detdn,kinetic,potential,aa
  double precision,dimension(:,:),allocatable::htmp,otmp
  double precision,dimension(:),allocatable::eigval
  !allocate(hamilt(L,L))
  allocate(htmp(L,L))
  allocate(otmp(L,L))
  allocate(eigval(L))
  allocate(tbasup(N,Nup,L))
  allocate(tbasdn(N,Ndn,L))

  do bit1=1,L

     do bit2=1,L

        call greens(BASUP(:,:,bit1),BASUP(:,:,bit2),GREENUP(:,:,bit1,bit2),GREEN2UP(:,:,bit1,bit2),N,nup,detup)
        call greens(BASDN(:,:,bit1),BASDN(:,:,bit2),GREENDN(:,:,bit1,bit2),GREEN2DN(:,:,bit1,bit2),N,ndn,detdn)
        call hamiltonian(GREENUP(:,:,bit1,bit2),GREENDN(:,:,bit1,bit2),en)
        HAMILT(bit1,bit2)=en*detup*detdn
        O(bit1,bit2)=(detup*detdn)
        !print*,bit1,bit2,hamilt(bit1,bit2),o(bit1,bit2)
     end do
  end do
  htmp=HAMILT
  otmp=O
  !print*,otmp
  !write(*,*) 'info for 2nd DSYGV call: before',info
  call DSYGV( 1, 'V' , 'U' , L , htmp , L , otmp , L , eigval , WORK1, LWORK1, INFO )
  !write(*,*) 'info for 2nd DSYGV call: after',info
  if (info/=0) then
     write(*,*) 'info=',info,' in dsygv in egr'
  endif
  eham=eigval(1)
  evec=htmp
  print*,eham
  write(OUNIT,*)'Energy in EGR', eham
! print*,o
!  o=otmp
 

end subroutine egr

end module energy
