module measurement
 
 use variable
 use util
 use inv
 
 implicit none
 private
 public::ham2,measure

contains

subroutine corr(basup,basdn,rho,sig,d_occ,ch_ch,sp_sp,e_tot,h2,L,Nup,Ndn,N)
  double precision,dimension(:,:,:),intent(in)::basup,basdn
  double precision,dimension(:,:),intent(out)::ch_ch,sp_sp
  double precision,dimension(:),intent(out)::rho,sig,d_occ
  double precision,intent(out)::e_tot,h2
  integer,intent(in)::L,Nup,Ndn,N
  integer::bit1,bit2,i,j
  double precision,dimension(:,:,:),allocatable::tbasup,tbasdn
  double precision::en,detup,detdn,kinetic,potential,aa,numer,fac,hsq
  double precision,dimension(:,:),allocatable::tch_ch,tsp_sp
  double precision,dimension(:),allocatable::trho,tsig,td_occ
 
 allocate(trho(N))
 allocate(tsig(N))
 allocate(td_occ(N))
 allocate(tch_ch(N,N))
 allocate(tsp_sp(N,N))

  rho=0.0d0
  sig=0.0d0
  ch_ch=0.0d0
  sp_sp=0.0d0
 

  numer=0.0d0
  e_tot=0.0d0
  hsq=0.0d0
  do bit1=1,L

     do bit2=1,L

        call greens(basup(:,:,bit1),basup(:,:,bit2),greenup,green2up,N,nup,detup)
        call greens(basdn(:,:,bit1),basdn(:,:,bit2),greendn,green2dn,N,ndn,detdn)

        call ham2(greenup,green2up,greendn,green2dn,h2,N)
        call measure(greenup,greendn,trho,tsig,td_occ,tch_ch,tsp_sp,N)
        fac=detup*detdn
        do i=1,N
         rho(i)=fac*eigvec(bit1,1)*eigvec(bit2,1)*trho(i)+rho(i)
         sig(i)=fac*eigvec(bit1,1)*eigvec(bit2,1)*tsig(i)+sig(i)
         d_occ(i)=fac*eigvec(bit1,1)*eigvec(bit2,1)*td_occ(i)+d_occ(i)
        end do
        do i=1,N
        do j=1,N
         ch_ch(i,j)=fac*eigvec(bit1,1)*eigvec(bit2,1)*tch_ch(i,j)+ch_ch(i,j)
         sp_sp(i,j)=fac*eigvec(bit1,1)*eigvec(bit2,1)*tsp_sp(i,j)+sp_sp(i,j)
        end do
        end do
        call hamiltonian(greenup,greendn,en,kinetic,potential,N)
        e_tot=e_tot+en*eigvec(bit1,1)*eigvec(bit2,1)*detup*detdn
        hsq=hsq+h2*eigvec(bit1,1)*eigvec(bit2,1)*detup*detdn
        numer=numer+eigvec(bit1,1)*eigvec(bit2,1)*detup*detdn
        !print*,bit1,bit2,hamilt(bit1,bit2),o(bit1,bit2)
     end do
  end do
 
   rho=rho/numer
   sig=sig/numer
   d_occ=d_occ/numer
   ch_ch=ch_ch/numer
   sp_sp=sp_sp/numer
 
  e_tot=e_tot/numer
  hsq=hsq/numer
  h2=hsq

 end subroutine corr

 
 subroutine measure(greenup,greendn,rho,sig,docc,ch_ch,sp_sp,N)
 double precision,dimension(:,:),intent(in)::greenup,greendn
 double precision,dimension(:,:),intent(out)::ch_ch,sp_sp
 double precision,dimension(:),intent(out)::rho,sig,docc
 integer,intent(in)::N
 double precision,dimension(:),allocatable::rhoup,rhodn
 double precision::sum1,sum2
 allocate(rhoup(N))
 allocate(rhodn(N))
 do i=1,N
 rhoup(i)=greenup(i,i)
 rhodn(i)=greendn(i,i)
 rho(i)=greenup(i,i)+greendn(i,i)
 sig(i)=greenup(i,i)-greendn(i,i)
 docc(i)=greenup(i,i)*greendn(i,i)
 end do

 do i=1,n
 do j=1,n
   sum1=rhoup(i)*rhoup(j) + greenup(i,j)*green2up(i,j)
   sum2=rhoup(i)*rhodn(j) + rhodn(i)*rhoup(j)
   ch_ch(i,j)=2.0d0*sum1 + sum2
   sp_sp(i,j)=2.0d0*sum1 - sum2
   !write(21,*)i,j,charge(i,j),spin(i,j)
 end do
 end do
 end subroutine measure

 subroutine ham2(greenup,green2up,greendn,green2dn,h2,N)
 double precision,dimension(:,:),intent(in)::greenup,green2up,greendn,green2dn
 double precision,intent(out)::h2
 integer,intent(in)::N
 double precision::h2a,h2b,h2c
 integer::i,j,i1,j1,i2,j2

 ! U^2 term
  h2a=0.0d0
  do i=1,N
     do j=1,N
        h2a=h2a+U(i)*U(j)*(&
             greenup(i,i)*greenup(j,j)*greendn(i,i)*greendn(j,j)&
             +greenup(j,i)*green2up(i,j)*greendn(i,i)*greendn(j,j)&
             +greenup(i,i)*greenup(j,j)*greendn(j,i)*green2dn(i,j)&
             +greenup(j,i)*green2up(i,j)*greendn(j,i)*green2dn(i,j))
     enddo

  enddo
  


 
   ! t^2 term
  h2b=0.0d0
  do i=1,NBND
     i1=ilbond(1,i)
     i2=ilbond(2,i)
     do j=1,NBND
        j1=ilbond(1,j)
        j2=ilbond(2,j)
        h2b=h2b+f0(i1,i2)*f0(j1,j2)*(greenup(i1,i2)*greenup(j1,j2)+greenup(j2,i1)*green2up(i2,j1)&
             +greenup(i1,i2)*greendn(j1,j2)&
             +greendn(i1,i2)*greenup(j1,j2)&
             +greendn(i1,i2)*greendn(j1,j2)+greendn(j2,i1)*green2dn(i2,j1)&

             +greenup(i2,i1)*greenup(j1,j2)+greenup(j2,i2)*green2up(i1,j1)&
             +greenup(i2,i1)*greendn(j1,j2)&
             +greendn(i2,i1)*greenup(j1,j2)&
             +greendn(i2,i1)*greendn(j1,j2)+greendn(j2,i2)*green2dn(i1,j1)&

             +greenup(i1,i2)*greenup(j2,j1)+greenup(j1,i1)*green2up(i2,j2)&
             +greenup(i1,i2)*greendn(j2,j1)&
             +greendn(i1,i2)*greenup(j2,j1)&
             +greendn(i1,i2)*greendn(j2,j1)+greendn(j1,i1)*green2dn(i2,j2)&

             +greenup(i2,i1)*greenup(j2,j1)+greenup(j1,i2)*green2up(i1,j2)&
             +greenup(i2,i1)*greendn(j2,j1)&
             +greendn(i2,i1)*greenup(j2,j1)&
             +greendn(i2,i1)*greendn(j2,j1)+greendn(j1,i2)*green2dn(i1,j2))
     enddo
  enddo


 ! tU term
  h2c=0.0d0
  do i=1,NBND
     i1=ilbond(1,i)
     i2=ilbond(2,i)
     do k=1,N
        h2c=h2c+U(k)*f0(i1,i2)*2.0d0*( &
              greenup(i2,i1)*greenup(k,k)*greendn(k,k) + greenup(k,i1)*green2up(i2,k)*greendn(k,k) &
             +greenup(i1,i2)*greenup(k,k)*greendn(k,k) &
             +greenup(k,i2)*green2up(i1,k)*greendn(k,k) &
             +greendn(i2,i1)*greendn(k,k)*greenup(k,k) &
             +greendn(k,i1)*green2dn(i2,k)*greenup(k,k) &
             +greendn(i1,i2)*greendn(k,k)*greenup(k,k) &
             +greendn(k,i2)*green2dn(i1,k)*greenup(k,k) )
       enddo
     enddo
 !print*,h2a,h2b,h2c
   h2=h2a+h2b+h2c

 end subroutine ham2
end module measurement
