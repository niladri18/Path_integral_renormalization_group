module variable

! define all the global variables


implicit none

 !character::V,U,L
double precision,dimension(:,:,:,:),allocatable::GREENUP,GREENDN,GREEN2UP,GREEN2DN,kgreenup,kgreendn,kgreen2up,kgreen2dn
double precision,dimension(:,:,:),allocatable::BASUP,BASDN,BAS,NBASUP,NBASDN,tgreenup,tgreendn,tgreen2up,tgreen2dn
double precision,dimension(:,:),allocatable::ARR,F,v_n,F0,v_ij,h_s,vtmp,vsc,basis,prod,overlap,overlapinv,green1,green2,green,charge,spin,newbasis,kexp,O,h0,htmp,evec,vec,t_o,x_o
double precision,dimension(:,:),allocatable::KPROJECT,basisup,basisdn,green1up,green1dn,HAMILT,basis2up,basis2dn,tbas1up,tbas1dn,tbas2up,tbas2dn,ham,xham,HS
integer,dimension(:,:),allocatable::s_index
double precision,dimension(:),allocatable::d,p_on,d_hf,x_pos,y_pos,z_pos,eigval,eps,U,z,d_n,tmp,tsc,work,rhoup,rhodn,rho,d_occ,norm,uprojectup,uprojectdn,work1,sig
double precision::sum1,sum2,kinetic,potential,del,detup,detdn,xx
double precision::e_save,alphaup,alphadn,a,field,ran,conv,e1,aa,eup,edn,e_t,e_tot,e_x,aa1,aa2,esave
double precision::omega,det,tau,hsq,var,aup,adn,ent,EGND
integer::ne,i,i1,i2,j,k,k1,l,x,xc,it,k2,i3,i4,nn,no_states,nbonds1,nbonds2,nbonds3,nbonds4,nbonds5,numst,Nup,Ndn,N,NA,ANNEAL
integer::ff,site,mu,nu,tstate,nbnd,hf,sulp,ii,lwork,info,v,nwork,iter,slater,sp,lwork1,bit1,bit2,OUNIT,BUNIT,EVUNIT
integer,dimension(:,:),allocatable::ilbond

!character(len=80)::FILENAME
!logical::iexist


end module variable
