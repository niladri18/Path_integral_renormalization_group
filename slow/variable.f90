module variable

! define all the global variables


implicit none

 !character::V,U,L
double precision,dimension(:,:,:),allocatable::basup,basdn,bas,tbasup,tbasdn
double precision,dimension(:,:),allocatable::arr,f,v_n,f0,v_ij,h_s,vtmp,vsc,basis,prod,overlap,overlapinv,green1,green2,green,charge,spin,newbasis,kexp,o,h0,htmp,eigvec,vec
double precision,dimension(:,:),allocatable::kpart,kproject,basisup,basisdn,greenup,greendn,green1up,green1dn,green2up,green2dn,hamilt,basis2up,basis2dn,tbas1up,tbas1dn,tbas2up,tbas2dn
integer,dimension(:,:),allocatable::s_index,SYM_TABLE
double precision,dimension(:),allocatable::d,p_on,d_hf,x_pos,y_pos,z_pos,eigval,eps,U,z,d_n,tmp,tsc,work,rhoup,rhodn,rho,d_occ,norm,uprojectup,uprojectdn,work1,sig
double precision::sum1,sum2,kinetic,potential,del,detup,detdn,xx
double precision::e_save,alphaup,alphadn,a,field,ran,conv,e1,aa,eup,edn,e_t,e_tot
double precision::omega,det,tau,hsq,var
integer::ne,n,i,i1,i2,j,k,k1,l,x,xc,it,k2,i3,i4,nn,no_states,nbonds1,nbonds2,nbonds3,nbonds4,nbonds5,numst,N_SYM
integer::ff,site,mu,nu,tstate,nbnd,hf,sulp,ii,lwork,info,v,nwork,iter,slater,sp,nup,ndn,lwork1,bit1,bit2,gbit,IUNIT,OUNIT,BUNIT,EVUNIT
integer,dimension(:,:),allocatable::ilbond

character(len=80)::FILENAME
logical::iexist

end module variable
