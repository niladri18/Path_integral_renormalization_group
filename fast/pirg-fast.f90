module pirg


     use setupmod
     use variable
     use lapack
     use inv
     use util
     use energy

        implicit none
        PRIVATE
        PUBLIC :: pirg_run

contains

subroutine pirg_run
!****************************************************************
!      Contains annealing and call of main PIRG projection
!****************************************************************
double precision :: esave


   NE=Nup + Ndn
   LWORK=3*N-1
   LWORK1=3*L-1

 call allocation 

 call write_setup
 
        f=0.0


        f0=ARR
  call keproject(f0)
  print*,'Size of KPROJECT',size(KPROJECT,1),size(KPROJECT,2)
  !print,KPROJECT
  call make_basis
 
  call make_HS 
  do ANNEAL = 1,NA
  print*,'anneal',anneal
  call pirg_anneal_a
  end do
  write(OUNIT,*)'Completed annealing'
  !print*,'Completed annealing'

!do ITER = 1,10
! ITER = 10 
esave = EGND 
ITER = 0
do  
   call fast_update
   del = abs(esave-EGND)
   ITER = ITER + 1
   write(OUNIT,'(A,I5,A,F12.8)')'Energy at',ITER,'=',EGND
   if (del<conv) exit
end do
   write(OUNIT,'(A,F12.8)')'Ground state energy = ',EGND

! do i=1,L
! do j=1,L
! print*,i,j,HAMILT(i,j)
! end do
! end do


!  do    
!  esave = EGND
!  call fast_update 
!  DEL = dabs(EGND-esave)
!  if(DEL<CONV)exit
!  end do
! print*,'Final energy',EGND
end subroutine pirg_run



!*****************************************************************************************



subroutine pirg_anneal_a

double precision :: esave,x,elow,aa
integer :: field,i1,j1,k,site,bit,i,j
double precision,allocatable :: psiup(:,:),psidn(:,:),tmpup(:,:),tmpdn(:,:)
! modifies the basis states with random H-S fields
! and accepts if it is lower in energy

  if(ANNEAL == 1)then
  call egr(esave,EVEC)
  end if
  
  write(OUNIT,*)'Initial energy',esave
  


 
  allocate(psiup(N,Nup))  
  allocate(psidn(N,Ndn))
  allocate(tmpup(N,Nup))
  allocate(tmpdn(N,Ndn))
  
  do k = 1,L
        psiup= BASUP(:,:,k)
        psidn= BASDN(:,:,k)

        !KE        
    call dgemm('N','N',N,Nup,N,1.0d0,KPROJECT,N,BASUP(:,:,k),N,0.0d0,tmpup,N)
    !call dcopy(N*Nup,tmpup,1,BASUP(:,:,k),1)
    tmpup = BASUP(:,:,k) 
    call dgemm('N','N',N,Ndn,N,1.0d0,KPROJECT,N,BASDN(:,:,k),N,0.0d0,tmpdn,N)
    !call dcopy(N*Ndn,tmpdn,1,BASDN(:,:,k),1)
    tmpdn = BASDN(:,:,k) 
 
!print*,'Done KE'       
        !PE
    do site=1,N
        call random_number(x)
        if(x<0.5)then
         field=1
        else
         field=0
        end if
!print*,HS(site,field)
!print*,site,field,HS(site,field)
        call dscal(Nup,HS(site,field),BASUP(site,:,k),1)

        call dscal(Ndn,HS(site,ieor(field,1)),BASDN(site,:,k),1)
    end do
        
!print*,'Done PE'
        call egr(elow,EVEC)
!        print*,esave,elow
        if(elow<esave)then
        esave = elow
        else
        psiup = BASUP(:,:,k)
        psidn = BASDN(:,:,k)
 
        !call dcopy(N*Nup,psiup,1,BASUP(:,:,k),1)
        !call dcopy(N*Nup,psidn,1,BASDN(:,:,k),1) 
        end if
        
  end do
 EGND = esave

 deallocate(psiup)
 deallocate(psidn) 
 deallocate(tmpup)
 deallocate(tmpdn)

end subroutine pirg_anneal_a


subroutine fast_update

double precision,allocatable ::gup(:,:),gdn(:,:),gupp(:,:),gdnn(:,:),otmp(:,:),eigval(:)
double precision,allocatable ::gtilup(:,:),gtildn(:,:),htmp(:,:),psiup(:,:),psidn(:,:)
double precision,allocatable ::tmpup(:,:),tmpdn(:,:)
integer :: i,j,k,bit1,bit2,f,t,i1,j1
double precision :: elow,esave,x,y,emat,aa

 
 call egr(elow,EVEC)
 !elow = EGND
 !EGND = elow

                allocate(gup(N,N))
                allocate(gdn(N,N))
                allocate(gupp(N,N))
                allocate(gdnn(N,N))
                allocate(otmp(L,L))
                allocate(htmp(L,L))
                allocate(eigval(L))
                allocate(gtilup(N,N))
                allocate(gtildn(N,N))
                allocate(tmpup(N,Nup))
                allocate(tmpdn(N,Ndn))
                allocate(psiup(N,Nup))
                allocate(psidn(N,Ndn))

    do k=1,L
        call modgs(BASUP(1,1,k),NUP,aa)
        call modgs(BASDN(1,1,k),NDN,aa)
    end do

!print*,'enterd fast update'
 
  do bit1 = 1,L

                psiup= BASUP(:,:,bit1)
                psidn= BASDN(:,:,bit1)
                !call dcopy(N*N,GREENUP(1,1,bit1,bit2),1,gup,1)
                !call dcopy(N*N,GREENDN(1,1,bit1,bit2),1,gdn,1)
                !call dcopy(L*L,O,1,otmp,1)
        !KE       
!print*,'Size of BASUP',size(BASUP(:,:,bit1),1),size(BASUP(:,:,bit1),2) 
!print*,'Size of BASUP',size(BASUP(:,:,bit1),1)
!print*,'Size of Kproject', size(KPROJECT, 1), size(KPROJECT, 2)

!print*,'Done!'
                call dgemm('N','N',N,Nup,N,1.0d0,KPROJECT,N,BASUP(:,:,bit1),N,0.0d0,tmpup,N)
!print*,'Size of tmpup', size(tmpup,1),size(tmpup,2)
                !call dcopy(N*Nup,tmpup,1,BASUP(:,:,bit1),1)
                BASUP(:,:,bit1) = tmpup
                call dgemm('N','N',N,Ndn,N,1.0d0,KPROJECT,N,BASDN(:,:,bit1),N,0.0d0,tmpdn,N)
                !call dcopy(N*Ndn,tmpdn,1,BASDN(:,:,bit1),1)
                BASDN(:,:,bit1) = tmpdn
                call egr(elow,EVEC)
!                if(elow<EGND)then
!print*,'Kinetic energy projection done!'

        !PE
        do site = 1,N

         do f = 0,1
                t=ieor(f,1)

                do bit2 = 1,L

               ! call dcopy(N*N,GREENUP(1,1,bit1,bit2),1,gup,1)
               ! call dcopy(N*N,GREENDN(1,1,bit1,bit2),1,gdn,1)
               ! call dcopy(L*L,O,1,otmp,1)

                gup = GREENUP(:,:,bit1,bit2) 
                gdn = GREENDN(:,:,bit1,bit2) 
                otmp = O
                !Calculation of G-tildae
                if(bit2==bit1)then
                        x = -(HS(site,f)-1)/(1.0d0+gup(site,site)*(HS(site,f)-1))
                        y = -(HS(site,t)-1)/(1.0d0+gup(site,site)*(HS(site,t)-1))
 
                        !calculating G-tildae        
                       ! call dcopy(N*N,gup,1,gtilup,1)
                        gtilup = gup
                        call dger(N,N,x,gup(:,site),1,gup(site,:),N,&
                            gtilup,N)
                        !print*,'Gtil',gtilup
                        !print*,'Size of gtil',size(gtilup,1),size(gtilup,2)
                        call dscal(N,HS(site,f),gtilup(:,site),1)
                        


                        call dcopy(N*N,gdn,1,gtildn,1)
                        call dger(N,N,x,gdn(1,site),1,gdn(site,1),N,&
                            gtilup,N)
                        call dscal(N,HS(site,t),gtildn(1,site),1)
                       
                        ! updating overlap
                        otmp(bit1,bit2) = O(bit1,bit2)*(HS(site,f))*(HS(site,t))*(HS(site,f))*(HS(site,t))*&
                                       gtilup(site,site)*gup(site,site)*gtildn(site,site)*gdn(site,site)
                        ! updating green's function


                        x = -(HS(site,f)-1)/(1.0d0+gtilup(site,site)*(HS(site,f)-1))               
                        y = -(HS(site,f)-1)/(1.0d0+gtildn(site,site)*(HS(site,f)-1))    
  
                        call dcopy(N*N,gtilup,1,gupp,1)
                        call dger(N,N,x,gupp(1,site),1,gupp(site,1),N,&
                                        gupp,N)
                        call dscal(N,HS(site,f),gupp(site,1),1)


                        call dcopy(N*N,gtildn,1,gdnn,1)
                        call dger(N,N,x,gdnn(1,site),1,gdnn(site,1),N,&
                                        gdnn,N)
                        call dscal(N,HS(site,t),gdnn(site,1),1)

                else 

                        x = -(HS(site,f)-1)/(1.0d0+gup(site,site)*(HS(site,f)-1))
                        y = -(HS(site,t)-1)/(1.0d0+gup(site,site)*(HS(site,t)-1))
                
                        ! updating overlap
                        otmp(bit2,bit1) = O(bit2,bit1)*(HS(site,f))*(HS(site,t))*&
                                       gup(site,site)*gdn(site,site)
                        otmp(bit1,bit2) = otmp(bit2,bit1)                  
                        ! updating green's function
                        
                        x = -(HS(site,f)-1)/(1.0d0+gup(site,site)*(HS(site,f)-1))
                        y = -(HS(site,f)-1)/(1.0d0+gdn(site,site)*(HS(site,f)-1))

                        call dcopy(N*N,gup,1,gupp,1)
                        call dger(N,N,x,gupp(1,site),1,gupp(site,1),N,&
                                        gupp,N)
                        call dscal(N,HS(site,f),gupp(1,site),1)



                        call dcopy(N*N,gdn,1,gdnn,1)
                        call dger(N,N,x,gdnn(1,site),1,gdnn(site,1),N,&
                                        gdnn,N)
                        call dscal(N,HS(site,f),gdnn(1,site),1)

                        
                end if


                call dcopy(L*L,HAMILT,1,htmp,1)
        
!do i1=1,L
!do j1=1,L
!print*,i1,j1,otmp(i1,j1),htmp(i1,j1)
!end do
!end do        
!print*,'Entered fast update'

!print*,'copied H to htmp'
                call hamiltonian(gupp,gdnn,emat)
!print*,emat,bit1,bit2,otmp(bit2,bit1)
                htmp(bit2,bit1) = emat*otmp(bit2,bit1)
                htmp(bit1,bit2) = htmp(bit2,bit1)
                !print*,'H-matrix elements',bit1,bit2,htmp(bit1,bit2)

!print*,'HTMP'
!                print*,htmp                !print*,'Copied overlaps to otmp'
                call DSYGV( 1, 'V' , 'U' , L , htmp , L , otmp , L , eigval , WORK1, LWORK1, INFO )
!print*,'esave',eigval(1),elow
                if(eigval(1)<elow)then
                        
                        esave = eigval(1)
                        elow = esave
!                        print*,'Eigval',eigval(1)
                        call dcopy(N*N,gupp,1,GREENUP(1,1,bit1,bit2),1)
!print*,'Copied Greens function'
                        call dcopy(N*N,gdnn,1,GREENDN(1,1,bit1,bit2),1)
                        call dcopy(L*L,htmp,1,HAMILT,1)
                        call dcopy(L*L,otmp,1,O,1)
                else 
                cycle
                end if
                !print*,'Eigval',eigval(1)
                EGND = esave
!                write(OUNIT,*)'Energy at iteration step ', ITER
!                write(OUNIT,*)EGND


                end do
         end do
       end do
  end do

                deallocate(gup)
                deallocate(gdn)
                deallocate(gupp)
                deallocate(gdnn)
                deallocate(gtilup)
                deallocate(gtildn)
                deallocate(otmp)
                deallocate(htmp)
                deallocate(eigval)




end subroutine fast_update
!******************************************************************************************

subroutine make_HS
integer :: i,j
double precision :: x,a,s 

 do i = 0,1
        do j=1,N
                x = dsqrt(dtanh(TAU*U(j)*0.25d0))
                a = 0.5d0*log((1+x)/(1-x))
                s=i*2-1
                HS(j,i)=dexp(2.0d0*a*s - TAU*U(j)*0.5d0)
        end do
 end do

end subroutine make_HS 
                
!*****************************************************************************************

subroutine write_setup

integer :: i,j,k

! writes the input on the output file

write(OUNIT,*)"Number of electrons",ne
 write(OUNIT,*)"Total no of bonds",nbnd
 write(OUNIT,*)"The Hopping terms"
 do i=1,N
 do j=1,N
        if((ARR(i,j).ne.0.0).and.(i<j))then
        write(OUNIT,'(2I5,F12.5)')i,j,ARR(i,j)
        end if
 end do
 end do

 write(OUNIT,*)'Hubbard U'
 do i=1,N
 write(OUNIT,'(I5,F12.5)')i,U(i)
 end do




end subroutine write_setup
!*********************************************************************************************88

subroutine write_wf
integer :: i,j,k



 !write(BUNIT)N
 !write(BUNIT)Nup
 !write(BUNIT)Ndn
 !write(BUNIT)L 
 do k=1,L

 do i=1,N

 do j=1,Nup
 write(BUNIT)BASUP(i,j,K)
! print*,K,i,j,basup(i,j,k)
 end do
 end do
 end do

 do k=1,L
 do i=1,N
 do j=1,Ndn
 write(BUNIT)BASDN(i,j,K)
! print*,K,i,j,basdn(i,j,k)
 end do
 end do
 end do



end subroutine write_wf



!**************************************************************************************




subroutine make_basis
integer :: i,j,k
double precision :: ran,aa


 do k=1,L
 do j=1,Ne
 do i=1,N
 call random_number(ran)
 BAS(i,j,k)=1.0d0*ran-0.50d0
 end do
 end do
 end do

 do k=1,L
 call modgs(BAS(:,:,k),Ne,aa)
 end do

 
 do k=1,L
            do j=1,Nup
                do i=1,N
                BASUP(i,j,k)=BAS(i,j,k)
                end do
            end do
 end do

 do k=1,L
            do j=1,Ndn
                do i=1,N
                BASDN(i,j,k)=BAS(i,j+Nup,k)
                end do
            end do
 end do
end subroutine make_basis




!***************************************************************************************



subroutine keproject(f0)
!for calculating the KE-projector
double precision,dimension(:,:),intent(inout)::f0
double precision,dimension(:),allocatable::tmp
double precision,dimension(:,:),allocatable::kexp,kpart,v_n

integer :: i,j,info



 print*,'Now diagonalizing',N
 allocate(tmp(N))        
 call DSYEV('V','U', N, f0, N, tmp, WORK, LWORK, info )
 
 
! do i=1,n
! print*,i,tmp(i)
! end do

 allocate(v_n(N,N))
 
 v_n=f0

! print*,"Calculating the Kinectic energy projector"

 allocate(kexp(1:N,1:N))
 kexp=0.0d0
 do i=1,N
  kexp(i,i)=dexp((-TAU*tmp(i))/1.0d0)
 end do

 allocate(kpart(N,N))
 call dgemm('N','N',N,N,N,1.0d0,v_n,N,kexp,N,0.0d0,kpart,N)

 call dgemm('N','T',N,N,N,1.0d0,kpart,N,v_n,N,0.0d0,KPROJECT,N)
 
 deallocate(kexp)
 deallocate(kpart)
 deallocate(v_n)
 

end subroutine keproject





end module pirg 
