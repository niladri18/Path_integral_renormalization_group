program noninteracting

     use variable
     use lapack
     use inv
     use util
     use energy
     use measurement

 tau=0.01
 IUNIT=91
 OUNIT=22
 BUNIT=77  !for the binary file 
 EVUNIT=45 !for the coefficient of the eigen-vector

call get_command_argument(1,FILENAME)

    inquire(file=FILENAME,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'File not found'
     stop
  endif

  print*,FILENAME
  open(unit=IUNIT, file=FILENAME, action="read")
  open(unit=OUNIT, file=trim(FILENAME)//'.out', status='replace', &
          action='write')
  open(unit = BUNIT, status='replace',file=trim(FILENAME)//'.bin',form='unformatted')

 open(unit=EVUNIT, file=trim(FILENAME)//'.wgt',status='replace', &
                        action = 'write')
  ! open(unit=91,file="1d.inp")
   open(unit=21,file="greens.out")
   read(IUNIT,*) n,nup,ndn,L,slater,conv
   

	ne=nup+ndn
   lwork=3*n-1 
   lwork1=3*L-1
	allocate(arr(1:n,1:n))
	allocate(f(1:n,1:n))
	allocate(f0(1:n,1:n))
	allocate(U(1:n)) 
	allocate(v_n(1:n,1:n))
	allocate(d(1:n))
	allocate(d_n(1:n))
	allocate(tmp(1:n))
	allocate(work(1:lwork))
	allocate(work1(1:lwork1))
 
 print*,"lwork",lwork,lwork1

allocate(ilbond(2,200))

 !do i=1,2
! do j=1,200
! ilbond(i,j)=0.0
! end do
! end do


! ====== Reading Main File ====
 
 !read(91,*) n
 read(IUNIT,*) nbnd 
 print*,"total no of bonds",nbnd
 !Define this array variable: ilbnd(2,100)
 read(IUNIT,*) (ilbond(1,k),ilbond(2,k),arr(ilbond(1,k),ilbond(2,k)),arr(ilbond(2,k),ilbond(1,k)),k=1,nbnd)
 !read(91,*)(arr(ilbond(2,k)),arr(ilbond(2,k),ilbond(1,k))k=1,nbnd)
 !write(*,*) 'links',(ilbond(1,k),ilbond(2,k),k=1,nbnd)
 !Define siten, Hub-U, and chem-pot arrays
 !read(91,*) (eps(i),i=1,n)
 read(IUNIT,*) (U(i),i=1,n)


 print*,"Number of electrons",ne
 
 !arr=-arr
 write(OUNIT,*)"Number of electrons",ne
 write(OUNIT,*)"total no of bonds",nbnd


 do i=1,n
 do j=1,n
      if((arr(i,j).ne.0.0).and.(i<j))then
       write(OUNIT,*)i,j,arr(i,j)
      end if
 end do
 end do




  


 f=0.0
!v_ij=0.0



 f0=arr




 call DSYEV('V','U', N, arr, N, tmp, WORK, lwork, info )
 
 
 do i=1,n
 write(OUNIT,*)i,U(i)
 end do



 print*,"Huckel calculation ends"
 
v_n=arr
 print*,"Calculating the Kinectic energy projector"

 allocate(kexp(1:N,1:N))
 kexp=0.0d0
 do i=1,N
  kexp(i,i)=exp(-tau*tmp(i)/1.0d0)
 end do

 allocate(kpart(N,N))
 call dgemm('N','N',N,N,N,1.0d0,v_n,N,kexp,N,0.0d0,kpart,n)

 allocate(kproject(N,N))
 call dgemm('N','T',N,N,N,1.0d0,kpart,N,v_n,N,0.0d0,kproject,n)
 
!do i=1,n
!do j=1,n
!  print*,i,j,kproject(i,j)
!end do
!end do 
 
 !! ========== allocation section ==========
 allocate(green(N,N))
 allocate(greenup(N,N))
 allocate(greendn(N,N)) 
 allocate(green2up(N,N))
 allocate(green2dn(N,N))
 allocate(d_occ(n))
 allocate(rho(n))
 allocate(sig(n))
 allocate(rhoup(n))
 allocate(rhodn(n))
 allocate(charge(n,n))
 allocate(spin(n,n))
 allocate(uprojectup(n))
 allocate(uprojectdn(n))
 allocate(eigvec(L,L))
 allocate(vec(L,L))
 !constructs the slater basis here
 allocate(bas(N,Ne,L))
 do k=1,L
 do j=1,Ne
 do i=1,N
 call random_number(ran)
 bas(i,j,k)=2.0d0*ran-1.0d0
 end do
 end do 
 end do

 !print*,"Hi!"
 !if (slater==0)then
 do k=1,L
 call modgs(bas(:,:,k),Ne,aa)
 end do
 !print*,"Hi"
 
 call modgs(bas,Ne,aa)

 allocate(basup(N,nup,L))
 allocate(basdn(N,ndn,L))
 allocate(norm(Nup))
 do k=1,L
 do j=1,Nup
 norm(j)=0.0d0
 do i=1,N
 basup(i,j,k)=bas(i,j,k)
 norm(j)=norm(j)+basup(i,j,k)**2
 end do
 do i=1,N
 basup(i,j,k)=basup(i,j,k)/norm(j)
 end do
 end do
 end do
 
 deallocate(norm)
 allocate(norm(Ndn))
 do k=1,L
 do j=1,Ndn
 norm(j)=0.0d0
 do i=1,N
 basdn(i,j,k)=bas(i,j+Nup,k)
 norm(j)=norm(j)+basdn(i,j,k)**2
 end do
 !print*,norm(j)
 do i=1,N
 basdn(i,j,k)=basdn(i,j,k)/norm(j)
 end do
 end do
 end do


    do i=1,L
       call modgs(basup(:,:,i),Nup,aa)
       call modgs(basdn(:,:,i),Ndn,aa)
    enddo 

 call eground(basup,basdn,e_save,vec,L,Nup,Ndn,N)

 print*,e_save
 write(OUNIT,*)'Initial energy',e_save

 !stop

 ! ====== end calculating the initial energy ====== !
 
 







! ================starting the iteration process ===================== !







 iter=0
 allocate(tbasup(N,Nup,L))
 allocate(tbasdn(N,Ndn,L))
 print*,"Starting the iteration"
 write(OUNIT,*)"Starting the iteration"



 do 
    e1=e_save
    iter=iter+1

 !if(iter>143)stop
! if(iter==1)then
! print*,o
! end if

    ! ====== loop for each basis =======!

    do bit1=1,L
    gbit=bit1
       ! ======= the interaction part ====== !


!       do ff=0,1 ! each HS field

       do site = 1,n
          
          xx = (dsqrt(dtanh(tau*U(site)/4.0d0)))
          a = 0.5d0*dlog((1+xx)/(1-xx)) 
          !field=-1.d0 ! each of the field
        


          do ff=0,1
             field=2*ff-1.d0
             alphaup = 2.d0*a*field-(tau*U(site)/2.0d0)
             alphadn = -2.d0*a*field-(tau*U(site)/2.0d0)
             !uproject=0.0d0
             uprojectup(site) = dexp(alphaup)
             uprojectdn(site) = dexp(alphadn)
          
   
             tbasup=basup
             tbasdn=basdn 
             do i=1,Nup
                tbasup(site,i,bit1)=uprojectup(site)*tbasup(site,i,bit1)
                !print*,'tbasup',tbasup(site,i,bit1) 
             end do
             
             do i=1,Ndn
                tbasdn(site,i,bit1)=uprojectdn(site)*tbasdn(site,i,bit1)
                !print*,'tbasdn',tbasdn(site,i,bit1)
             end do

    do i=1,L
       call modgs(tbasup(:,:,i),Nup,aa)
       call modgs(tbasdn(:,:,i),Ndn,aa)
    enddo

!*******************************************************************

!       if(iter==140)then
!        print*,'before',site,bit1
!        do bit2=1,L
!        print*,bit2
!        print*,basup(:,:,bit2)!,tgreendn
!        end do
!       end if
!
!****************************************************************

              
             call eground(tbasup,tbasdn,e_t,vec,L,Nup,Ndn,N) 


!           print*,'From slow update',e_t,site,bit1,iter
!       if(iter==140)then
!        print*,'after',site,bit1
!        do bit2=1,L
!        print*,bit2
!        print*,tbasup(:,:,bit2)!,tgreendn(:,:,bit2)
!        end do
!       end if






!**************************************************
!						  *
 
 !      if(iter==140)then
          !print*,e_t,e_save
 !         if(e_t<e_save)then
 !            print*,'accepted',e_t,e_save,site,bit1
 !         else
 !            print*,'not accepted',e_t,e_save,site,bit1
 !         end if
 !       end if

!						  *
!**************************************************

             if(e_t<e_save)then
                basup=tbasup
                basdn=tbasdn
                e_save=e_t
                EIGVEC=vec
  !  do i=1,L
  !     call modgs(basup(:,:,i),Nup,aa)
  !     call modgs(basdn(:,:,i),Ndn,aa)
  !  enddo
         !    else
        !     cycle
             end if

     ! if(iter==140)then
     !   print*,'after',site,bit1
     !   do bit2=1,L
     !   print*,bit2
     !   print*,basup(:,:,bit2)!,tgreendn(:,:,bit2)
     !   end do
     !  end if


!             field=field+2.d0
             
          end do !for the site
       end do ! for each HS field
 !    print*,"=============bit==========",bit1
!      print*, 'PE',e_t,iter,e_save
 ! =========== the kinetic energy part ========== !

 !  if(iter==140)then
 !   print*,'leaving the U-projection loop'
 ! end if


       tbasup=basup
       tbasdn=basdn
       
       allocate(newbasis(N,nup))
       call dgemm('N','N',N,nup,N,1.0d0,kproject,N,basup(:,:,bit1),N,0.0d0,newbasis,N)
       tbasup(:,:,bit1) = newbasis
       deallocate(newbasis)
       allocate(newbasis(N,ndn))
       call dgemm('N','N',N,ndn,N,1.0d0,kproject,N,basdn(:,:,bit1),N,0.0d0,newbasis,N)
       tbasdn(:,:,bit1) = newbasis
       deallocate(newbasis)

    do i=1,L
       call modgs(tbasup(:,:,i),Nup,aa)
       call modgs(tbasdn(:,:,i),Ndn,aa)
    enddo
       
       call  eground(tbasup,tbasdn,e_t,vec,L,Nup,Ndn,N)

       if(e_t<e_save)then
          basup=tbasup
          basdn=tbasdn
          e_save=e_t
          EIGVEC=vec
!   do i=1,L
!      call modgs(basup(:,:,i),Nup,aa)
!      call modgs(basdn(:,:,i),Ndn,aa)
!   enddo
!      else
      cycle 
      end if
   !    print*, 'KE',e_save
    enddo

    !do i=1,L
    !   call modgs(basup(:,:,i),Nup,aa)
    !   call modgs(basdn(:,:,i),Ndn,aa)
    !enddo

    del=dabs(e_save-e1)
!    del=(e_save-e1)**2
     
    print*,"step",iter,e_save
    write(22,*)"step",iter,e_save
    if(del.le.conv)exit
 end do !iteration ends

! do K=1,L
! write(78,*)'basis no',k
! write(78,*),'#upspin'
! do i=1,N
! write(78,*)'#site',i
! do j=1,Nup
! write(78,*)basup(i,j,K)
! end do
! end do
! end do 
! write(78,*)'#dnspin'
! do K=1,L
! write(78,*)'basis no',k
! do i=1,N
! write(78,*)'#site',i
! do j=1,Ndn
! write(78,*)basdn(i,j,K)
! end do
! end do
! end do


 do K=1,L

 do i=1,N

 do j=1,Nup
 write(BUNIT)basup(i,j,K)
! print*,K,i,j,basup(i,j,k)
 end do
 end do
 end do

 do k=1,L
 do i=1,N
 do j=1,Ndn
 write(BUNIT)basdn(i,j,K)
! print*,K,i,j,basdn(i,j,k)
 end do
 end do
 end do


 print*,"Energy",e_save,"Converged in steps",iter
 write(OUNIT,*)"Energy",e_save,"Converged in steps",iter

 do i=1,L
 write(EVUNIT,*)EIGVEC(i,1)
 end do

!do i=1,N
!print*,uprojectdn(i)
!end do


 ! ============The measurements========== !

 call corr(basup,basdn,rho,sig,d_occ,charge,spin,e_tot,hsq,L,Nup,Ndn,N)

  var=(hsq-e_tot**2)/(e_tot**2)  
 print*,"Mean energy",e_tot,var
 write(OUNIT,*)"Mean energy",e_tot,var 
 print*,"Squared H",hsq
        write(OUNIT,*)"Bond order"
  do i=1,n-1
        write(OUNIT,*)i,i+1,greenup(i,i+1)+greendn(i,i+1)
  end do

! calculating the charge densities
 write(OUNIT,*)"Charge density and spin density"
 do i=1,n
 write(OUNIT,*)i,rho(i),sig(i)
 end do

! calculating the double occupancies

write(OUNIT,*)"Double occupancy"
do i=1,n
!d_occ(i)=greenup(i,i)*greendn(i,i)
!rhoup(i)=greenup(i,i)
!rhodn(i)=greendn(i,i)
!rho(i)=rhoup(i)+rhodn(i)
!write(21,*)"Double occupancy"
write(OUNIT,*)i,d_occ(i)!,rho(i)
end do


! calculates the charge charge correlation


   write(OUNIT,*)"Charge and spin correlation"
 do i=1,n
 do j=1,n
!   sum1=rhoup(i)*rhoup(j) + greenup(i,j)*green2up(i,j)
!   sum2=rhoup(i)*rhodn(j) + rhodn(i)*rhoup(j)
!   charge(i,j)=2.0d0*sum1 + sum2
!   spin(i,j)=2.0d0*sum1-sum2
   write(OUNIT,*)i,j,charge(i,j),spin(i,j)
 end do
 end do

 



end program noninteracting
 
