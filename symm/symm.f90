program symmetry

!use makefile to compile
     use mpi
     use variable
     use lapack
     use inv
     use util
     use energy
     use measurement



!*********************************************************************
!                                                                    *       
!                       Initialize MPI
!                                                                    *
!*********************************************************************
!                                                                    *

        call MPI_INIT ( ierr )

      
        !Get the number of processors this job is using
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)

        ! Get the rank of the processor this thread is running on.
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

!                                                                    *
!*********************************************************************

 IUNIT1=91
 IUNIT2=92
 IUNIT3=93


call get_command_argument(1,FILENAME)

    inquire(file=FILENAME,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'Input file not found'
     stop
  endif

  print*,FILENAME


call get_command_argument(2,FNAME)
 print*,FNAME
    inquire(file=FNAME,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'Symmetry file not found'
     stop
  endif

call get_command_argument(3,WF)

    inquire(file=WF,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'Wavefunction file not found'
     stop
  endif

open(unit=IUNIT1,file=FILENAME, action = 'read')
open(unit=IUNIT2, status = 'old', file = WF, form ='unformatted')
open(unit=IUNIT3,file=FNAME, action = 'read')
open(unit=OUNIT, file = trim(FILENAME)//'.corr', status = 'replace', &
                action = 'write')
!open(unit=93,file="translation")
 

read(IUNIT1,*) n,nup,ndn,L,slater,conv

print*,'Number of proc',numproc
write(OUNIT,*)'Number of proc',numproc

LWORK = 3*L-1

!*********************************************************************
!                                                                    *
!                       Allocation section I
!                                                                    *
!*********************************************************************
!                                                                    *
 allocate(basup(N,nup,L))
 allocate(basdn(N,ndn,L))

 allocate(tbasup(N,nup,L))
 allocate(tbasdn(N,ndn,L))

 allocate(hamilt(L,L))
 allocate(h_part(L,L))
 allocate(o(L,L))
 allocate(o_part(L,L))


 allocate(greenup(N,N))
 allocate(greendn(N,N))
 allocate(green2up(N,N))
 allocate(green2dn(N,N))




 allocate(WORK(1:LWORK))
 allocate(eigval(1:L))
 allocate(vec(L,L))
!                                                                    *
!*********************************************************************



 print*,n,nup,ndn

 pi = 4.0d0*datan(1.0d0)
 

!*********************************************************************
!                                                                    *
!                      Reading saved wavefunction
!                                                                    *
!*********************************************************************
!                                                                    *

 do K=1,L
 
 do i=1,N
 
 do j=1,Nup
 read(IUNIT2)basup(i,j,K)
! print*,K,i,j,basup(i,j,k)
 end do
 end do
 end do

 do k=1,L
 do i=1,N
 do j=1,Ndn
 read(IUNIT2)basdn(i,j,K)
! print*,K,i,j,basdn(i,j,k)
 end do
 end do
 end do

!***********************************************************************
!                           Reading symmetries 
!**********************************************************************

 read(IUNIT3,*)N_SYM
 allocate(SYM_TABLE(N_SYM,N))
 do k=1,N_SYM
 read(IUNIT3,*)(SYM_TABLE(i,k), i=1,N)
 end do

! print*,'In processor',rank,'the symmetry table is'
! do i=1,N
! do k=1,N_SYM
!  write(OUNIT,*),i,k,SYM_TABLE(i,k)
! end do
! end do

!                                                                    * 
!*********************************************************************







 read(IUNIT1,*) nbnd
 !print*,"total no of bonds",nbnd




!*********************************************************************
!                                                                    *
!                       Allocation section II
!                                                                    *
!*********************************************************************
!                                                                    *
        allocate(arr(1:n,1:n))
        allocate(f0(1:n,1:n))
        allocate(U(1:n))
        allocate(ilbond(2,nbnd))
!                                                                    *
!*********************************************************************











!*********************************************************************
!                    Reading the input file
!*********************************************************************
!                                                                    *

 read(IUNIT1,*)(ilbond(1,k),ilbond(2,k),arr(ilbond(1,k),ilbond(2,k)),arr(ilbond(2,k),ilbond(1,k)),k=1,nbnd)
 read(IUNIT1,*) (U(i),i=1,n)
!                                                                    *
!*********************************************************************


        f0=arr




if(rank==0)then
 
 write(OUNIT,*)'Bonds'
 do i=1,N
 do j=1,N
 if((i>j).and.arr(i,j).ne.0.0d0)then
  write(OUNIT,*)i,j,arr(i,j)
 end if
 end do
 end do

 
 write(OUNIT,*)'Hubbard U'
 do i=1,N
 write(OUNIT,*)i,U(i)
 end do

 
 write(OUNIT,*)'The symmetry table'
 do i=1,N_SYM
 do k=1,N
  write(OUNIT,*),i,k,SYM_TABLE(i,k)
 end do
 end do

end if
 
  !tbasup = basup
  !tbasdn = basdn

  count = N_SYM/numproc
  res = mod(N_SYM,numproc)
  
!**********************************************************
! Choosing the last CPU to put the residual entries there
!**********************************************************
!                                                         *
  if(rank==numproc)then
    mpires=1
  else
    mpires=0
  end if
!                                                         *
!**********************************************************


 ! if (rank.ne.0) then
 !       start = (rank-1)*count + 1
 !       finish = rank*count + mpires*res
 ! 
 ! end if
       start = rank*count + 1
       finish = (rank+1)*count + mpires*res
         momnt = 2.0d0*3.14d0/N

 !      print*,'start and finish count',start,finish,'in proc',rank
 !do k_val=0,N_SYM-1
 k_val = 0
 h_part=0.0
 o_part=0.0

 do bit1=1,L
   
  !if(rank==0)cycle
 print*,'I am processor',rank
 write(OUNIT,*)'I am processor',rank

  do symc=start,finish

 !       print*,start,finish
       !print*,'I am processor',rank 
        tbasup = basup
        tbasdn = basdn
        !print*,'Up-spin'
        do j=1,Nup
        do i=1,N
                
                kk = SYM_TABLE(i,symc)
!                print*,kk,symc                
                tbasup(kk,j,bit1)=basup(i,j,bit1)
                !print*,kk,j,bit1,sync
         !       print*,tbasup(kk,j,bit1)
        enddo
        enddo

        !print*,'Dn-spin'

        do j=1,Ndn
        do i=1,N
                kk = SYM_TABLE(i,symc)
                tbasdn(kk,j,bit1)=basdn(i,j,bit1)
                !print*,kk,j,bit1 
                !print*,tbasdn(kk,j,bit1)
        enddo
        enddo
       
  !      print*,'Translation done successfully in',rank
  write(OUNIT,*)'Symmetry done done successfully in',rank

        do bit2=1,L
!        print*,'entered loop',size(greendn),size(green2dn)
        call greens(tbasup(:,:,bit1),tbasup(:,:,bit2),greenup,green2up,N,nup,detup)
!        print*,'entered loop',size(greendn),size(green2dn)
 !       print*,nup,detup
        call greens(tbasdn(:,:,bit1),tbasdn(:,:,bit2),greendn,green2dn,N,ndn,detdn)
!        print*,'Hi! entered loop',size(greendn),size(green2dn)
  !      print*,ndn,detdn
        call hamiltonian(greenup,greendn,en,kinetic,potential,N)
        !print*,'After calculating Greens function',en
        e_t = en*detup*detdn*dcos(momnt*k_val*symc) 
   !     print*,en,detup,detdn,momnt
        h_part(bit1,bit2) = e_t + h_part(bit1,bit2)
        h_part(bit2,bit1) = h_part(bit1,bit2)
        o_t = detup*detdn*dcos(momnt*k_val*symc)
        o_part(bit1,bit2) = o_t + o_part(bit1,bit2)
        o_part(bit2,bit1) = o_part(bit1,bit2)

        end do

  end do ! end summing over the symmetries


  

 end do ! This is the loop for each basis


    ! do i=1,L
    !    do j=1,L
    !    o_part(i,j) = o_part(i,j)/N_sym
    !    h_part(i,j) = h_part(i,j)/N_sym
    !    end do
    ! end do
  !   do i=1,L
   !     do j=1,L
    !    print*,i,j,o_part(i,j),h_part(i,j)
     !   end do
     !end do
        print*,'In processor',rank,'the hamiltonian is'
        print*,h_part
        dim=L*L
 !if(rank/=0)then
  call MPI_REDUCE(h_part,hamilt,dim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(o_part,o,dim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 ! print*,'Finish processor in:',rank
 !else
  !call MPI_REDUCE(MPI_IN_PLACE,h_part,dim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  !call MPI_REDUCE(MPI_IN_PLACE,o_part,dim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! print*,'Finish processor in:',rank
 !end if

  hamilt = hamilt/N_SYM
  o = o/N_SYM

  !h_part = h_part/N_sym
  !o_part = o_part/N_sym
 if(rank==0)then
  print*,'full hamiltonian',hamilt
   write(OUNIT,*)'full Many Body hamiltonian'
   write(OUNIT,*)hamilt
 end if

  if (rank==0)then
        print*,'I am the Master',rank
         
          !call DSYGV( 1, 'V' , 'U' , L , h_part , L , o_part , L , eigval , WORK, LWORK, INFO )
          call DSYGV( 1, 'V' , 'U' , L , hamilt , L , o , L , eigval , WORK, LWORK, INFO )

          if(info/=0)then
                 print*,'info neq 0',info
          end if

        print*,'Energy for K=',k_val,'is',eigval(1)
        write(OUNIT,*)'Energy for K=',k_val,'is',eigval(1)

  end if
 !end do ! This is the loop for each k-value




 call MPI_Finalize(ierr)


end program symmetry



 





