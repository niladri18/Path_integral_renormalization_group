module setupmod
        use variable

        implicit none
        PRIVATE
        PUBLIC :: read_config,FILENAME,allocation
        CHARACTER(len=80) :: FILENAME

contains

subroutine read_config

!******************************************************************
!                    Read the input files 
!******************************************************************




integer, parameter :: nnames = 8
integer, parameter :: IUNIT = 10
integer :: status,i

character(len=20), parameter :: names(1:nnames) = &
        (/'nsites        ','n_up         ','n_dn       ',&
          'lmax          ','conv        ' ,'bond_file ',&   
          'U_file        ','NA          '/)    

character(len=150) :: cmd, val
character(len=150) :: buffer
logical :: iexist


  inquire(file=FILENAME,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'File not found'
     stop
  endif

  print*,FILENAME
  open(unit=IUNIT, file=FILENAME, action="read")


    do
     read(iunit, '(a100)', IOSTAT=status) buffer
     !print*,buffer
     if (status/=0) exit

     buffer = adjustl(buffer) !trim white-space left and right
     buffer = trim(buffer)

    if(index(buffer,'#')==0)then
        i = index(buffer,'=')
        if (i==0) then
                i=len_trim(buffer)+1
        else
                val = buffer(i+1:len_trim(buffer))
                cmd = trim(cmd)
        end if
!         print*,'val',val,cmd
        cmd = buffer(1:i-1)
        cmd = trim(cmd)


        i=1
        do while (i<=nnames)
                !print*,names(i)
           if (trim(names(i)) == cmd) then
              exit
           endif
           i=i+1
        enddo


        select case(i)
        case(1)
                read(val,*)n
      !          print*,val
        case(2)
                read(val,*)nup
      !          print*,val
        case(3)
                read(val,*)ndn
      !          print*,val
        case(4)
                read(val,*)l
       !         print*,val
        case(5)
                read(val,*)conv
        !        print*,val
        case(6)
                call read_bond(val)
        !        print*,val
        case(7)
                call onsite(val)
        !        print*,val
        case(8) 
                read(val,*)NA
        end select
       end if

    end do


    close(IUNIT)

end subroutine read_config



subroutine read_bond(fname)

! read hopping terms as:
! i j  t_ij

  character(len=80), intent(in) :: fname
  character(len=80) :: filename

  integer :: i,j,k,ios
  double precision :: t
  logical :: iexist


  filename=adjustl(trim(fname))
  
   print*,'bond file',fname
  ! print*,'No of sites',N

  inquire(file=filename,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'Bond file',filename,' not found'
     stop
  endif
  open(unit=20, file=filename, action="read")

  read(20,*)NBND
        allocate (ARR(N,N))
        allocate (ilbond(2,200))
        
        ARR=0.0d0

        do i = 1, NBND
                read(20,*)j,k,t

                if (j>N.or.j<1.or.k>N.or.k<1) then
                write(*,*) 'Illegal site number on line ',i,' of bond file'
                stop
                endif

               ilbond(1,i) = j
               ilbond(2,i) = k 
                arr(j,k) = -t
                arr(k,j) = -t
        end do
     print*,'Completed assigning hopping terms.'
        close(20)

end subroutine read_bond




subroutine onsite(fname)

! reads the value of Hubbard U

character(len=80),intent(in) :: fname
character(len=80) :: filename

   integer :: i,j
   logical :: iexist

  allocate (U(N))
  U=0.0d0

  filename=adjustl(trim(fname))
  inquire(file=filename,exist=iexist)
  if (iexist.eqv. .false.) then
     write(*,*) 'Site file ',filename,' not found'
     stop
  endif

  open(unit=20, file=filename, action="read")

  do i=1,N
     read(20,*) U(i)
     if (U(i)<0.0d0) then
        write(*,*) 'Error: this version does not handle negative U!'
        stop
     endif
  enddo
  close(20)
end subroutine onsite 



subroutine allocation

 allocate(F(1:N,1:N))
 allocate(F0(1:N,1:N))
 allocate(TMP(1:N))
 allocate(WORK(1:3*N-1))
 allocate(WORK1(1:3*L-1))
 allocate(KPROJECT(N,N))
 allocate(HS(N,0:1))

 allocate(BAS(N,Ne,L))
 allocate(BASUP(N,nup,L))
 allocate(BASDN(N,nup,L))
 allocate(NBASUP(N,ndn,L))
 allocate(NBASDN(N,ndn,L))

 allocate(GREENUP(N,N,L,L))
 allocate(GREENDN(N,N,L,L))
 allocate(GREEN2UP(N,N,L,L))
 allocate(GREEN2DN(N,N,L,L))
 allocate(O(L,L))
 allocate(HAMILT(L,L)) 
 allocate(EVEC(L,L))

end subroutine allocation

end module setupmod
