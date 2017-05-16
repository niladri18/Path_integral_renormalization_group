program noninteracting


     use setupmod
     use pirg
     use variable
     use lapack
     use inv
     use util
     use energy

 tau=0.01
 OUNIT=22
 BUNIT=23 ! for the binary data
 EVUNIT=45 ! for the eigenvector

  call get_command_argument(1,FILENAME)

  open(unit=OUNIT, file=trim(FILENAME)//'.out', status='replace', &
          action='write')

  call cpu_time(start)

  call read_config 
  print*,'files read'

  open(unit=BUNIT, file=trim(FILENAME)//'.bin', form ='UNFORMATTED', &
                        status='replace', action = 'write')
  open(unit=EVUNIT, file=trim(FILENAME)//'.wgt',status='replace', &
                        action = 'write')
 

   
        call pirg_run
  call cpu_time(finish)
  write(OUNIT, '(A,F12.5)')'Time elapsed = ',finish-start




end program noninteracting
 
