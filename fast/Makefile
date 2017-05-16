version = 05-21-2014


FC = ifort
FCFLAGS = -O1 -ipo
LIB = -L/opt/intel/ictce/3.1/intel64/mkl/10.0.011/lib/em64t -mkl



PROGRAM = pirg

all = $(PROGRAM)

all: variable1.o setup.o inv.o lapack.o utils.o energy-rtc.o pirg-fast.o main.o 
	$(FC) $(FCFLAGS) -o pirg setup.o lapack.o utils.o inv.o energy-rtc.o variable1.o  pirg-fast.o main.o ${LIB}

variable1.o: variable1.f90
	$(FC) $(FCFLAGS) -c variable1.f90 $(LIB)
setup.o: setup.f90
	$(FC) $(FCFLAGS) -c setup.f90
inv.o: inv.f90
	$(FC) $(FCFLAGS) -c inv.f90 $(LIB)
lapack.o: lapack.f90
	$(FC) $(FCFLAGS) -c lapack.f90 $(LIB)
utils.o: utils.f90
	$(FC) $(FCFLAGS) -c utils.f90 $(LIB)
energy-rtc.o: energy-rtc.f90
	$(FC) $(FCFLAGS) -c energy-rtc.f90 $(LIB)
pirg-fast.o: pirg-fast.f90
	$(FC) $(FCFLAGS) -c pirg-fast.f90  $(LIB)
main.o: main.f90
	$(FC) $(FCFLAGS) -c main.f90  $(LIB)
#pirg-rtc.o: pirg-rtc.f90
#
#$(FC) $(FCFLAGS) -c pirg-rtc.f90 $(LIB)






clean: 
	rm -f pirg setup.o variable.o inv.o lapack.o utils.o energy.o utype.o pirg10.o *.o, *.mod,*~
                                                                                                                                                                                                                                                                                                                                      
