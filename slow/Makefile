version = 05-21-2014


FC = ifort
FCFLAGS = -O1 -ipo
LIB = -L/opt/intel/ictce/3.1/intel64/mkl/10.0.011/lib/em64t -mkl



PROGRAM = pirg

all = $(PROGRAM)

all: variable.o inv.o lapack.o utils.o energy.o utype.o pirg10.o 
	$(FC) $(FCFLAGS) -o pirg lapack.o utils.o inv.o energy.o variable.o  utype.o pirg10.o ${LIB}

variable.o: variable.f90
	$(FC) $(FCFLAGS) -c variable.f90 $(LIB)
inv.o: inv.f90
	$(FC) $(FCFLAGS) -c inv.f90 $(LIB)
lapack.o: lapack.f90
	$(FC) $(FCFLAGS) -c lapack.f90 $(LIB)
utils.o: utils.f90
	$(FC) $(FCFLAGS) -c utils.f90 $(LIB)
energy.o: energy.f90
	$(FC) $(FCFLAGS) -c energy.f90 $(LIB)
utype.o: utype.f90
	$(FC) $(FCFLAGS) -c utype.f90  $(LIB)
pirg10.o: pirg10.f90
	$(FC) $(FCFLAGS) -c pirg10.f90  $(LIB)
#pirg-rtc.o: pirg-rtc.f90
#
#$(FC) $(FCFLAGS) -c pirg-rtc.f90 $(LIB)

#variable.f90 inv.f90 lapack.f90 utils.f90 energy.f90 utype.f90 pirg10.f90




clean: 
	rm -f pirg setup.o variable.o inv.o lapack.o utils.o energy.o utype.o pirg10.o *.o, *.mod,*~
                                                                                                                                                                                                                                                                                                                                      
