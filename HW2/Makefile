FC=gfortran
FFLAGS=-O3 -Wall -Wextra
SRC1=seed.f90 HW_2_halfsphere.f90
OBJ=${SRC1:.f90=.o}

%.o: %.f90	
	$(FC) $(FFLAGS) -o $@ -c $<
halfsphere: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	rm *.o *.mod 