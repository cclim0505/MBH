PROGRAM = run.out
### gfortran compiler options
#FC = gfortran
#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

VPATH = src

### mpifortran compilre options
FC = mpifort
FFLAGS = -llapack 


### source names
DRIVER = main
SIMUL = simulation

MPIVAR = mpi_var
DIR = directory
CONST = constants
CGE = coord_grad_ene
GUPTA = gupta
DFTB = ene_dftb
RANDOM = random_coord
POTENT = potential
BASIN = basin_hopping
MONTE = monte
INERTIA = inertia
SPLICE = cut_splice
GEOMETRIC = geometric_drive
MOVES = moves

PARAM = param			#module to compliment senior's code
GRAD = grad			#gradient subroutine developed by senior
ARRMAT = array_matrix
#ENERGIES = energies

INIT = initialise

BLAS = blas
LINPACK = linpack
TIMER = timer
LBFGS = lbfgsb
OPTIM = optimization


### OBJECT LIST
OBJS = 	$(CONST).o\
	$(DIR).o\
	$(MPIVAR).o\
	$(CGE).o\
	$(GUPTA).o\
	$(DFTB).o\
	$(RANDOM).o\
	$(POTENT).o\
	$(BASIN).o\
	$(MONTE).o\
	$(INERTIA).o\
	$(SPLICE).o\
	$(GEOMETRIC).o\
	$(MOVES).o\
	$(INIT).o\
	$(ARRMAT).o\
	$(BLAS).o\
	$(LINPACK).o\
	$(TIMER).o\
	$(LBFGS).o\
	$(OPTIM).o\
	$(SIMUL).o\
	$(DRIVER).o


all :  main

main: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(PROGRAM)

$(MPIVAR).o: $(MPIVAR).f
	$(FC) -c $(FFLAGS) $< 

$(CONST).o: $(CONST).f
	$(FC) -c $(FFLAGS) $< 

$(DIR).o: $(DIR).f
	$(FC) -c $(FFLAGS) $< 

$(CGE).o: $(CGE).f $(CONST).o
	$(FC) -c $(FFLAGS) $< 

$(GUPTA).o: $(GUPTA).f $(CONST).o $(DIR).o $(CGE).o
	$(FC) -c $(FFLAGS) $< 

$(DFTB).o: $(DFTB).f $(CONST).o $(CGE).o
	$(FC) -c $(FFLAGS) $< 

$(RANDOM).o: $(RANDOM).f $(CONST).o $(DIR).o $(CGE).o $(MPIVAR).o
	$(FC) -c $(FFLAGS) $< 

$(POTENT).o: $(POTENT).f $(CONST).o $(GUPTA).o $(DFTB).o
	$(FC) -c $(FFLAGS) $< 

$(BASIN).o: $(BASIN).f $(CONST).o $(DIR).o $(CGE).o $(RANDOM).o $(GUPTA).o $(POTENT).o
	$(FC) -c $(FFLAGS) $< 

$(MONTE).o: $(MONTE).f $(CONST).o $(DIR).o $(CGE).o
	$(FC) -c $(FFLAGS) $< 

$(INERTIA).o: $(INERTIA).f $(CONST).o $(CGE).o
	$(FC) -c $(FFLAGS) $< 

$(SPLICE).o: $(SPLICE).f $(CONST).o $(CGE).o $(INERTIA).o
	$(FC) -c $(FFLAGS) $< 

$(GEOMETRIC).o: $(GEOMETRIC).f $(CONST).o $(CGE).o $(RANDOM).o $(BASIN).o
	$(FC) -c $(FFLAGS) $< 

$(MOVES).o: $(MOVES).f $(BASIN).o $(SPLICE).o $(MONTE).o $(GEOMETRIC).o
	$(FC) -c $(FFLAGS) $< 

$(ARRMAT).o: $(ARRMAT).f $(CONST).o $(CGE).o
	$(FC) -c $(FFLAGS) $< 

$(INIT).o: $(INIT).f $(DIR).o  $(CGE).o $(POTENT).o $(GUPTA).o $(RANDOM).o $(MONTE).o $(BASIN).o
	$(FC) -c $(FFLAGS) $< 

$(BLAS).o: $(BLAS).f 
	$(FC) -c $(FFLAGS) $< 

$(LINPACK).o: $(LINPACK).f 
	$(FC) -c $(FFLAGS) $< 

$(TIMER).o: $(TIMER).f 
	$(FC) -c $(FFLAGS) $< 

$(LBFGS).o: $(LBFGS).f $(BLAS).o $(LINPACK).o $(TIMER).o
	$(FC) -c $(FFLAGS) $< 

$(OPTIM).o: $(OPTIM).f $(CONST).o $(CGE).o $(POTENT).o $(ARRMAT).o $(DFTB).o
	$(FC) -c $(FFLAGS) $< 

$(SIMUL).o: $(SIMUL).f $(INIT).o $(CGE).o $(CONST).o $(SPLICE).o $(INERTIA).o $(RANDOM).o\
		$(OPTIM).o $(MONTE).o $(MOVES).o $(POTENT).o
	$(FC) -c $(FFLAGS) $< 

$(DRIVER).o: $(DRIVER).f $(SIMUL).o $(MPIVAR).o
	$(FC) -c $(FFLAGS) $< 

clean:
	rm *.o
	rm *.mod
cleaner:
	rm run.out

