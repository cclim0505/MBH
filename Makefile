### gfortran compiler options
#FC = gfortran
#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

### mpifortran compilre options
FC = mpifort
FFLAGS = -llapack 

DRIVER = main.f

SIMULATION = simulation.f
MPIVAR = mpi_var.f

CONSTANTS = constants.f
CGE = coord_grad_ene.f
GUPTA = gupta.f
DFTB = ene_dftb.f
RANDOM = random_coord.f
POTENTIAL = potential.f
BASIN = basin_hopping.f
MONTE = monte.f
INERTIA = inertia.f
SPLICE = cut_splice.f

PARAM = param.f		#module to compliment senior's code
GRAD = grad.f		#gradient subroutine developed by senior
ARRMAT = array_matrix.f
#ENERGIES = energies.f

INITIALISE = initialise.f

BLAS = blas.f
LINPACK = linpack.f
TIMER = timer.f
LBFGS = lbfgsb.f
OPTIM = optimization.f


## Building drivers
DRIVER_OPEN = 


all :  main


main : $(CONSTANTS) $(MPIVAR) $(CGE) $(GUPTA) $(DFTB) $(RANDOM) $(POTENTIAL) $(BASIN) $(MONTE) $(INERTIA) $(SPLICE) $(INITIALISE) $(ARRMAT) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(SIMULATION) $(DRIVER) 
	$(FC) $(FFLAGS) $(CONSTANTS) $(MPIVAR) $(CGE) $(GUPTA) $(DFTB) $(RANDOM) $(POTENTIAL) $(BASIN) $(MONTE) $(INERTIA) $(SPLICE) $(INITIALISE) $(ARRMAT) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(SIMULATION) $(DRIVER)  -o run.out



clean :
	rm run.out
	rm *.o
	rm *.mod

