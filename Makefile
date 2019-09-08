FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

DRIVER = main.f

CONSTANTS = constants.f
CGE = coord_grad_ene.f
GUPTA = gupta.f
DFTB = ene_dftb.f
RANDOM = random_coord.f
POTENTIAL = potential.f
BASIN = basin_hopping.f
MONTE = monte.f
INERTIA = inertia.f

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


main : $(CONSTANTS) $(CGE) $(GUPTA) $(DFTB) $(RANDOM) $(POTENTIAL) $(BASIN) $(MONTE) $(INERTIA) $(INITIALISE) $(ARRMAT) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(DRIVER) 
	$(FC) $(FFLAGS) $(CONSTANTS) $(CGE) $(GUPTA) $(DFTB) $(RANDOM) $(POTENTIAL) $(BASIN) $(MONTE) $(INERTIA) $(INITIALISE) $(ARRMAT) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(DRIVER)  -o run.out



clean :
	rm run.out
	rm *.o
	rm *.mod

