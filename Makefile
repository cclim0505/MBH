FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

DRIVER = main.f

CONSTANTS = constants.f
TYPES = types.f		#not initialised
GUPTA = gupta.f
INITIALISE = initialise.f
PARAM = param.f
GRAD = grad.f		#gradient subroutine developed by senior
RANDOM = random_coord.f
ARRMAT = array_matrix.f
BASIN = basin_hopping.f

BLAS = blas.f
LINPACK = linpack.f
TIMER = timer.f
LBFGS = lbfgsb.f
OPTIM = optimization.f


## Building drivers
DRIVER_OPEN = 


all :  main





main : $(CONSTANTS) $(TYPES) $(GUPTA) $(INITIALISE) $(RANDOM) $(ARRMAT) $(BASIN) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(DRIVER) 
	$(FC) $(FFLAGS) $(CONSTANTS) $(TYPES) $(GUPTA) $(INITIALISE) $(RANDOM) $(ARRMAT) $(BASIN) $(PARAM) $(GRAD) $(BLAS) $(LINPACK) $(TIMER) $(LBFGS) $(OPTIM) $(DRIVER)  -o run.out



clean :
	rm run.out
	rm *.o
	rm *.mod

