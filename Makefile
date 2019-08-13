FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

DRIVER = main.f

CONSTANTS = constants.f
TYPES = types.f
GUPTA = gupta.f
INITIALISE = initialise.f
PARAM = param.f
GRAD = grad.f
RANDOM = random_coord.f
BASIN = basin_hopping.f

## Building drivers
DRIVER_OPEN = 


all :  main





main : $(CONSTANTS) $(TYPES) $(GUPTA) $(INITIALISE) $(RANDOM) $(BASIN) $(PARAM) $(GRAD) $(DRIVER) 
	$(FC) $(FFLAGS) $(CONSTANTS) $(TYPES) $(GUPTA) $(INITIALISE) $(RANDOM) $(BASIN) $(PARAM) $(GRAD) $(DRIVER)  -o run.out



clean :
	rm run.out
	rm *.o

