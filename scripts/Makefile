FC=ifort

FCFLAGS = -L$(NETCDF)/lib -lnetcdff -I$(NETCDF)/include
OBJ = netcdf_routines_mod.o

# List of Fortran programs
TRACK=track_analysis track_forecast track_forecast_cases
ID=identification_algorithm_globe identification_algorithm_globe_cases

# Default target: Build all programs
all: $(TRACK) $(ID)

# Rule to build each program

$(TRACK): %: %.f90
	$(FC) $(FCFLAGS) $^ -o $@

identification_algorithm_globe: identification_algorithm_global_noDisambigSteps.f90
	$(FC) $(FCFLAGS) $^ -o $@

identification_algorithm_globe_cases: identification_algorithm_global_noDisambigSteps_cases.f90 $(OBJ)
	$(FC) $(FCFLAGS) $^ -o $@

$(OBJ): netcdf_routines_mod.f90
	$(FC) -c $^

clean:
	rm *.o *.mod $(TRACK) $(ID)
