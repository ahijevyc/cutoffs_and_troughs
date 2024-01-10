module reset
rm *.o *.mod
set FCFLAGS="-ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl "
ifort $FCFLAGS -I${NETCDF}/include -c netcdf_routines_mod.f90
ifort $FCFLAGS -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o identification_algorithm_global_noDisambigSteps.f90 -o identification_algorithm_globe
ifort $FCFLAGS -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o identification_algorithm_global_noDisambigSteps_cases.f90 -o identification_algorithm_globe_cases
ifort -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o track_analysis.f90 -o track_analysis
unset FCFLAGS
