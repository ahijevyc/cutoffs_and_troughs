rm *.o *.mod
ifort -L$MPI_ROOT/lib -lmpi -ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include -c netcdf_routines_mod.f90
ifort  -L$MPI_ROOT/lib -lmpi -ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o identification_algorithm_global_noDisambigSteps_cases.f90 -o identification_algorithm_globe_cases

