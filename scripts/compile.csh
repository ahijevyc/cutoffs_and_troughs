module reset
rm *.o *.mod
ifort -ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl -I${NETCDF}/include -c netcdf_routines_mod.f90
ifort -ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o identification_algorithm_global_noDisambigSteps.f90 -o identification_algorithm_globe
ifort -ip -fp-model precise -w -ftz -align all -fno-alias -FR -convert big_endian -assume byterecl -lnetcdff -I${NETCDF}/include netcdf_routines_mod.o identification_algorithm_global_noDisambigSteps_cases.f90 -o identification_algorithm_globe_cases
