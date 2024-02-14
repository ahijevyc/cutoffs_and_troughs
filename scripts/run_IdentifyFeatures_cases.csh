#!/bin/csh
#
#run_IdentifyFeatures_cases.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Wrapper file to....
# 	---extract 500hPa height from gfs grib2 files
#	---call fortran code to identify cutoff lows using kasugaEA21 algorithm 
#29 Sep 2021
###########
#8 Nov 2021
# 	--Code is submitted via driver->dispatch to casper & submits each month as a separate job
#15 Dec 2021
#	--Improved batch submission format (driver_IdentifyFeatures.csh --> run_IdentifyFeatures.csh --> ./identification_algorithm_globe
###########
#15 Sep 2023
#	--This version of the code is used for UFS case studies. Extract desired fields from .nc stitched output
#######################################################

#Set some general vars
set start_time 	= `date +%s` 	# Record the start time for debugging/optimizing hours request
set DEBUG 	= "debug" 	# Debugging mode? [Outputs netcdf file & binary gridded features if DEBUG="debug", netcdf only if DEBUG="debugonly"]	
set SMOOTH 	= "smth9" 	# Use a 9-point smoother in the fortran code to smooth out mesoscale features smaller than those of interest
set INCU 	= "0"		# Probably not needed anymore. Increment the unit number used by the below-executed fortran code to write output files

#Load the grib modules (needed to use wgrib2 commands)
source /etc/profile.d/z00_modules.csh

module load wgrib2
module load nco
#Set datetime lists and fhour
set ITIME 	= ${1}		# ITIME (YYYYMMDDhh, from the driver)
set FHOUR 	= ${2}		# FHOUR (fhhh, from the driver)

#Set some model info (for future adaptability)
set PARENT 	= /glade/scratch/klupo/UFS-MRW/UFS_PRODRUNS_OUTPUT 	# Where the UFS data lives 
set PARENT 	= /glade/campaign/mmm/parc/mwong/ufs-mrw 	# Where the UFS data lives 
set FLEN 	= "F240"						# Forecast length (F240)
set RES 	= "C768"				# Model res (C768, output is on 0.25latlon grid)
set EXT 	= "nc"					# Model data extension (nc)

#Set some domain info (for -small_grib)
set SLAT 	= "-90"		# Southern domain boundary
set NLAT 	= "90"		# Northern domain boundary
set WLON 	= "0.0"		# Western domain boundary
set ELON 	= "359.75"	# Eastern domain boundary
set CYCLIC 	= "yes"		# Is the domain cyclic? "yes" or "no"

#Set some output file info

  
#Data management tasks
set DDIR 	= $PARENT/$ITIME.$FLEN.$RES					# Locate the Model data for a the correct date

#Set input and output names
set INFILE  = $DDIR/interp/interp_fv3_history2d_${ITIME}_$FHOUR.$EXT		# The original model data file (formatted initialtime.fhour)    
set ODIR=.	# outdir for identification_algorithm_global_noDisambigSteps_cases
#set OUTFILE = $TMPDIR/$ITIME"."$FHOUR"."$EXT	    	# The the subset file
#ncks -O -v z500,t500,u500,v500,rh500 $INFILE $OUTFILE

# Set Kasuga et al (2021) parameters
set MAXR  = "2100000"	  # maximum radius (meters; 2,100,000 meters = 2,100 kilometers)
set MINR  = "100000"	  # minimum radius (meters; 100,000 meters = 100 kilometers)
set NUMR  = "21"	  # how many radii to check (for incrememnts of 100 km) 
set SRMAX = "2.25"	  # maximum slope ratio (used to restrict weak features in strong background flow)
set SOMIN = "10.0"
      

./identification_algorithm_globe_cases $INFILE $ITIME $FHOUR $MAXR $MINR $NUMR $SRMAX $CYCLIC $DEBUG $ODIR $SMOOTH $SOMIN

#rm $OUTFILE  

set end_time = `date +%s`
echo "total time for "$ITIME $FHOUR" was "`expr $end_time - $start_time`" seconds"
