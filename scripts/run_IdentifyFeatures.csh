#!/bin/csh
#
#run_IdentifyFeatures.csh
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
#######################################################

#Set some general vars
set start_time 	= `date +%s` 	# Record the start time for debugging/optimizing hours request
set DEBUG 	= "debug" 	# Debugging mode? [Outputs netcdf file & binary gridded features if DEBUG="debug", netcdf only if DEBUG="debugonly"]	
set SMOOTH 	= "smth9" 	# Use a 9-point smoother in the fortran code to smooth out mesoscale features smaller than those of interest
set INCU 	= "0"		# Probably not needed anymore. Increment the unit number used by the below-executed fortran code to write output files

#Load the grib modules (needed to use wgrib2 commands)
module load wgrib2

#Set some model info (for future adaptability)
set PARENT 	= /glade/campaign/collections/rda/data/ds084.1 	# Where the GFS RDA data lives
set MODEL 	= "gfs"						# Which model (gfs)
set RES 	= "0p25"					# Model res (0p25 deg)
set EXT 	= "grib2"					# Model data extension (grib2)

#Set some domain info (for -small_grib)
set SLAT 	= "-90"		# Southern domain boundary
set NLAT 	= "90"		# Northern domain boundary
set WLON 	= "0.0"		# Western domain boundary
set ELON 	= "359.75"	# Eastern domain boundary
set CYCLIC 	= "yes"		# Is the domain cyclic? "yes" or "no"

#Set some output file info
set EXT2 	= "nc"		# File extenstion of the data subset submitted to fortran code
set SUB 	= "sub"		# Note that it is a "sub"set
set DIR4MISSING = $SCRATCH/cutofflow/data	# This is eventually the same as the output directory, but need to populate "missing data files" 


#Set datetime lists and fhour
set YYYY 	= ${1}		# Year (from the driver)
set MM 		= ${2}		# Month (from the driver)
set FHOUR 	= ${3}		# Fhour (From the driver)

#Options for month length
set D31 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
set D30 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
set D29 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")
set D28 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28")
set HOURS 	= ("00" "06" "12" "18")
set mm 		= "00"		# "minutes" string

#Get the appropriate month length for given year/month
if( $YYYY == "2012" || $YYYY == "2016" || $YYYY == "2020" )then
  set DL = ($D29)
else
  set DL = ($D28)
endif
if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
  set DX = ($D30)
endif
if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
  set DX = ($D31)
endif
if( $MM == "02" )then
  set DX = ($DL)
endif
  
foreach DD ($DX) 
  foreach hh ($HOURS) 
    set DDIR 	= $PARENT/$YYYY/$YYYY$MM$DD						# Locate the Model data for a the correct date
    set DUM 	= $TMPDIR/$YYYY/$MM/$FHOUR/temp$YYYY$MM$DD$hh$FHOUR.grb    	# Set dummy filename [this MUST be unique if running concurrent wrapper scripts]
    
    #Set input and output names
    set INFILE 	= $DDIR/$MODEL.$RES.$YYYY$MM$DD$hh.$FHOUR.$EXT			# The original model data (formatted initialtime.fhour)	
    set OUTFILE = $TMPDIR/$MODEL.$RES.$YYYY$MM$DD$hh.$FHOUR.$SUB.$EXT2		# The output file of the wgrib2 steps (formatted initialtime.fhour.sub)

    #Set wgrib2 extraction parameters
    set VAR 	= "HGT"
    set TVAR 	= "TMP"
    set UVAR 	= "UGRD"
    set VVAR 	= "VGRD"
    set RVAR 	= "RH"
    set LEV 	= "500 mb"	# This var (and the below string) denote the level at which cutoff lows and troughs are identified at (could use 200 hPa, if so chose)
    set VARSTR 	= "HGT_500mb"
    set LEV0 	= "850 mb"	# These levels are for extracting additional data about the identified cutoff or trough
    set LEV2 	= "200 mb"

    # Long string of levels and variables to extract from the grib2 file 
    set GREPSTR = ':'$VAR':'"$LEV0"':\|:'$VAR':'"$LEV"':\|:'$VAR':'"$LEV2"':\|:'$TVAR':'"$LEV0"':\|:'$TVAR':'"$LEV"':\|:'$TVAR':'"$LEV2"':\|:'$UVAR':'"$LEV0"':\|:'$UVAR':'"$LEV"':\|:'$UVAR':'"$LEV2"':\|:'$VVAR':'"$LEV0"':\|:'$VVAR':'"$LEV"':\|:'$VVAR':'"$LEV2"':\|:'$RVAR':'"$LEV0"':\|:'$RVAR':'"$LEV"':\|:'$RVAR':'"$LEV2"':'

    #Create the dummy file in case something is wrong with the raw model data. This file MUST BE UPDATED if the output format (order, missing values, new vars, etc) of the identification scheme is changed
    #THIS MUST BE UPDATED IF THE OUTPUT FROM THE TRACKING ALGORITHM IS MODIFIED!!!!
    set MISSINGFILE = $DIR4MISSING/$VARSTR/$MODEL.$RES.$YYYY$MM$DD$hh.$FHOUR.dat
    mkdir -p `dirname $MISSINGFILE`
    echo "ITIME,FHOUR,So(m/100km),LAT(N),LON(E),Ro(km),SR,BGo(m/100km),BGo-lat(m/100km),BGo-lon(m/100km),ZMIN(m),ZLAT(N),ZLON(E),CutClosedTrof,Z850(m),Z500(m),Z200(m),T850(K),T500(K),T200(K),U850(m/s),U500(m/s),U200(m/s),V850(m/s),V500(m/s),V200(m/s),MR850(g/kg),MR500(g/kg),MR200(g/kg),600kmZ500(m),600kmT500(K),600kmU500(m/s),600kmV500(m/s),600kmMR500(g/kg)," > $MISSINGFILE
    echo "-9999,-9999,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9" >> $MISSINGFILE

    if ( ! -f $INFILE ) then
      echo $INFILE was missing, do nothing
    else 
      # Set Kasuga et al (2021) parameters
      set MAXR 	= "2100000" 	# maximum radius (meters; 2,100,000 meters = 2,100 kilometers)
      set MINR 	= "100000"	# minimum radius (meters; 100,000 meters = 100 kilometers)
      set NUMR 	= "21"		# how many radii to check (for incrememnts of 100 km) 
      set SRMAX = "2.25" 	# maximum slope ratio (used to restrict weak features in strong background flow)
      set SOMIN = "10.0"
      
      # Subset grib2 file
      wgrib2 $INFILE -s | grep "$GREPSTR"  | wgrib2 -i $INFILE -small_grib $WLON":"$ELON $SLAT":"$NLAT $DUM	# Do the subsetting
      wgrib2 $DUM -netcdf $OUTFILE										# convert the file to netcdf
      rm $DUM
  
      ./identification_algorithm_globe $OUTFILE $MAXR $MINR $NUMR $SRMAX $VARSTR $CYCLIC $DEBUG $INCU $SMOOTH $SOMIN

      rm $OUTFILE  
    endif

  end # hours
end # days
#rm $TMPDIR"/"$YYYY"/"$MM"/"$FHOUR"/"*	# clean the temporaray scratch directory

set end_time = `date +%s`
echo "total time for "$YYYY" "$MM" "$FHOUR" was "`expr $end_time - $start_time`" seconds"
 

