#!/bin/csh
#
#kasugaEA21_wrapper.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Wrapper file to....
# 	---extract 500hPa height from gfs grib2 files
#	---call fortran code to identify cutoff lows using kasugaEA21 algorithm 
#29 Sep 2021
###########
#
#Debugging mode? [Outputs netcdf file & binary gridded features if DEBUG="debug"]
set DEBUG = "binary"	

#Increment the unit number used by the below-executed fortran code to write output files
set INCU = "0"

#Set some model info
set PARENT = "/glade/collections/rda/data/ds084.1/"
set MODEL = "gfs"
set RES = "0p25"
set EXT = "grib2"

#Set some domain info
set SLAT = "-90"
set NLAT = "90"
set WLON = "0.0"
set ELON = "359.75"
set CYCLIC = "yes"	# "yes" or "no"

#Set some output file info
set EXT2 = "nc"
set SUB = "sub"
set DIR4MISSING = "/glade/work/klupo/postdoc/kasugaEA21/version7/"
set TMPDIR = "/glade/scratch/klupo/ks21_temp/"

#Set datetime lists
set YEARS = ("2018")
set MONTHS = ("01")# "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

set D31 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
set D30 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
set D29 =  ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")
set D28 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28")
set D31 = ("01")

set HOURS = ("00" "06" "12" "18")
set mm = "00"

#Set initial time info and directory
foreach YYYY ($YEARS) #= "2020"
  if( $YYYY == "2012" || $YYYY == "2016" || $YYYY == "2020" )then
  set DL = ($D29)
  else
  set DL = ($D28)
  endif
foreach MM ($MONTHS) #= "09"
  if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
  set DX = ($D30)
  endif
  if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
  set DX = ($D31)
  endif
  if( $MM == "02" )then
  set DX = ($DL)
  endif
  
foreach DD ($DX) #DD = "09"
foreach hh ($HOURS) # hh = "00"
set DDIR = $PARENT$YYYY"/"$YYYY$MM$DD"/"

#Set forecast hour
foreach FHOUR ("f000")# "f024" "f048" "f072")

#Set dummy filename [this MUST be unique if running concurrent wrapper scripts]
set DUM = $TMPDIR"temp"$YYYY$MM$DD$hh$FHOUR".grb"

#Set input and output names
set INFILE = $DDIR$MODEL"."$RES"."$YYYY$MM$DD$hh"."$FHOUR"."$EXT
set OUTFILE = $TMPDIR$MODEL"."$RES"."$YYYY$MM$DD$hh"."$FHOUR"."$SUB"."$EXT2

#Set wgrib2 extraction parameters
set VAR = "HGT"
set TVAR = "TMP"
set UVAR = "UGRD"
set VVAR = "VGRD"
set RVAR = "RH"
set LEV = "500 mb"
set VARSTR = "HGT_500mb"


set LEV0 = "850 mb"
set LEV2 = "200 mb"

set GREPSTR = ':'$VAR':'"$LEV0"':\|:'$VAR':'"$LEV"':\|:'$VAR':'"$LEV2"':\|:'$TVAR':'"$LEV0"':\|:'$TVAR':'"$LEV"':\|:'$TVAR':'"$LEV2"':\|:'$UVAR':'"$LEV0"':\|:'$UVAR':'"$LEV"':\|:'$UVAR':'"$LEV2"':\|:'$VVAR':'"$LEV0"':\|:'$VVAR':'"$LEV"':\|:'$VVAR':'"$LEV2"':\|:'$RVAR':'"$LEV0"':\|:'$RVAR':'"$LEV"':\|:'$RVAR':'"$LEV2"':'

#Create a dummy file if the model data is missing.
#THIS MUST BE UPDATED IF THE OUTPUT FROM THE TRACKING ALGORITHM IS MODIFIED!!!!
if ( ! -f "$INFILE" ) then
  #echo "$INFILE does not exist. Generating dummy .dat file"
  #set MISSINGFILE = $DIR4MISSING$VARSTR"/"$MODEL"."$RES"."$YYYY$MM$DD$hh"."$FHOUR".dat"
  #echo "ITIME,FHOUR,So(m/100km),LAT(N),LON(E),Ro(km),SR,BGo(m/100km),BGo-lat(m/100km),BGo-lon(m/100km),ZMIN(m),ZLAT(N),ZLON(E),Z850(m),Z500(m),Z200(m),T850(K),T500(K),T200(K),U850(m/s),U500(m/s),U200(m/s),V850(m/s),V500(m/s),V200(m/s),RH850(%),RH500(%),RH200(%)" > $MISSINGFILE
  #echo "-9999,-9999,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9" >> $MISSINGFILE

else 
  #Set Kasuga et al (2021) parameters
  set MAXR = "2100000" 	# maximum radius
  set MINR = "100000"	# minimum radius
  set NUMR = "21"		# how many radii to check (for incrememnts of 100 km) 
  set SRMAX = "3.0" 	# maximum slope ratio

  #Subset grib2 file
  wgrib2 $INFILE -s | grep "$GREPSTR"  | wgrib2 -i $INFILE -small_grib $WLON":"$ELON $SLAT":"$NLAT $DUM
  wgrib2 $DUM -netcdf $OUTFILE
  rm $DUM
  
  ./identification_algorithm $OUTFILE $MAXR $MINR $NUMR $SRMAX $VARSTR $CYCLIC $DEBUG $INCU
  #./kasugaEA21_algorithm $OUTFILE $MAXR $MINR $NUMR $SRMAX $VARSTR $CYCLIC $DEBUG

  rm $OUTFILE  
endif

end 
end
end 
end 
end
