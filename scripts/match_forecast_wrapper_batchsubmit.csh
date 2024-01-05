#!/bin/csh
#
#match_forecast_wrapper.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Wrapper file to....
# 	match features identified by KS21 algorithm at fhours to tracked analysis features
#15 Feb 2022
#	
###########
#

set f = ${1}
set FHOUR = ${2}
#set HEXT = ${3}

#if ( $HEXT == "NH" ) then
#  set HEMIS = "NHEM"
#endif
#if ( $HEXT == "SH" ) then
#  set HEMIS = "SHEM"
#endif

set VARSTR = "HGT_500mb"
set PARENT = "/glade/work/klupo/postdoc/kasugaEA21/version8/"$VARSTR"/"
set MODEL = "gfs"
set RES = "0p25"
set EXT	= "track"
set EXT2 = "dat"

#Set datetime lists
#KL - don't actually need to set these, since the wildcard character in the echo > LNAME command below will include all years/months in the list file in old->new order

set DDIR = $PARENT

#Set verification  file list for fortran code
set VNAME = $PARENT"f000.trackfiles.list"
#echo `ls $PARENT$MODEL"."$RES"."*".f000."$EXT.$HEXT > $VNAME`

#Refresh the top level list file
set MULTINAME = $PARENT"fxlists"$f".list"
if ( -f "$MULTINAME" ) then
  rm $MULTINAME
endif
touch $MULTINAME

#Set fhour dat file list for fortran code
set LNAME = $PARENT"."$FHOUR".list"
echo `ls $PARENT$MODEL"."$RES"."*"."$FHOUR"."$EXT2 > $LNAME`
echo $LNAME >> $MULTINAME


./match_forecast $VNAME $MULTINAME $VARSTR $f
