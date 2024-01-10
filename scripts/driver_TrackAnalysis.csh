#!/bin/csh
#
#driver_TrackAnalysis.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Wrapper file to....
# 	track features identified by KS21 algorithm at the analysis time only
#21 Oct 2021
#16 Mar 2022
#	Expanded dataset to include 2021-Feb 2022 [15 Jan 2015 -- 28 Feb 2022 ]
#3 May 2022
#	Removed hemisphere extensions. While this was convenient, features are now permitted over globe, and are allowed to cross the equator (if any exist)
###########
#

#set HEXT = ${1}
#if ( $HEXT == "NH" ) then
#  set HEMIS = "NHEM"
#endif
#if ( $HEXT == "SH" ) then
#  set HEMIS = "SHEM"
#endif

set VARSTR = "HGT_500mb"
#set PARENT = /glade/work/klupo/postdoc/kasugaEA21/version9/$VARSTR
set MODEL = "gfs"
set RES = "0p25"
set EXT = "sub.dat"
set DEFAULT_PMAX = "1.5"
set DEFNORM_So = "15.0"
set DEFNORM_BGo = "40.0"
set DEFNORM_Ro = "1000.0"
set DEFPREDERR = "1500.0"
set DEFOPPX    = "1000.0"
set DEFOPPY    = "1000.0"
set DEFMAXDIST = "1500.0"

set USECONFIG = "pmax1.5_2stdnorms_munozDmax1200_oppmax700" #${1}

if( $USECONFIG == "default" ) then
  set PMAX = $DEFAULT_PMAX
  set NORM_So = $DEFNORM_So
  set NORM_BGo = $DEFNORM_BGo
  set NORM_Ro = $DEFNORM_Ro
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.0" ) then
  set PMAX = "1.0"
  set NORM_So = $DEFNORM_So
  set NORM_BGo = $DEFNORM_BGo
  set NORM_Ro = $DEFNORM_Ro
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.0_meannorms" ) then
  set PMAX = "1.0"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax2.0_meannorms" ) then
  set PMAX = "2.0"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.0_stdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "6.8"
  set NORM_BGo = "8.1"
  set NORM_Ro = "316.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.0_2stdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "13.6"
  set NORM_BGo = "16.2"
  set NORM_Ro = "632.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.0_meanPLSstdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "24.4"
  set NORM_BGo = "22.4"
  set NORM_Ro = "904.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $USECONFIG == "pmax1.5_meannorms_munozDmax1200_oppmax700" ) then
  set PMAX = "1.5"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = "1200.0"
  set OXMAX = "700.0"
  set OYMAX = "700.0"
  set DMAX = "1200.0"
endif
if( $USECONFIG == "pmax1.5_2stdnorms_munozDmax1200_oppmax700" ) then
  set PMAX = "1.5"
  set NORM_So = "13.6"
  set NORM_BGo = "16.2"
  set NORM_Ro = "632.0"
  set EMAX = "1200.0"
  set OXMAX = "700.0"
  set OYMAX = "700.0"
  set DMAX = "1200.0"
endif
#Set forecast hour
set FHOUR = "f000"

#Set file list for fortran code
#set LNAME = $PARENT/AllInits.$FHOUR.list
#set LNAME = "2020Inits.f000.list"
#echo `ls $PARENT$MODEL"."$RES"."*"."$FHOUR"."$EXT"."$HEXT > $LNAME`
set LNAME=$TMPDIR/$$.txt

ls $TMPDIR/????/??/f???/$MODEL.$RES.*.$FHOUR.$EXT > $LNAME

./track_analysis $LNAME $VARSTR $PMAX $NORM_So $NORM_BGo $NORM_Ro $EMAX $OXMAX $OYMAX $DMAX


