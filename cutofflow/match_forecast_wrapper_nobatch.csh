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

set HEXT = ${1}

if ( $HEXT == "NH" ) then
  set HEMIS = "NHEM"
endif
if ( $HEXT == "SH" ) then
  set HEMIS = "SHEM"
endif

set VARSTR = "HGT_500mb"
set PARENT = "/glade/work/klupo/postdoc/kasugaEA21/version7/"$VARSTR"/"$HEMIS"/"
set MODEL = "gfs"
set RES = "0p25"
set EXT	= "track"
set EXT2 = "dat"

#Set datetime lists
#KL - don't actually need to set these, since the wildcard character in the echo > LNAME command below will include all years/months in the list file in old->new order
#set YEARS = ("2015" "2016" "2017" "2018" "2019" "2020")
#set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

set DDIR = $PARENT

#Set forecast hour
set FHOUR = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")

#Set verification  file list for fortran code
set VNAME = $PARENT".f000.trackfiles.list"
echo `ls $PARENT$MODEL"."$RES"."*".f000."$EXT.$HEXT > $VNAME`

#Refresh the top level list file
set MULTINAME = $PARENT".fxlists.list"
if ( -f "$MULTINAME" ) then
  rm $MULTINAME
endif
touch $MULTINAME

#Set fhour dat file list for fortran code
foreach FF ($FHOUR)
  set LNAME = $PARENT"."$FF".list"
  echo `ls $PARENT$MODEL"."$RES"."*"."$FF"."$EXT2.$HEXT > $LNAME`
  echo $LNAME >> $MULTINAME
end



./match_forecast $VNAME $MULTINAME $VARSTR $HEMIS $HEXT

#rm $OUTFILE
#end 
#end
#end 
#end 
#end
