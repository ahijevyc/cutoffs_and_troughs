#!/bin/csh

#Kevin Lupo
#klupo@ucar.edu
#
#checktotals.csh
#Make sure all .dat files were created (includes placeholders). 30 day months should have 120 per fhour, 31 day months should have 124, leap febs should have 116, non-leap febs should have 112
#########################################


set PARENT = "~/work_new/postdoc/kasugaEA21/version9/HGT_500mb/"
set MODEL = "gfs"
set RES = "0p25"
set BASE = "$PARENT$MODEL.$RES."
set EXT = ".dat"
set YEARS = ("2015" "2017") #"2016" "2017" "2018" "2019" "2021" "2022")
#set FHOURS = ("f000" "f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
set FHOURS = ("f066" "f072" "f126" "f132" "f138" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
foreach YYYY ($YEARS)

if( $YYYY == "2022" )then
  set MONTHS = ("01" "02" "03")
else
  set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
endif

foreach MM ($MONTHS)

if( $YYYY == "2012" || $YYYY == "2016" || $YYYY == "2020" )then
  set FTOT = "116"
else
  set FTOT = "112"
endif
if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
  set TOT = "120"
endif
if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
  set TOT = "124"
endif
if( $MM == "02" )then
  set TOT = $FTOT
endif

foreach FHOUR ($FHOURS)
  set CHECK = `ls $BASE$YYYY$MM*"."$FHOUR$EXT | wc -l`
  if( $CHECK != $TOT )then
    echo "$YYYY $MM $FHOUR -- $CHECK / $TOT - ERROR"
  else
    echo "$YYYY $MM $FHOUR -- $CHECK / $TOT"
  endif  
end 
end
end 

