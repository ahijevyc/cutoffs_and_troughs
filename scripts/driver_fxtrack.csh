#!/bin/csh

#Kevin Lupo
#klupo@ucar.edu
#
#driver_fxtrack.csh batch submits the queue_fxtrack script for given itimes & fhours
#27 Apr 2022
#Submit in groups of 2 years; performance issues if too many submited at once
#########################################

#set HEXT = ${1}
set TCONFIGS = ("default" "pmax1.0" "pmax1.0_meannorms" "pmax2.0_meannorms" "pmax1.0_stdnorms" "pmax1.0_2stdnorms" "pmax1.0_meanPLSstdnorms" "pmax1.5_meannorms_munozDmax1200_oppmax700" "pmax1.5_2stdnorms_munozDmax1200_oppmax700")
set MCONFIGS = ("default" "pmax2.0_fixedradius_2stdnorms_equalweights_munozEmax800")
set SELECTTCONFIG = 9 #${1}
set SELECTMCONFIG = 2
set TCONFIG = $TCONFIGS[$SELECTTCONFIG]
set MCONFIG = $MCONFIGS[$SELECTMCONFIG]
echo "Tracking using the $TCONFIG configuration"
echo "Matching using the $MCONFIG configuration"
set YEARS = ("2022") #"2017" "2018" "2019" "2020" "2021")
set MONTHS =  ("01" "02" "03" "04" "05" "06" "07" "08")# "09" "10" "11" "12")

set dispatch = "dispatch_fxtrack.csh"

foreach YEAR ($YEARS)
foreach MM ($MONTHS)
  cp queue_fxtrack_template queue_fxtrack
  echo "./"$dispatch" "$YEAR" "$MM" "$TCONFIG" "$MCONFIG >> queue_fxtrack
  qsubcasper queue_fxtrack
  rm queue_fxtrack
end 
end

