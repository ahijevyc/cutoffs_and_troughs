#!/bin/csh

#Kevin Lupo
#klupo@ucar.edu
#
#driver_kasugaEA21.csh batch submits the queue_cutoffID script for given itimes & fhours
#8 Dec 2021
#keep # of jobs sent to queue relatively low, submit in groups of 4
#
#9 Dec 2021
#can go larger...htc limits are 468 CPUs/4680 GB memory per user at any one time
#...so, could technically submit all fhours simulataneously...
#...lets do 168 and 180 jobs (12 months x 14 and 15 hours)
#########################################



#set YEARS = ("2015" "2016" "2017" "2018" "2019" "2021" "2022")
#set YEARS = ("2017" "2018" "2019")
set YEARS = ("2022")
#set MONTHS = ("01")# ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")


#set FHOURS = ("f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
#set FHOURS = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168")
#set FHOURS = ("f000") 
#set FHOURS = ("f066" "f072" "f126" "f132" "f138")
set FHOURS = ("f000" "f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
set dispatch = "dispatch_identification.csh"

foreach FHOUR ($FHOURS)
foreach YEAR ($YEARS)

if( $YEAR == "2022" )then
  set MONTHS = ("04" "05" "06" "07" "08")
#else if ($YEAR == "2015" ) then
#  set MONTHS = ("03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
else
  set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
endif

foreach MM ($MONTHS)
  cp queue_identification_template queue_identification
  echo "./"$dispatch" "$YEAR" "$MM" "$FHOUR >> queue_identification
  qsubcasper queue_identification
  rm queue_identification
end 
end
end 

