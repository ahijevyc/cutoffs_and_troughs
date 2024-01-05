#!/bin/csh

#Kevin Lupo
#klupo@ucar.edu
#
#
#15 Feb 2022
#driver_matching.csh batch submits the queue_matching script for given fhours
#Because of how the .f90 code is configured to interpret "f", this code MUST be submitted such that seq first interval last matches the submitted forecast hours (f006 must be first=1, starting with f054 must be first=9), otherwise, the matching algorithm may attempt to match, for example, a 54 hour forecast with an analysis only at fhour 6.
#########################################

#set HEXT = ${1}

set ALLFHOURS = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
#set ALLFHOURS = ("f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")

set dispatch = "dispatch_matching.csh"

foreach f (`seq 1 1 40`) 
  set FHOUR = $ALLFHOURS[$f]  
  cp queue_matching_template queue_matching
  echo "./"$dispatch" "$f" "$FHOUR >> queue_matching
  echo $FHOUR
  qsubcasper queue_matching
  rm queue_matching
end  

