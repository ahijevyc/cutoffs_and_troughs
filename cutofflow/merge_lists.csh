#!/bin/csh
#Kevin Lupo klupo@ucar.edu
#4 May 2022
#merge_lists.csh
#utility to merge output lists from cutoff tracking or matching scheme into a single list for each fhour
#Generally helps to run this on an interactive casper node in the background with five instances (one for each list option)
#Submit interactive job as: execcasper -A NMMM0021 -l select=1:ncpus=5:mem=50GB
#
set LIST = ${1}			# See below for forecast hour grouping
set TYPE = ${2}			# track or match
set WDIR = "/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/longlists/"

#Note that this merge doesn't need to be done for match type, since f000 is by definition analysis tracking only. Only use would be for file name consistency in analysis code 
if( $LIST == "VALID" ) then
set FHOURS = ("f000")
endif

if( $LIST == "SHORT" ) then
set FHOURS = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048")
endif

if( $LIST == "MEDSHORT" ) then
set FHOURS = ("f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096")
endif

if( $LIST == "MEDIUM" ) then
set FHOURS = ("f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144")
endif

if( $LIST == "MEDLONG" ) then
set FHOURS = ("f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192")
endif

if( $LIST == "LONG" ) then
set FHOURS = ("f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
endif

if( $LIST == "092020TEST" ) then
set FHOURS = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168")
endif

cd $WDIR
foreach FHOUR ($FHOURS)
  echo "Merging $FHOUR $TYPE"
  cat "../"*"$FHOUR.$TYPE" > "./AllFeatures.$FHOUR.$TYPE.orig"
  grep -vwE "(ITIME)" "AllFeatures.$FHOUR.$TYPE.orig" > "AllFeatures.$FHOUR.$TYPE"
  #cat "../gfs.0p25.2020"*"$FHOUR.$TYPE" > "./AllFeatures.$FHOUR.$TYPE.orig"
  #grep -vwE "(ITIME)" "AllFeatures.$FHOUR.$TYPE.orig" > "AllFeatures.$FHOUR.$TYPE.2020.PSUMonly_pmax1.5_2stdnorms_munozDmax1200_oppmax700"
  rm "AllFeatures.$FHOUR.$TYPE.orig"
end


