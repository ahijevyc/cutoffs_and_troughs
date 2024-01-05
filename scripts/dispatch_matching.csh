#!/bin/csh

set PARENT = "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set f = ${1}
set FHOUR = ${2}
#set HEXT = ${3}
set SDIR = "/glade/scratch/klupo/ks21_temp/"$FHOUR"/"

if ( ! -d $SDIR) then
  mkdir -p $SDIR
endif

cd $SDIR
cp $PARENT"match_forecast_wrapper_batchsubmit.csh" .
cp $PARENT"match_forecast" .

./match_forecast_wrapper_batchsubmit.csh $f $FHOUR 
