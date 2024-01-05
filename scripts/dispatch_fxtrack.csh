#!/bin/csh

set PARENT = "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set YEAR = ${1}
set MM = ${2}
set TCONFIG = ${3}
set MCONFIG = ${4}
set SDIR = "/glade/scratch/klupo/ks21_temp/"$YEAR"/"$MM"/"

if ( ! -d $SDIR) then
  mkdir -p $SDIR
endif

cd $SDIR
cp $PARENT"track_forecast_wrapper.csh" .
cp $PARENT"track_forecast" .

./track_forecast_wrapper.csh $YEAR $MM $TCONFIG $MCONFIG
