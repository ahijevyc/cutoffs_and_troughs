#!/bin/csh

set PARENT = "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set YEAR = ${1}
set MM = ${2}
set FHOUR = ${3}
set SDIR = "/glade/scratch/klupo/ks21_temp/"$YEAR"/"$MM"/"$FHOUR"/"


module load "grib-api"
module load "grib-bins"
module load "grib-libs"

if ( ! -d $SDIR) then
  mkdir -p $SDIR
endif

cd $SDIR
cp $PARENT"identification_wrapper_batchsubmit.csh" .
cp $PARENT"identification_algorithm_globe" .

./identification_wrapper_batchsubmit.csh $YEAR $MM $FHOUR
