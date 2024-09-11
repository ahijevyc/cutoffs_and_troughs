#!/bin/csh
#
# Kevin Lupo
# klupo@ucar.edu
#
# driver_TrackForecast.csh 
# 
#########################################

# =========== user set admin vars ========== #
set SCRIPTDIR 	= `pwd`

# ======== user set track/match params ===== #
set TCONFIGS = ("default" "pmax1.0" "pmax1.0_meannorms" "pmax2.0_meannorms" "pmax1.0_stdnorms" "pmax1.0_2stdnorms" "pmax1.0_meanPLSstdnorms" "pmax1.5_meannorms_munozDmax1200_oppmax700" "pmax1.5_2stdnorms_munozDmax1200_oppmax700")
set MCONFIGS = ("default" "pmax2.0_fixedradius_2stdnorms_equalweights_munozEmax800")
set SELECTTCONFIG = 9 #${1}
set SELECTMCONFIG = 2
set TCONFIG = $TCONFIGS[$SELECTTCONFIG]
set MCONFIG = $MCONFIGS[$SELECTMCONFIG]
echo "Tracking using the $TCONFIG configuration"
echo "Matching using the $MCONFIG configuration"

# ========== user set analysis vars ======== #
set YEARS = `seq 2020 2020`
set MONTHS =  `seq -w 09 09`

# ========================================== #

foreach YEAR ($YEARS)							# For each user selected YEARS
    foreach MM ($MONTHS)						# For each user selected MONTHS
        set WORKDIR = $SCRATCH/ks21_tmp/$YEAR$MM# Set the working directory (in scratch space)

        if ( ! -d $WORKDIR ) mkdir -p $WORKDIR
        cd $WORKDIR								# Enter the working directory
        echo WORKDIR=$WORKDIR
        ln -sf $SCRIPTDIR/track_forecast .		# Symlink the identifation script (compiled fortran code) to the working directory
        cd $SCRIPTDIR
        set echo
        $SCRIPTDIR/run_TrackForecast.csh $YEAR $MM $TCONFIG $MCONFIG
        unset echo
    end #MM
end #YEAR
