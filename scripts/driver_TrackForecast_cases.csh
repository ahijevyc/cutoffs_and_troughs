#!/bin/csh
#
# Kevin Lupo
# klupo@ucar.edu
#
# driver_TrackForecast_cases.csh 
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

# UFS cases 
set CASESDIR=/glade/campaign/mmm/parc/mwong/ufs-mrw

set isensemble=0

set FLENS = (024 048 072 240)  # zero-pad to match path name
if ($isensemble) set FLENS = (192 120)
foreach FLEN ($FLENS)
    set DDIRS = (`ls -d $CASESDIR/??????????.F$FLEN.C768`) # deterministic
    if ($isensemble) set DDIRS = (`ls -d $CASESDIR/E??????????.p??.F$FLEN.C768`) # ensemble

    foreach DDIR ($DDIRS)
        set WORKDIR = $SCRATCH/ks21_tmp/`basename $DDIR` # Set the working directory (in scratch space)

        if ( ! -d $WORKDIR ) mkdir -p $WORKDIR
        cd $WORKDIR								# Enter the working directory
        echo WORKDIR=$WORKDIR
        ln -sf $SCRIPTDIR/track_forecast_cases .  # Symlink the identification script (compiled fortran code) to the working directory
        $SCRIPTDIR/run_TrackForecast_cases.csh $DDIR $TCONFIG $MCONFIG F$FLEN
    end #ITIME
end #FLEN
