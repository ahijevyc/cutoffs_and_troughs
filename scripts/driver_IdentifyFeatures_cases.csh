#!/bin/csh
#
#klupo@ucar.edu
#
#driver_IdentifyFeatures_cases.csh 
#14 Dec 2022
#batch submits the queue_cutoffID script for given itimes & fhours
#8 Dec 2021
#keep # of jobs sent to queue relatively low, submit in groups of 4
#
#9 Dec 2021
#can go larger...htc limits are 468 CPUs/4680 GB memory per user at any one time
#...so, could technically submit all fhours simulataneously...
#...lets do 168 and 180 jobs (12 months x 14 and 15 hours)
#15 Sep 2023
#This version of the code is used to submit initial dates/fhours for UFS case study simulations
#########################################



# =========== user set admin vars ========== #
set SCRIPTDIR 	= $SCRATCH/cutofflow/scripts
set SCRIPT	= run_IdentifyFeatures_cases.csh
set CASESDIR=/glade/campaign/mmm/parc/mwong/ufs-mrw
# ========================================== #

set RES 	= C768				# Model res (C768, output is on 0.25latlon grid)
foreach FLEN (192 120) # zero-pad to match path name
    #set DDIRS = (`ls -d $CASESDIR/??????????.F$FLEN.$RES`) # deterministic
    set DDIRS = (`ls -d $CASESDIR/E??????????.p??.F$FLEN.$RES`) # ensemble members

    foreach DDIR ($DDIRS)
        set FHOURS=`seq -w 006 6 $FLEN`
        foreach FHOUR ($FHOURS)
            set WORKDIR = $SCRATCH/ks21_tmp/`basename $DDIR`
            mkdir -vp $WORKDIR
            cd $WORKDIR
            echo WORKDIR=$WORKDIR

            ln -sf $SCRIPTDIR/identification_algorithm_globe_cases $SCRIPTDIR/$SCRIPT .

            set cmd="./$SCRIPT $DDIR f$FHOUR F$FLEN"
            echo $cmd
            $cmd

        end
    end
end
