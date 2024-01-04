#!/bin/csh
#
# Kevin Lupo
# klupo@ucar.edu
#
# driver_TrackForecast_cases.csh 
# 
# Created 14 July 2023
# 	--batch submits the run_TrackForecast script for given itimes & fhours
#
# Note...
#	--htc limits are 468 CPUs/4680 GB memory per user at any one time
#	
#	
#########################################

# =========== user set admin vars ========== #
set RUN_IN_PBS  = no
set PROJ_NUMBER = NMMM0021
set QUEUE       = casper
set N_NODES  	= 1
set N_CPUS      = 1
set MEMORY 	= "100GB"
set WALLTIME 	= "00:30:00"
set SCRIPTDIR 	= "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set SCRATCHDIR 	= "/glade/scratch/klupo/ks21_tmp/"
set SCRIPT	= "run_TrackForecast_cases"

# ======== user set track/match params ===== #
set TCONFIGS = ("default" "pmax1.0" "pmax1.0_meannorms" "pmax2.0_meannorms" "pmax1.0_stdnorms" "pmax1.0_2stdnorms" "pmax1.0_meanPLSstdnorms" "pmax1.5_meannorms_munozDmax1200_oppmax700" "pmax1.5_2stdnorms_munozDmax1200_oppmax700")
set MCONFIGS = ("default" "pmax2.0_fixedradius_2stdnorms_equalweights_munozEmax800")
set SELECTTCONFIG = 9 #${1}
set SELECTMCONFIG = 2
set TCONFIG = $TCONFIGS[$SELECTTCONFIG]
set MCONFIG = $MCONFIGS[$SELECTMCONFIG]
echo "Tracking using the $TCONFIG configuration"
echo "Matching using the $MCONFIG configuration"

set ITIMES = ("2019102206") 

# ========================================== #


foreach ITIME ($ITIMES)							# For each user selected YEARS
 
  set JOBNAME = "TrackForecast_"$ITIME				  # Set the job name			  
  set WORKDIR = $SCRATCHDIR"/"$ITIME  				  # Set the working directory (in scratch space)
    
  if ( ! -d $WORKDIR ) then						  # Make the working directory if necessary
    mkdir -p $WORKDIR
  endif
  cd $WORKDIR								  # Enter the working directory
  ln -sf $SCRIPTDIR/"track_forecast_cases" .  				  # Symlink the identifation script (compiled fortran code) to the working directory
  
  if( -e ${WORKDIR}/"$JOBNAME".log ) then				  # Reset the log file if necessary
    rm ${WORKDIR}/"$JOBNAME".log
  endif
  
  if ( $RUN_IN_PBS == "yes" ) then								  # Run in PBS queuing system
    echo "2i\"  								  >! FTrack.sed
    echo "#==================================================================\"   >> FTrack.sed
    echo "#PBS -N "$JOBNAME"\"  						  >> FTrack.sed
    echo "#PBS -j oe\"  							  >> FTrack.sed
    echo "#PBS -o ${WORKDIR}/"$JOBNAME".log\"					  >> FTrack.sed
    echo "#PBS -A ${PROJ_NUMBER}\"						  >> FTrack.sed
    echo "#PBS -q ${QUEUE}\"							  >> FTrack.sed
    echo "#PBS -l walltime=${WALLTIME}\"					  >> FTrack.sed
    echo "#PBS -l select=${N_NODES}:ncpus=${N_CPUS}:mem=${MEMORY}\"		  >> FTrack.sed
    echo "#=================================================================="    >> FTrack.sed
    echo 's%${1}%'"${ITIME}%g"							  >> FTrack.sed   # Pass the year, month, day and hour to the "run" script
    echo 's%${2}%'"${TCONFIG}%g"						  >> FTrack.sed  
    echo 's%${3}%'"${MCONFIG}%g"						  >> FTrack.sed  
    echo 's%${4}%'"${WORKDIR}%g"						  >> FTrack.sed 

    sed -f FTrack.sed $SCRIPTDIR"/"$SCRIPT".csh" >! $SCRIPT".pbs"
    set jobid = `qsubcasper $SCRIPT".pbs"`
    echo "${JOBNAME}:  ${jobid}"
    rm -rf $SCRIPT".pbs" FTrack.sed
    sleep 1
  else
    #cd $SCRIPTDIR
    ln -sf $SCRIPTDIR/"run_TrackForecast_cases.csh" . 
    ./run_TrackForecast_cases.csh $ITIME $TCONFIG $MCONFIG $WORKDIR
  endif
end #ITIME


