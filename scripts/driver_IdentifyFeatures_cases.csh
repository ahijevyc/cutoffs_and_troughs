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
set RUN_IN_PBS  = no
set PROJ_NUMBER = NMMM0021
set QUEUE       = casper
set N_NODES  	= 1
set N_CPUS      = 1
set MEMORY 	= "5GB"
set WALLTIME 	= "00:21:00"
set SCRIPTDIR 	= $SCRATCH/cutofflow/scripts
set SCRIPT	= run_IdentifyFeatures_cases
# ========== user set analysis vars ======== #
set IYYYYMMDDhh	= (2019102206 2019112306 2019121900 2020020806 2020021912 2020022900 2020040512 2020040812 2020051512 2020082818 2020090600 2020102418 2020112406 2021031318 2021040800 2021041012 2021042706 2021101518 2021102506 2021123006 2022021212 2022041006 2022050112 2022061018)
set FHOURS=`seq -w 006 6 240`
# ========================================== #

foreach FHOUR ($FHOURS)
  set FHOUR = f$FHOUR
  foreach ITIME ($IYYYYMMDDhh)
    set JOBNAME		= IdentifyFeatures_${ITIME}_$FHOUR		
    set WORKDIR 	= $SCRATCH/ks21_tmp/$ITIME # Used to have $FHOUR subdirectory too. Why? driver_TrackForecast_cases.csh can't find them.
    mkdir -vp $WORKDIR
    cd $WORKDIR
    echo WORKDIR=$WORKDIR
    ln -sf $SCRIPTDIR/identification_algorithm_globe_cases .
  
    if( -e ${WORKDIR}/$JOBNAME.log ) rm ${WORKDIR}/$JOBNAME.log
    
    if ( $RUN_IN_PBS == "yes" ) then  #  PBS queuing system
      echo "2i\"							      >! FIdent.sed
      echo "#==================================================================\" >> FIdent.sed
      echo "#PBS -N "$JOBNAME"\"					      >> FIdent.sed
      echo "#PBS -j oe\"						      >> FIdent.sed
      echo "#PBS -A ${PROJ_NUMBER}\"					      >> FIdent.sed
      echo "#PBS -q ${QUEUE}\"  					      >> FIdent.sed
      echo "#PBS -l walltime=${WALLTIME}\"				      >> FIdent.sed
      echo "#PBS -l select=${N_NODES}:ncpus=${N_CPUS}:mem=${MEMORY}\"	      >> FIdent.sed
      echo "#=================================================================="  >> FIdent.sed
      echo 's%${1}%'"${ITIME}%g" 					      >> FIdent.sed  
      echo 's%${2}%'"${FHOUR}%g"						      >> FIdent.sed  

      sed -f FIdent.sed $SCRIPTDIR/$SCRIPT.csh >! $SCRIPT.pbs
      set jobid = `qsub $SCRIPT.pbs`
      echo "${JOBNAME}:  ${jobid}"
      rm $SCRIPT.pbs FIdent.sed
      sleep 1

    else
    ln -sf $SCRIPTDIR/$SCRIPT.csh .
    set echo
    ./run_IdentifyFeatures_cases.csh $ITIME $FHOUR
    unset echo
    #/glade/u/home/klupo/scripts/Fall2021/run_stoch_traj.csh ${n} ${PERT_SCHEMED} ${CONFIG} ${CASE} ${MODE} >! ${WORK_DIR}/stoch_${n}.log

    endif
  end #ITIME
end  #f
