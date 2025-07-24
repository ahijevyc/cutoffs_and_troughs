#!/bin/csh
#
#klupo@ucar.edu
#
#driver_IdentifyFeatures.csh 
#14 Dec 2022
#batch submits the queue_cutoffID script for given itimes & fhours
#8 Dec 2021
#keep # of jobs sent to queue relatively low, submit in groups of 4
#
#9 Dec 2021
#can go larger...htc limits are 468 CPUs/4680 GB memory per user at any one time
#...so, could technically submit all fhours simulataneously...
#...lets do 168 and 180 jobs (12 months x 14 and 15 hours)
#########################################


# =========== user set admin vars ========== #
set RUN_IN_PBS  = yes
set PROJ_NUMBER = NMMM0021
set QUEUE       = casper
set N_NODES  	= 1
set N_CPUS      = 1
set MEMORY 	= "5GB"
set WALLTIME 	= "00:21:00"
set SCRIPTDIR 	= `pwd`
set SCRIPT	= run_IdentifyFeatures
# ========== user set analysis vars ======== #
set YEARS	= ("2015" "2016" "2017" "2018" "2019" "2020" "2021" "2022" "2023")
set FHOURS	= (`seq -w 006 6 240`)
# ========================================== #

foreach FHOUR ($FHOURS) 
  set FHOUR=f$FHOUR
  foreach YEAR ($YEARS)
    
    if( $YEAR == "2023" )then
      set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08")
    else 
      set MONTHS = `seq -w 1 12`
    endif
  
    foreach MM ( $MONTHS )
      set JOBNAME	= IdentifyFeatures_$YEAR${MM}_$FHOUR		
      set WORKDIR 	= $SCRATCH/ks21_tmp/$YEAR/$MM/$FHOUR
      mkdir -vp $WORKDIR
      cd $WORKDIR
      echo WORKDIR=$WORKDIR
      ln -sf $SCRIPTDIR/identification_algorithm_globe .
  
      if( -e ${WORKDIR}/$JOBNAME.log ) rm ${WORKDIR}/$JOBNAME.log
    
      if ( $RUN_IN_PBS == "yes" ) then  #  PBS queuing system
        cat <<END | qsub
#!/bin/csh
#PBS -N $JOBNAME
#PBS -j oe
#PBS -A $PROJ_NUMBER
#PBS -q $QUEUE
#PBS -l walltime=$WALLTIME
#PBS -l select=${N_NODES}:ncpus=${N_CPUS}:mem=${MEMORY}
$SCRIPTDIR/$SCRIPT.csh $YEAR $MM $FHOUR
END
        echo $JOBNAME
        sleep 1

      else

      #/glade/u/home/klupo/scripts/Fall2021/run_stoch_traj.csh ${n} ${PERT_SCHEMED} ${CONFIG} ${CASE} ${MODE} >! ${WORK_DIR}/stoch_${n}.log

      endif
    end #MM
  end #YEAR
end  #f
