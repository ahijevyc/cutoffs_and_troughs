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
set SCRIPTDIR 	= "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set SCRIPT	= "run_IdentifyFeatures"
# ========== user set analysis vars ======== #
set YEARS	= ("2015" "2016" "2017" "2018" "2019" "2020" "2021" "2022" "2023")
#set YEARS	= ("2020")
#set FHOURS = ("f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
#set FHOURS = ("f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168")
set FHOURS = ("f000" "f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
set FHOURS	= ("f000")
# ========================================== #

foreach FHOUR ($FHOURS) 
  foreach YEAR ($YEARS)
    
    if( $YEAR == "2023" )then
      set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08")
    else 
      set MONTHS = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
    endif
  
    foreach MM ( $MONTHS )
      set JOBNAME	= "IdentifyFeatures_"$YEAR$MM"_"$FHOUR		
      set WORKDIR 	= "/glade/scratch/klupo/ks21_tmp/"$YEAR"/"$MM"/"$FHOUR"/"
      if ( ! -d $WORKDIR ) then
        mkdir -p $WORKDIR
      endif
      cd $WORKDIR
      ln -sf $SCRIPTDIR/"identification_algorithm_globe" .
  
      if( -e ${WORKDIR}/"$JOBNAME".log ) then
        rm ${WORKDIR}/"$JOBNAME".log
      endif
    
      if ( $RUN_IN_PBS == "yes" ) then  #  PBS queuing system
        echo "2i\"  								>! FIdent.sed
        echo "#==================================================================\" >> FIdent.sed
        echo "#PBS -N "$JOBNAME"\"						>> FIdent.sed
        echo "#PBS -j oe\"  							>> FIdent.sed
        echo "#PBS -o ${WORKDIR}/"$JOBNAME".log\"  				>> FIdent.sed
        echo "#PBS -A ${PROJ_NUMBER}\"						>> FIdent.sed
        echo "#PBS -q ${QUEUE}\"						>> FIdent.sed
        echo "#PBS -l walltime=${WALLTIME}\"					>> FIdent.sed
        echo "#PBS -l select=${N_NODES}:ncpus=${N_CPUS}:mem=${MEMORY}\"		>> FIdent.sed
        echo "#=================================================================="  >> FIdent.sed
        echo 's%${1}%'"${YEAR}%g" 						>> FIdent.sed  
        echo 's%${2}%'"${MM}%g" 						>> FIdent.sed  
        echo 's%${3}%'"${FHOUR}%g" 						>> FIdent.sed  

        sed -f FIdent.sed $SCRIPTDIR"/"$SCRIPT".csh" >! $SCRIPT".pbs"
        set jobid = `qsubcasper $SCRIPT".pbs"`
        echo "${JOBNAME}:  ${jobid}"
        rm -rf $SCRIPT".pbs" FIdent.sed
        sleep 1

      else

      #/glade/u/home/klupo/scripts/Fall2021/run_stoch_traj.csh ${n} ${PERT_SCHEMED} ${CONFIG} ${CASE} ${MODE} >! ${WORK_DIR}/stoch_${n}.log

      endif
    end #MM
  end #YEAR
end  #f