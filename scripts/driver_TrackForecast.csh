#!/bin/csh
#
# Kevin Lupo
# klupo@ucar.edu
#
# driver_TrackForecast.csh 
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
set SCRIPTDIR 	= /glade/u/home/klupo/postdoc/scripts/kasugaEA21
set SCRIPTDIR 	= $SCRATCH/cutofflow/scripts
set SCRATCHDIR 	= $SCRATCH/ks21_tmp
set SCRIPT	= run_TrackForecast

# ======== user set track/match params ===== #
set TCONFIGS = ("default" "pmax1.0" "pmax1.0_meannorms" "pmax2.0_meannorms" "pmax1.0_stdnorms" "pmax1.0_2stdnorms" "pmax1.0_meanPLSstdnorms" "pmax1.5_meannorms_munozDmax1200_oppmax700" "pmax1.5_2stdnorms_munozDmax1200_oppmax700")
set MCONFIGS = ("default" "pmax2.0_fixedradius_2stdnorms_equalweights_munozEmax800")
set SELECTTCONFIG = 9 #${1}
set SELECTMCONFIG = 2
set TCONFIG = $TCONFIGS[$SELECTTCONFIG]
set MCONFIG = $MCONFIGS[$SELECTMCONFIG]
echo "Tracking using the $TCONFIG configuration"
echo "Matching using the $MCONFIG configuration"

# ======= currently unused date lists ====== #
# Options for month length (For this test period, May uses $D31 and only 00 UTC) 
#set D31 	= ( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" )
#set D30 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
#set D29 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")
#set D28 	= ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28")


# ========== user set analysis vars ======== #
#set FSET	= "EnVAR"			# For efficiency, still submit ERA5 dates as individual jobs (ncks command in run script takes time)
#set FORECAST	= "MPAS15-3_15kmICs"
#set YEARS	= ("2019")
#set MONTHS	= ("05")
#set DAYS	= ( $D31 ) #("01" "02" "03" "04" "05" "06" "07" "08" "09" "10")
#set HOURS	= ( "00" )
#set FHOURS 	= ("f000")# "f006" "f012" "f018" "f024" "f030" "f036" "f042" "f048" "f054" "f060" "f066" "f072" "f078" "f084" "f090" "f096" "f102" "f108" "f114" "f120" "f126" "f132" "f138" "f144" "f150" "f156" "f162" "f168" "f174" "f180" "f186" "f192" "f198" "f204" "f210" "f216" "f222" "f228" "f234" "f240")
set YEARS = ("2020") #"2017" "2018" "2019" "2020" "2021")
set MONTHS =  ("09") # "02" "03" "04" "05" "06" "07" "08")# "09" "10" "11" "12")

# ========================================== #


foreach YEAR ($YEARS)							# For each user selected YEARS
  foreach MM ($MONTHS)							# For each user selected MONTHS
      
    #if( $YEAR == "2012" || $YEAR == "2016" || $YEAR == "2020" )then   # Get the appropriate month length for given year/month (May 2019 for this experimental period. Other options available. User will need to add other leap years if necessary)
    #  set DL = ($D29)
    #else
    #  set DL = ($D28)
    #endif
    #if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
    #  set DX = ($D30)
    #endif
    #if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
    #  set DX = ($D31)
    #endif
    #if( $MM == "02" )then
    #  set DX = ($DL)
    #endif
    #foreach DD ($DAYS)
    #  foreach hh ($HOURS)
	
        set JOBNAME = TrackForecast_$YEAR$MM					# Set the job name			
        set WORKDIR = $SCRATCHDIR/$YEAR$MM					# Set the working directory (in scratch space)
	  
        if ( ! -d $WORKDIR ) then						# Make the working directory if necessary
          mkdir -p $WORKDIR
        endif
        cd $WORKDIR								# Enter the working directory
        echo WORKDIR=$WORKDIR
        ln -sf $SCRIPTDIR/track_forecast .					# Symlink the identifation script (compiled fortran code) to the working directory
  
        if( -e ${WORKDIR}/$JOBNAME.log ) then					# Reset the log file if necessary
          rm ${WORKDIR}/$JOBNAME.log
        endif
    
        if ( $RUN_IN_PBS == "yes" ) then  								# Run in PBS queuing system
          echo "2i\"  									>! FTrack.sed
          echo "#==================================================================\" 	>> FTrack.sed
          echo "#PBS -N "$JOBNAME"\"						    	>> FTrack.sed
          echo "#PBS -j oe\"  						 		>> FTrack.sed
          echo "#PBS -o ${WORKDIR}/"$JOBNAME".log\"  				    	>> FTrack.sed
          echo "#PBS -A ${PROJ_NUMBER}\"						>> FTrack.sed
          echo "#PBS -q ${QUEUE}\"						   	>> FTrack.sed
          echo "#PBS -l walltime=${WALLTIME}\"						>> FTrack.sed
          echo "#PBS -l select=${N_NODES}:ncpus=${N_CPUS}:mem=${MEMORY}\"		>> FTrack.sed
          echo "#=================================================================="  	>> FTrack.sed
          echo 's%${1}%'"${YEAR}%g" 						    	>> FTrack.sed  	# Pass the year, month, day and hour to the "run" script
          echo 's%${2}%'"${MM}%g" 						    	>> FTrack.sed 
   	  echo 's%${3}%'"${TCONFIG}%g" 						    	>> FTrack.sed  
          echo 's%${4}%'"${MCONFIG}%g" 						    	>> FTrack.sed  
          echo 's%${5}%'"${WORKDIR}%g" 					    		>> FTrack.sed 
	
          sed -f FTrack.sed $SCRIPTDIR"/"$SCRIPT".csh" >! $SCRIPT".pbs"
          set jobid = `qsub $SCRIPT".pbs"`
          echo "${JOBNAME}:  ${jobid}"
          rm -rf $SCRIPT".pbs" FTrack.sed
          sleep 1
        else
	  cd $SCRIPTDIR
      set echo
	  ./run_TrackForecast.csh $YEAR $MM $TCONFIG $MCONFIG
      unset echo
	endif
  #    end #hh
 #   end #DD
  end #MM
end #YEAR


