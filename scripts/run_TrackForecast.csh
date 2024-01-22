#!/bin/csh
#
###########
# Kevin Lupo
# klupo@ucar.edu
#
# run_TrackForecast.csh
#
# Created 14 July 2023
#######################################################

# Set some general vars
set start_time 	= `date +%s` 	# Record the start time for debugging/optimizing hours request

# Read the passed vars
set YYYY	= ${1}		# initial year
set MM		= ${2} 		# initial month
set TCONFIG	= ${3}  	# Configuration of the tracking scheme
set MCONFIG	= ${4}  	# Configuration of the matching scheme


########################### Block to configure the tracking scheme##############
set DEFAULT_PMAX = "1.5"	 # Default maximum normalized tracking penalty term
set DEFNORM_So = "15.0" 	 # Default normalization for So tracking penalty
set DEFNORM_BGo = "40.0"	 # Default normalization for BGo tracking penalty
set DEFNORM_Ro = "1000.0"	 # Default normalization for Ro tracking penalty
set DEFPREDERR = "1500.0"	 # Default maximum "predicted distance" error for tracking
set DEFOPPX    = "1000.0"	 # Default maximum zonal motion in opposite direction of predicted motion
set DEFOPPY    = "1000.0"	 # Default maximum merid motion in opposite direction of predicted motion
set DEFMAXDIST = "1500.0"	 # Default maximum distance travelled

if( $TCONFIG == "default" ) then
  set PMAX = $DEFAULT_PMAX
  set NORM_So = $DEFNORM_So
  set NORM_BGo = $DEFNORM_BGo
  set NORM_Ro = $DEFNORM_Ro
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.0" ) then
  set PMAX = "1.0"
  set NORM_So = $DEFNORM_So
  set NORM_BGo = $DEFNORM_BGo
  set NORM_Ro = $DEFNORM_Ro
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.0_meannorms" ) then
  set PMAX = "1.0"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax2.0_meannorms" ) then
  set PMAX = "2.0"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.0_stdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "6.8"
  set NORM_BGo = "8.1"
  set NORM_Ro = "316.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.0_2stdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "13.6"
  set NORM_BGo = "16.2"
  set NORM_Ro = "632.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.0_meanPLSstdnorms" ) then
  set PMAX = "1.0"
  set NORM_So = "24.4"
  set NORM_BGo = "22.4"
  set NORM_Ro = "904.0"
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif
if( $TCONFIG == "pmax1.5_meannorms_munozDmax1200_oppmax700" ) then
  set PMAX = "1.5"
  set NORM_So = "17.6"
  set NORM_BGo = "14.3"
  set NORM_Ro = "588.0"
  set EMAX = "1200.0"
  set OXMAX = "700.0"
  set OYMAX = "700.0"
  set DMAX = "1200.0"
endif

######## USE THIS ONE
if( $TCONFIG == "pmax1.5_2stdnorms_munozDmax1200_oppmax700" ) then
  set PMAX = "1.5"
  set NORM_So = "13.6"
  set NORM_BGo = "16.2"
  set NORM_Ro = "632.0"
  set EMAX = "1200.0"
  set OXMAX = "700.0"
  set OYMAX = "700.0"
  set DMAX = "1200.0"
endif
########################################################################



########################### Block to configure the MATCHING scheme##############
set DEF_MATCH_PMAX = "2.0"			# Default maximum matching penalty
set DEF_MATCHNORM_So = "6.0"			# Default normalization for So matching penalty
set DEF_MATCHNORM_BGo = "15.0"			# Default normalization for BGo matching penalty
set DEF_MATCHNORM_Ro = "750.0"			# Default normalization for Ro matching penalty
set DEF_MATCHNORM_DIST = "1500.0"		# Default normalization for distance matching penalty
set DEF_MATCHNORM_Z500    = "1000.0"		# Default normalization for Z500 matching penalty
set DEF_MATCHNORM_T500    = "1000.0"		# Default normalization for T500 matching penalty
set DEF_MATCHNORM_U500    = "1000.0"		# Default normalization for U500 matching penalty
set DEF_MATCHNORM_V500    = "1000.0"		# Default normalization for V500 matching penalty
set DEF_MATCHNORM_Q500    = "1000.0"		# Default normalization for Q500 matching penalty

if( $MCONFIG == "default" ) then
  set PMAX = $DEFAULT_PMAX
  set NORM_So = $DEFNORM_So
  set NORM_BGo = $DEFNORM_BGo
  set NORM_Ro = $DEFNORM_Ro
  set EMAX = $DEFPREDERR
  set OXMAX = $DEFOPPX
  set OYMAX = $DEFOPPY
  set DMAX = $DEFMAXDIST
endif

####### USE THIS ONE
if( $MCONFIG == "pmax2.0_fixedradius_2stdnorms_equalweights_munozEmax800" ) then
  set MATCH_PMAX = "2.0"
  set MATCH_NORM_So = "13.6"
  set MATCH_NORM_BGo = "16.2"
  set MATCH_NORM_Ro = "632.0"
  set MATCH_NORM_DIST = "800.0"
  set MATCH_NORM_Z500	= "543.2"
  set MATCH_NORM_T500	= "17.6"
  set MATCH_NORM_U500	= "16.6"
  set MATCH_NORM_V500	= "14.2"
  set MATCH_NORM_Q500	= "1.14"
endif

########################################################################

set VARSTR = "HGT_500mb"
set PARENT = /glade/work/klupo/postdoc/kasugaEA21/version9/$VARSTR
set MODEL = "gfs"
set RES = "0p25"
set EXT	= "track"
set EXT2 = "dat"
set DDIR = $PARENT
set VALIDLIST = $PARENT/f000.trackfiles.list


set D31 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
set D30 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
set D29 =  ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")
set D28 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28")
#set D31 = ("01" "02")
set D30 	= ("01")
set HOURS = ("00") # "06" "12" "18")
set mm = "00"

#Set initial time info and directory
if( $YYYY == "2012" || $YYYY == "2016" || $YYYY == "2020" )then
  set DL = ($D29)
else
  set DL = ($D28)
endif
if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
  set DX = ($D30)
endif
if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
  set DX = ($D31)
endif
if( $MM == "02" )then
  set DX = ($DL)
endif

foreach DD ( $DX )  
foreach hh ( $HOURS ) # hh = "00"

set VALIDSUBSET = $YYYY$MM$DD$hh.f000.V240list
set FXSUBSET = $YYYY$MM$DD$hh.fxlist
set TINDEX_ST = `grep -n $YYYY$MM$DD$hh.f000 $VALIDLIST | cut -f1 -d:`
set TINDEX_FN = `expr $TINDEX_ST + 40` #40`
set TINDEX_QT = `expr $TINDEX_FN + 1`

sed -n "${TINDEX_ST},${TINDEX_FN}p;${TINDEX_QT}q" $VALIDLIST > $VALIDSUBSET
ls $PARENT/$MODEL.$RES.$YYYY$MM$DD$hh.f*.dat > $FXSUBSET
set NV = `wc -l $VALIDSUBSET | cut -f1 -d" "`

echo Running $YYYY$MM$DD$hh
set echo
./track_forecast $VALIDSUBSET $FXSUBSET $NV $VARSTR $PMAX $NORM_So $NORM_BGo $NORM_Ro $EMAX $OXMAX $OYMAX $DMAX $MATCH_PMAX $MATCH_NORM_So $MATCH_NORM_BGo $MATCH_NORM_Ro $MATCH_NORM_DIST $MATCH_NORM_Z500 $MATCH_NORM_T500 $MATCH_NORM_U500 $MATCH_NORM_V500 $MATCH_NORM_Q500
unset echo

end 
end
set end_time = `date +%s`
echo "total time for "$YYYY" "$MM" was "`expr $end_time - $start_time`" seconds"

