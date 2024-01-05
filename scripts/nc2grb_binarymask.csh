#!/bin/csh
#
#nc2grb_binarymask.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Auxilliary file to....
# 	Convert binary object .nc files to .grb files for MODE compatability
###########

set VARSTR = "HGT_500mb"
set DDIR = "/glade/work/klupo/postdoc/kasugaEA21/version7/"$VARSTR"/gridded/"
set SDIR = "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set MODEL = "gfs"
set RES = "0p25"
set EXT = "nc"
set EXT2 = "grb"

cd $DDIR     
      
set YEARS = ("2018")
set MONTHS = ("01")# "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

set D31 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
set D30 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
set D29 =  ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")
set D28 = ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28")
#set D31 = ("01")

set HOURS = ("00" "06" "12" "18")
set mm = "00"

foreach YYYY ($YEARS) #= "2020"
  if( $YYYY == "2012" || $YYYY == "2016" || $YYYY == "2020" )then
  set DL = ($D29)
  else
  set DL = ($D28)
  endif
foreach MM ($MONTHS) #= "09"
  if( $MM == "04" || $MM == "06" || $MM == "09" || $MM == "11" )then
  set DX = ($D30)
  endif
  if( $MM == "01" || $MM == "03" || $MM == "05" || $MM == "07" || $MM == "08" || $MM == "10" || $MM == "12" )then
  set DX = ($D31)
  endif
  if( $MM == "02" )then
  set DX = ($DL)
  endif

foreach DD ($DX) #DD = "09"
foreach hh ($HOURS) # hh = "00"

#Set forecast hour
foreach FHOUR ("f000" "f024" "f048" "f072")
 
#Set input and output names
set INFILE = $MODEL"."$RES"."$YYYY$MM$DD$hh"."$FHOUR".binary."$EXT
set OUTFILE = $MODEL"."$RES"."$YYYY$MM$DD$hh"."$FHOUR".binary."$EXT2

cdo -b 16 -f grb copy $INFILE $OUTFILE

end 
end
end 
end 
end       
    
