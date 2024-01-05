#!/bin/csh
#
#split_hemispheres.csh
###########
#Kevin Lupo
#klupo@ucar.edu
#Auxilliary file to....
# 	split output files from KS21 algorithm into northern and southern hemisphere files
#	should make tracking more efficient
###########

set VARSTR = "HGT_500mb"
set DDIR = "/glade/work/klupo/postdoc/kasugaEA21/version7/"$VARSTR"/"
set SDIR = "/glade/u/home/klupo/postdoc/scripts/kasugaEA21/"
set MODEL = "gfs"
set RES = "0p25"
set EXT = "dat"
set EXTN = "NH"
set EXTS = "SH"
set META = `head -n 1 $DDIR/$MODEL.$RES.2015010100.f000.$EXT`
set METAF = "Metadata.line"
echo $META > $METAF


#Set datetime lists
set ITIMES = `cat Iyyyymmddhh.list.202203`
set FHOURS = `cat Fhhh.list`

cd $DDIR
foreach itime ($ITIMES)
  foreach fhour ($FHOURS)
    set FULLFILE = $MODEL.$RES.$itime.$fhour.$EXT
    set FULLLEN = `wc -l < "$FULLFILE"`
    #echo $FULLLEN
    if( $FULLLEN == "2" ) then
      #echo $FULLFILE" is a placeholder empty file. Copying to .SH and .NH" >> $SDIR"split_hemispheres.log"
      cp $FULLFILE "SHEM/"$FULLFILE.$EXTS
      cp $FULLFILE "NHEM/"$FULLFILE.$EXTN
    else
      #echo "Splitting "$FULLFILE >> $SDIR"split_hemispheres.log"
      set SIGNS = `cut -b 25 $FULLFILE`
      set LINE = "${#SIGNS}"
      set HLINE = `expr  $LINE + 1`
      csplit -s -f $FULLFILE $FULLFILE $HLINE
      mv $FULLFILE"00" "SHEM/"$FULLFILE.$EXTS
      cat $SDIR$METAF $FULLFILE"01" > "NHEM/"$FULLFILE.$EXTN
      rm $FULLFILE"01"
    endif
  end
end
       
    
