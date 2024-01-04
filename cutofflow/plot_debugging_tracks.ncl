;debug_kasugaEA21.ncl
;;;;;;;;;;;;;
;Kevin Lupo
;klupo@ucar.edu
;verify output of KS21 algorithm against model height fields
;15 Oct 2021
;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Load;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"			; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"				; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"			; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"			; Standard load
load "../ncl/shared/CLIMO_VARS.ncl"							; Date lists for climo plots
external GRIDCOUNT "/glade/u/home/klupo/postdoc/scripts/fortran/update_grid_nclf90.so"

;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Begin;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

;Set model info
PARENT	= "/glade/collections/rda/data/ds084.1/"
MODEL	= "gfs"
RES	= "0p25"
EXT	= "grib2"
PLEV	= "500hPa"
DDIR	= "/glade/work/klupo/postdoc/kasugaEA21/version8/HGT_500mb/"
GDIR	= DDIR+"gridded/"
FDIR	= "/glade/u/home/klupo/postdoc/figures/noaa-jtti/kasugaEA21/debugmaps/"
;Set array subscripts
lat0	= 0
lat1	= 90
lon0	= 0
lon1	= 360
plevc	= 50000.

;Set initial time info & forecast hours
YYYY	= (/"2020"/)
MM	= (/"09"/);,"11"/);(/"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"/)
DD	= D31;(/"01","02","03"/);,"02"/);D31		; This is a dummy
hh	= HOURS
mm	= "00"
ny	= dimsizes(YYYY)
nm	= dimsizes(MM)
nd	= dimsizes(DD)
nh	= dimsizes(hh)
nt	= dimsizes(YYYY)*dimsizes(MM)*dimsizes(DD)*dimsizes(hh)
fhour	= (/000/)
nf	= dimsizes(fhour)

;Initialize arrays
DS	= new((/ny,nm,nd,nh,nf/),"integer",-9999)
ID	= new((/ny,nm,nd,nh,nf,200/),"integer",-1)
So	= new((/ny,nm,nd,nh,nf,200/),"float",-9999.9)
Slat	= So
Slon	= So
radius 	= So
Zmin	= So
Zlat	= So
Zlon	= So
sr	= So
vtime	= new((/ny,nm,nd,nh,nf,200/),"string","9999")
itime	= vtime

;TRUED(0) = 2
;;;;;;;;;;;;;;;Read data and set time-dependent parameters;;;;;;;;;;
ct	= 0

do y=0,ny-1
 
  if(mod(toint(YYYY(y)),4).eq.0)then
    TRUED(1) = 29
  else
    TRUED(1) = 28
  end if
  
  do m=0,nm-1
    print(YYYY(y)+MM(m))
    do d=25,TRUED(m)-1;nd-1
      do h=0,nh-1
        do f=0,nf-1
    
          valididx		= cd_inv_calendar(tointeger(YYYY(y)),tointeger(MM(m)),tointeger(DD(d)),tointeger(hh(h)),tointeger(mm),0,"hours after 1800-01-01 00:00:00",0)
          initidx 		= valididx-fhour(f)
          initidx@units 	= valididx@units
          VALID			= cd_calendar(valididx,3)
          INITIAL		= toint(cd_calendar(initidx,0))
          vtime(y,m,d,h,f,:) 	= VALID
          itime(y,m,d,h,f,:)	= cd_calendar(initidx,3)    

          FHOUR	= "f"+sprinti("%0.3i",fhour(f))
    
          ;KFILE needs to be the fhour from the initial time valid at vtime
          KFILE			= str_concat((/DDIR,MODEL,".",RES,".",sprinti("%0.4i",INITIAL(:,0)),sprinti("%0.2i",INITIAL(:,1)),sprinti("%0.2i",INITIAL(:,2)),sprinti("%0.2i",INITIAL(:,3)),".",FHOUR,".track"/))
          if(.not.fileexists(KFILE))then
	    print("File:"+KFILE+" is missing...skipping...")
            continue
          end if
          
	  data 				= asciiread(KFILE,-1,"string")
          ds				= dimsizes(data)-1
	  
          DS(y,m,d,h,f)			= ds
	  ID(y,m,d,h,f,0:ds-1)		= toint(str_get_field(data(1:dimsizes(data)-1),3," "))
          So(y,m,d,h,f,0:ds-1)		= tofloat(str_get_field(data(1:dimsizes(data)-1),4," "))	  
          Slat(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),5," "))
          Slon(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),6," "))
          radius(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),8," "))
          sr(y,m,d,h,f,0:ds-1)		= tofloat(str_get_field(data(1:dimsizes(data)-1),9," "))
          Zmin(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),13," "))
          Zlat(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),14," "))
          Zlon(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),15," "))

	  FHOUR	= "f"+sprinti("%0.3i",fhour(f))
          GDIR	= str_concat((/PARENT,sprinti("%0.4i",INITIAL(0,0)),"/",sprinti("%0.4i",INITIAL(0,0)),sprinti("%0.2i",INITIAL(0,1)),sprinti("%0.2i",INITIAL(0,2)),"/"/))
          FILE	= str_concat((/GDIR,MODEL,".",RES,".",sprinti("%0.4i",INITIAL(0,0)),sprinti("%0.2i",INITIAL(0,1)),sprinti("%0.2i",INITIAL(0,2)),sprinti("%0.2i",INITIAL(0,3)),".",FHOUR,".",EXT/))
          fi	= addfile(FILE,"r")
          if(ct.eq.0)then
            lon	= fi->lon_0({lon0:lon1})
            lat	= fi->lat_0({lat0:lat1})
	    z500			= new((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),"float","No_FillValue")
	    z500!5			= "lat"
	    z500!6			= "lon"
	    z500&lat			= lat
	    z500&lon			= lon
	   
	  end if
          z500(y,m,d,h,f,:,:)		= (/fi->HGT_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})/)
          ct = ct+1
	  delete([/data,fi/])
        end do
      end do
    end do
  end do
end do
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Data management tasks;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
z500 = z500/10.0	; Convert to dam

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Make Plots;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("plotting")
do y=0,ny-1
  if(mod(toint(YYYY(y)),4).eq.0)then
    TRUED(1) = 29
  else
    TRUED(1) = 28
  end if  
  do m=0,nm-1
    print(YYYY(y)+MM(m))
    do d=25,TRUED(m)-1;nd-1
      do h=0,nh-1
        do f=0,nf-1
	  dum  = new((/DS(y,m,d,h,f)/),"graphic")
	  dum2 = dum
	  dum3 = dum
	  dumc = dum

          load "../ncl/shared/MAPRES.ncl"
          wks_type	= poltype
          plot_type	= FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+itime(y,m,d,h,f,0)+"."+"f"+sprinti("%0.3i",fhour(f))+"_debugtracks_polar"
          wks		= gsn_open_wks(wks_type,plot_type)
          gsn_merge_colormaps(wks,"GMT_gray","precip2_15lev")
          res		= POLRES
          res2		= POLRES
	  res3		= POLRES
	  
	  res@gsnLeftString	= PLEV+" HGT (dam) | KS21 So [m 100km~S~-1~N~]~C~Init: "+itime(y,m,d,h,f,0)
          res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))+"~C~Valid: "+vtime(y,m,d,h,f,0)
          res@cnFillOn		= True
          res@cnLinesOn		= False
          res@cnLevelSelectionMode = "ExplicitLevels"
          res@cnLevels		= fspan(5.0,40.0,8)
          res@cnFillColors	= (/0,11,10,9,8,7,6,5,4/)
          res@cnConstFEnableFill= True
          res@lbLabelBarOn	= True
	  res@lbTitleString	= "So [m 100km~S~-1~N~]"
          plot			= gsn_csm_contour_map(wks,0.0,res)

          res2@cnLevelSelectionMode = "ExplicitLevels"
          res2@cnLevels		= fspan(420.0,660.0,41)
          res2@cnLineLabelsOn	= True
          res2@cnLineThicknessF	= 4.0
          res2@cnLineLabelFontHeightF = 0.015
          res2@cnLineDashSegLenF	= 0.25
          res2@cnLineDashPattern	= 0
          res2@cnLineLabelPlacementMode = "Constant"
          res2@cnLineLabelFont	= 22

          plot2			= gsn_csm_contour(wks,z500(y,m,d,h,f,:,:),res2)
          overlay(plot,plot2)
	  
	  
	  
          pmres			= True
          pmres@gsMarkerIndex	= 16
          pmres@gsMarkerColor 	= 1
          pmres@gsMarkerSizeF	= 0.01

          pmres2			= pmres
          pmres2@gsMarkerIndex	= 16
	  pmres2@gsMarkerSizeF	= 0.005
          pmres2@gsMarkerColor	= 14

          plres			= True
          plres@gsLineThicknessF= 6.0

          ;plot the base polymarker and the radial circles

          do i=0,DS(y,m,d,h,f)-1
	    if(ID(y,m,d,h,f,i).ne.32198)then
	    continue
	    end if
            dum(i) 		= gsn_add_polymarker(wks,plot,Slon(y,m,d,h,f,i),Slat(y,m,d,h,f,i),pmres)
            circle		= geolocation_circle(Slat(y,m,d,h,f,i),Slon(y,m,d,h,f,i),radius(y,m,d,h,f,i),1,180,False)
            circle_lat		= circle[0]
            circle_lon	 	= circle[1]
            if(.not.ismissing(Zmin(y,m,d,h,f,i)))then
              plres@gsLineColor       = "navy"
              dum3(i)	      = gsn_add_polymarker(wks,plot,Zlon(y,m,d,h,f,i),Zlat(y,m,d,h,f,i),pmres2)
            else
              plres@gsLineColor       = "forestgreen"
            end if
            dumc(i)	      = gsn_add_polyline(wks,plot,circle_lon(0,0,:),circle_lat(0,0,:),plres)
            delete([/circle,circle_lat,circle_lon/])
          end do

          pmres@gsMarkerSizeF	= 0.007
          do i=0,DS(y,m,d,h,f)-1
	    if(ID(y,m,d,h,f,i).ne.26852)then
	    continue
	    end if
            if(So(y,m,d,h,f,i).lt.10.0)then
              pmres@gsMarkerColor = 17
            else if(So(y,m,d,h,f,i).ge.10.0.and.So(y,m,d,h,f,i).lt.20.0)then
              pmres@gsMarkerColor = 22
            else if(So(y,m,d,h,f,i).ge.20.0.and.So(y,m,d,h,f,i).lt.30.0)then
              pmres@gsMarkerColor = 23
            else if(So(y,m,d,h,f,i).ge.30.0.and.So(y,m,d,h,f,i).lt.40.0)then
              pmres@gsMarkerColor = 26
            else if(So(y,m,d,h,f,i).ge.40.0)then
              pmres@gsMarkerColor = 27
            end if
            end if
            end if
            end if
            end if
            dum2(i) 		= gsn_add_polymarker(wks,plot,Slon(y,m,d,h,f,i),Slat(y,m,d,h,f,i),pmres)
          end do
	  draw(plot)
          frame(wks)
          system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
          delete([/res,res2,dum,dum2,dum3,dumc/])
         
	end do
      end do
    end do
  end do
end do



end
