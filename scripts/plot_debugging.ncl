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
DDIR	= "/glade/work/klupo/postdoc/kasugaEA21/version10/HGT_500mb/"
GDIR	= DDIR+"gridded/"
FDIR	= "/glade/u/home/klupo/postdoc/figures/noaa-jtti/kasugaEA21/version10/"
;Set array subscripts
lat0	= 0
lat1	= 90
lon0	= 0
lon1	= 360
plev	= 50000.

;Set initial time info & forecast hours
YYYY	= (/"2020"/)
MM	= (/"09"/);(/"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"/)
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
So	= new((/ny,nm,nd,nh,nf,200/),"float",-9999.9)
Slat	= So
Slon	= So
radius 	= So
Zmin	= So
Zlat	= So
Zlon	= So
sr	= So
Type	= So
vtime	= new((/ny,nm,nd,nh,nf,200/),"string","9999")
itime	= vtime

TRUED(0) = 2
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
    do d=8,8;0,TRUED(m)-1;nd-1
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
          KFILE			= str_concat((/DDIR,MODEL,".",RES,".",sprinti("%0.4i",INITIAL(:,0)),sprinti("%0.2i",INITIAL(:,1)),sprinti("%0.2i",INITIAL(:,2)),sprinti("%0.2i",INITIAL(:,3)),".",FHOUR,".dat"/))
          if(.not.fileexists(KFILE))then
	    print("File:"+KFILE+" is missing...skipping...")
            continue
          end if
          
	  data 				= asciiread(KFILE,-1,"string")
          ds				= dimsizes(data)-1
	  
          DS(y,m,d,h,f)			= ds
          So(y,m,d,h,f,0:ds-1)		= tofloat(str_get_field(data(1:dimsizes(data)-1),3," "))	  
          Slat(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),4," "))
          Slon(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),5," "))
          radius(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),6," "))
          sr(y,m,d,h,f,0:ds-1)		= tofloat(str_get_field(data(1:dimsizes(data)-1),7," "))
          Zmin(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),11," "))
          Zlat(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),12," "))
          Zlon(y,m,d,h,f,0:ds-1) 	= tofloat(str_get_field(data(1:dimsizes(data)-1),13," "))
          Type(y,m,d,h,f,0:ds-1)	= tofloat(str_get_field(data(1:dimsizes(data)-1),14," "))
	  DFILE 			= str_concat((/GDIR,MODEL,".",RES,".",sprinti("%0.4i",INITIAL(:,0)),sprinti("%0.2i",INITIAL(:,1)),sprinti("%0.2i",INITIAL(:,2)),sprinti("%0.2i",INITIAL(:,3)),".",FHOUR,".debug.nc"/))
          fi				= addfile(DFILE,"r")
	  if(ct.eq.0)then
	    dum				= fi->So(0,{lat0:lat1},{lon0:lon1})
            lon				= dum&longitude
            lat				= dum&latitude
	    z500			= new((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),"float","No_FillValue")
	    z500!5			= "lat"
	    z500!6			= "lon"
	    z500&lat			= lat
	    z500&lon			= lon
	    so500			= z500
	    bg500			= z500
	    roMSK			= z500
	  end if
          z500(y,m,d,h,f,:,:)		= (/fi->HGT_500mb(0,{lat0:lat1},{lon0:lon1})/)
	  so500(y,m,d,h,f,:,:)		= (/fi->So(0,{lat0:lat1},{lon0:lon1})/)
	  bg500(y,m,d,h,f,:,:)		= (/fi->BGo(0,{lat0:lat1},{lon0:lon1})/)
	  roMSK(y,m,d,h,f,:,:)		= (/fi->MASKRo(0,{lat0:lat1},{lon0:lon1})/)
          ct = ct+1
	  delete([/data/])
        end do
      end do
    end do
  end do
end do
delete(dum)
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
    do d=8,8;0,2;TRUED(m)-1;nd-1
      do h=0,nh-1
        do f=0,nf-1
	  ;Plot the optimal slope debugging map
	  
          load "../ncl/shared/MAPRES.ncl"
          wks_type	= poltype
          plot_type	= FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+itime(y,m,d,h,f,0)+"."+"f"+sprinti("%0.3i",fhour(f))+"_debugSo_polar"
          wks		= gsn_open_wks(wks_type,plot_type)
          gsn_define_colormap(wks,"BlueWhiteOrangeRed")
          res		= POLRES
          res2		= POLRES
	  res3		= POLRES
	  
	  res@gsnLeftString	= PLEV+" HGT (dam) | KS21 So [m 100km~S~-1~N~]~C~Init: "+itime(y,m,d,h,f,0)
          res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))+"~C~Valid: "+vtime(y,m,d,h,f,0)
          res@cnFillOn		= True
          res@cnLinesOn		= False
          res@cnLevelSelectionMode = "ExplicitLevels"
          res@cnLevels		= fspan(-25.0,25.0,21);fspan(10.0,24.0,8);fspan(5.0,40.0,8)/2.0
          res@cnFillColors	= (/2,12,22,32,42,52,62,72,82,92,102,155,165,175,185,195,205,215,225,235,245,255/)
          res@cnConstFEnableFill= True
          res@lbLabelBarOn	= True
	  res@lbTitleString	= "So [m 100km~S~-1~N~]"
	  res@lbTitleFont	= 22
	  res@lbTitleFontHeightF= 0.015
	  res@lbLabelStride	= 2
          plot			= gsn_csm_contour_map(wks,so500(y,m,d,h,f,:,:),res)

          res2@cnLevelSelectionMode = "ExplicitLevels"
          res2@cnLevels		= fspan(420.0,660.0,41)
          res2@cnLineLabelsOn	= True
          res2@cnLineThicknessF	= 6.0
          res2@cnLineLabelFontHeightF = 0.011
          res2@cnLineDashSegLenF	= 0.25
          res2@cnLineDashPattern	= 0
          res2@cnLineLabelPlacementMode = "Computed"
	  res2@cnLineLabelBackgroundColor = -1
          res2@cnLineLabelFont	= 22

          plot2			= gsn_csm_contour(wks,z500(y,m,d,h,f,:,:),res2)
          overlay(plot,plot2)	  
	  
	  
	  res3@cnLevelSelectionMode = "ExplicitLevels"
	  res3@cnLevels		= (/10.0/)
	  res3@cnFillOn		= False
	  res3@cnLinesOn	= True
	  delete(res3@cnLineColor)
	  res3@cnLineColor	= "darkviolet"
	  res3@cnLineThicknessF = 15.0
	  res3@cnLineLabelsOn	= False
	  res3@lbLabelBarOn	= False
	  plot3			= gsn_csm_contour(wks,so500(y,m,d,h,f,:,:),res3)
	  overlay(plot,plot3)
	  
	  draw(plot)
	  frame(wks)	  
          system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
          delete([/res,res2,res3/])
	  
	  
	  ;Plot the optimal background slope debugging map
	  
          load "../ncl/shared/MAPRES.ncl"
          wks_type	= poltype
          plot_type	= FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+itime(y,m,d,h,f,0)+"."+"f"+sprinti("%0.3i",fhour(f))+"_debugBG_polar"
          wks		= gsn_open_wks(wks_type,plot_type)
          gsn_define_colormap(wks,"BlueWhiteOrangeRed")
          res		= POLRES
          res2		= POLRES
	  res3		= POLRES
	  
	  res@gsnLeftString	= PLEV+" HGT (dam) | KS21 BGo [m 100km~S~-1~N~]~C~Init: "+itime(y,m,d,h,f,0)
          res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))+"~C~Valid: "+vtime(y,m,d,h,f,0)
          res@cnFillOn		= True
          res@cnLinesOn		= False
          res@cnLevelSelectionMode = "ExplicitLevels"
          res@cnLevels		= fspan(0.0,50.0,21);fspan(10.0,24.0,8);fspan(5.0,40.0,8)/2.0
	  cols			= ispan(150,255,5)
	  cols(0)		= 0
          res@cnFillColors	= cols
	  delete(cols)
          res@cnConstFEnableFill= True
          res@lbLabelBarOn	= True
	  res@lbTitleString	= "BGo [m 100km~S~-1~N~]"
	  res@lbTitleFont	= 22
	  res@lbTitleFontHeightF= 0.015
	  res@lbLabelStride	= 2
          plot			= gsn_csm_contour_map(wks,bg500(y,m,d,h,f,:,:),res)

          res2@cnLevelSelectionMode = "ExplicitLevels"
          res2@cnLevels		= fspan(420.0,660.0,41)
          res2@cnLineLabelsOn	= True
          res2@cnLineThicknessF	= 6.0
          res2@cnLineLabelFontHeightF = 0.011
          res2@cnLineDashSegLenF	= 0.25
          res2@cnLineDashPattern	= 0
          res2@cnLineLabelPlacementMode = "Computed"
	  res2@cnLineLabelBackgroundColor = -1
          res2@cnLineLabelFont	= 22

          plot2			= gsn_csm_contour(wks,z500(y,m,d,h,f,:,:),res2)
          overlay(plot,plot2)	  
	  draw(plot)
	  frame(wks)	  
          system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
          delete([/res,res2/])
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  ;Plot the diagnostic debugging map
	  
	  dum  = new((/DS(y,m,d,h,f)/),"graphic")
	  dum2 = dum
	  dum3 = dum
	  dumc = dum

          load "../ncl/shared/MAPRES.ncl"
          wks_type	= poltype
          plot_type	= FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+itime(y,m,d,h,f,0)+"."+"f"+sprinti("%0.3i",fhour(f))+"_debug_polar"
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
          res@cnLevels		= fspan(10.0,24.0,8);fspan(5.0,40.0,8)/2.0
          res@cnFillColors	= (/0,11,10,9,8,7,6,5,4/)
          res@cnConstFEnableFill= True
          res@lbLabelBarOn	= True
	  res@lbTitleString	= "So [m 100km~S~-1~N~]"
	  res@lbTitleFont	= 22
	  res@lbTitleFontHeightF= 0.015
          plot			= gsn_csm_contour_map(wks,so500(y,m,d,h,f,:,:),res)

          res2@cnLevelSelectionMode = "ExplicitLevels"
          res2@cnLevels		= fspan(420.0,660.0,41)
          res2@cnLineLabelsOn	= True
          res2@cnLineThicknessF	= 6.0
          res2@cnLineLabelFontHeightF = 0.011
          res2@cnLineDashSegLenF	= 0.25
          res2@cnLineDashPattern	= 0
          res2@cnLineLabelPlacementMode = "Computed"
	  res2@cnLineLabelBackgroundColor = -1
          res2@cnLineLabelFont	= 22

          plot2			= gsn_csm_contour(wks,z500(y,m,d,h,f,:,:),res2)
          overlay(plot,plot2)
	  
	  res3@cnLevelSelectionMode = "ExplicitLevels"
	  res3@cnLevels		= (/0.5/)
	  res3@cnFillOn		= True
	  res3@cnFillColors	= (/"transparent","red"/)
	  res3@cnFillOpacityF	= 0.1
	  res3@lbLabelBarOn	= False
	  res3@cnLinesOn	= False
	  plot3			= gsn_csm_contour(wks,roMSK(y,m,d,h,f,:,:),res3)
	  overlay(plot,plot3)
	  
          pmres			= True
          pmres@gsMarkerIndex	= 16
          pmres@gsMarkerColor 	= 1
          pmres@gsMarkerSizeF	= 0.01

          pmres2			= pmres
          pmres2@gsMarkerIndex	= 16
	  pmres2@gsMarkerSizeF	= 0.005
          pmres2@gsMarkerColor	= 14

          plres			= True
          plres@gsLineThicknessF= 10.0

          ;plot the base polymarker and the radial circles

          do i=0,DS(y,m,d,h,f)-1
            dum(i) 			= gsn_add_polymarker(wks,plot,Slon(y,m,d,h,f,i),Slat(y,m,d,h,f,i),pmres)
            circle			= geolocation_circle(Slat(y,m,d,h,f,i),Slon(y,m,d,h,f,i),radius(y,m,d,h,f,i),1,180,False)
            circle_lat			= circle[0]
            circle_lon	 		= circle[1]
            if(Type(y,m,d,h,f,i).eq.2.0)then
              plres@gsLineColor       	= "navy"
              dum3(i)	      		= gsn_add_polymarker(wks,plot,Zlon(y,m,d,h,f,i),Zlat(y,m,d,h,f,i),pmres2)
	    else if(Type(y,m,d,h,f,i).eq.3.0)then
              plres@gsLineColor       	= "darkviolet"
              dum3(i)	      		= gsn_add_polymarker(wks,plot,Zlon(y,m,d,h,f,i),Zlat(y,m,d,h,f,i),pmres2)
            else if(Type(y,m,d,h,f,i).eq.0.0)then
              plres@gsLineColor       	= "gray"
              dum3(i)	      		= gsn_add_polymarker(wks,plot,Zlon(y,m,d,h,f,i),Zlat(y,m,d,h,f,i),pmres2)
            else 
              plres@gsLineColor       	= "forestgreen"
            end if
	    end if
	    end if
            dumc(i)	      		= gsn_add_polyline(wks,plot,circle_lon(0,0,:),circle_lat(0,0,:),plres)
            delete([/circle,circle_lat,circle_lon/])
          end do
   
          mv			= (/10.0,15.0,20.0,25.0/)
          pmres@gsMarkerSizeF	= 0.007
          do i=0,DS(y,m,d,h,f)-1
            if(So(y,m,d,h,f,i).lt.mv(0))then
              pmres@gsMarkerColor = 17
            else if(So(y,m,d,h,f,i).ge.mv(0).and.So(y,m,d,h,f,i).lt.mv(1))then
              pmres@gsMarkerColor = 22
            else if(So(y,m,d,h,f,i).ge.mv(1).and.So(y,m,d,h,f,i).lt.mv(2))then
              pmres@gsMarkerColor = 23
            else if(So(y,m,d,h,f,i).ge.mv(2).and.So(y,m,d,h,f,i).lt.mv(3))then
              pmres@gsMarkerColor = 26
            else if(So(y,m,d,h,f,i).ge.mv(3))then
              pmres@gsMarkerColor = 27
            end if
            end if
            end if
            end if
            end if
            dum2(i) 		= gsn_add_polymarker(wks,plot,Slon(y,m,d,h,f,i),Slat(y,m,d,h,f,i),pmres)
          end do
	  draw(plot)
	  
	  
	  ; --- Resources for the colorbar/labelbar
	  LBW		 = 0.4 		; Width of the labelbar in NDC space	   
  	  LBH		 = 0.04 	   	; Height...
	  LBFontHeight	 = 0.015	   
	  LBFontFace	 = 22			; Helvetica bold
	  LBFontHeightTitle= 0.015
	  LBY		 = 0.08;0.18			; Increase, move label bar up 	   
	  LBX		 = 0.3;0.2125	   	; Increase, move label bar right
	  cnLevels        = mv          ; Coloring levels
  	  cnLabels        = sprintf("%3.0f",cnLevels)         ; Labels
	  lres                          = True
  	  lres@lbTitleString	      	= "Marker So [m 100km~S~-1~N~]"
	  lres@vpWidthF                 = LBW                    ; Width of Labelbar
  	  lres@vpHeightF                = LBH                   ; Height of Labelbar
 	  lres@lbPerimOn                = False                   ; Turn off perimeter.
  	  lres@lbOrientation            = "Horizontal"            ; Default is "Vertical"
	  lres@pmLabelBarOrthogonalPosF = -1.0
	  lres@lbLabelAlignment         = "InteriorEdges"         ; Default is "BoxCenters".
	  lres@lbFillColors             = (/17,22,23,26,27/)             ; Colors for boxes.
	  lres@lbMonoFillPattern        = True                    ; Fill them all solid.
	  lres@lbLabelFontHeightF       = LBFontHeight                  ; Font Height
	  lres@lbTitleFont		= LBFontFace
	  lres@lbLabelFont		= LBFontFace
	  lres@lbTitleFontHeightF 	= LBFontHeightTitle
	  lres@lbLabelAutoStride        = True                    ; Auto correct labels
	  lres@vpYF                     = LBY                   ; location of left edge
	  lres@vpXF                     = LBX                    ; location of top edge
	  gsn_labelbar_ndc(wks,dimsizes(cnLevels)+1,cnLabels,0.2,0.2,lres)
	  frame(wks)

	  
	  
          system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
          delete([/res,res2,dum,dum2,dum3,dumc/])
         
	end do
      end do
    end do
  end do
end do



end
