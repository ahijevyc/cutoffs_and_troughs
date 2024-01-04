;debug_kasugaEA21.ncl
;;;;;;;;;;;;;
;Kevin Lupo
;klupo@ucar.edu
;plot forecast errors during june 2019 before and after the transition to fv3
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
DDIR	= "/glade/work/klupo/postdoc/kasugaEA21/version7/HGT_500mb/"
GDIR	= DDIR+"gridded/"
FDIR	= "/glade/u/home/klupo/postdoc/figures/noaa-jtti/kasugaEA21/debugmaps/"
;Set array subscripts
lat0	= 0
lat1	= 90
lon0	= 0
lon1	= 360
plev	= 50000.

;Set initial time info & forecast hours
YYYY	= (/"2019"/)
MM	= (/"06"/);(/"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"/)
DD	= D31		; This is a dummy
hh	= HOURS
mm	= "00"
ny	= dimsizes(YYYY)
nm	= dimsizes(MM)
nd	= dimsizes(DD)
nh	= dimsizes(hh)
nt	= dimsizes(YYYY)*dimsizes(MM)*dimsizes(DD)*dimsizes(hh)
fhour	= ispan(0,6,6)
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
vtime	= new((/ny,nm,nd,nh,nf,200/),"string","9999")
itime	= vtime

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
    do d=0,TRUED(m)-1;nd-1
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
          ;if(toint(itime(y,m,d,h,f,0)).ge.2019061200.0)then
	  ;  print("fv3")
	  ;else
	    ;print(itime(y,m,d,h,f,0)+"")
	  ;  print("pre")
	  ;end if
	    
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

	  DFILE 			= str_concat((/GDIR,MODEL,".",RES,".",sprinti("%0.4i",INITIAL(:,0)),sprinti("%0.2i",INITIAL(:,1)),sprinti("%0.2i",INITIAL(:,2)),sprinti("%0.2i",INITIAL(:,3)),".",FHOUR,".debug.nc"/))
          if(.not.fileexists(DFILE))then
	    print("File:"+DFILE+" is missing...skipping...")
            delete(data)
	    continue
          end if
	  fi				= addfile(DFILE,"r")
	  if(ct.eq.0)then
	    dum				= fi->So(0,{lat0:lat1},{lon0:lon1})
            lon				= dum&longitude
            lat				= dum&latitude
	    z500			= new((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),"float",dum@_FillValue)
	    z500!5			= "lat"
	    z500!6			= "lon"
	    z500&lat			= lat
	    z500&lon			= lon
	    so500			= z500
	    bg500			= z500
	    ro500			= z500
	    roMSK			= z500
	  end if
          z500(y,m,d,h,f,:,:)		= (/fi->HGT_500mb(0,{lat0:lat1},{lon0:lon1})/)
	  so500(y,m,d,h,f,:,:)		= (/fi->So(0,{lat0:lat1},{lon0:lon1})/)
	  bg500(y,m,d,h,f,:,:)		= (/fi->BGo(0,{lat0:lat1},{lon0:lon1})/)
	  ro500(y,m,d,h,f,:,:)		= (/fi->Ro(0,{lat0:lat1},{lon0:lon1})/)
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
ro500 = ro500/1000.0	; convert to km
print("gridding itime")
itime_grid_str				= conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),itime(:,:,:,:,:,0),(/0,1,2,3,4/))
itime_grid			 	= toint(itime_grid_str)
print("computing errors")
Ez500					= z500-conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),z500(:,:,:,:,0,:,:),(/0,1,2,3,5,6/))
Eso500					= so500-conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),so500(:,:,:,:,0,:,:),(/0,1,2,3,5,6/))
Ebg500					= bg500-conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),bg500(:,:,:,:,0,:,:),(/0,1,2,3,5,6/))
Ero500					= ro500-conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),ro500(:,:,:,:,0,:,:),(/0,1,2,3,5,6/))
EroMSK					= roMSK-conform_dims((/ny,nm,nd,nh,nf,dimsizes(lat),dimsizes(lon)/),roMSK(:,:,:,:,0,:,:),(/0,1,2,3,5,6/))

print("averaging errors over appropriate itimes")
Ez500_pre				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.lt.2019061200,Ez500,z500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)
Ez500_fv3				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.ge.2019061200,Ez500,z500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)

Eso500_pre				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.lt.2019061200,Eso500,so500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)
Eso500_fv3				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.ge.2019061200,Eso500,so500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)

Ebg500_pre				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.lt.2019061200,Ebg500,bg500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)
Ebg500_fv3				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.ge.2019061200,Ebg500,bg500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)

Ero500_pre				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.lt.2019061200,Ero500,ro500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)
Ero500_fv3				= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.ge.2019061200,Ero500,ro500@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)

EroMSK_pre			     	= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.lt.2019061200,EroMSK,roMSK@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)
EroMSK_fv3			     	= smth9_Wrap(smth9_Wrap(dim_avg_n_Wrap(where(itime_grid.ge.2019061200,EroMSK,roMSK@_FillValue),(/0,1,2,3/)),0.5,0.25,True),0.5,0.25,True)

;print(EroMSK_fv3)
;print(Ero500_fv3)
;print(Ebg500_fv3)

;print(Ez500_pre)
;print(Ez500_fv3)

Ez500_pre!1				= "lat"
Ez500_pre!2				= "lon"
Ez500_pre&lat				= lat
Ez500_pre&lon				= lon
copy_VarMeta(Ez500_pre,Ez500_fv3)

copy_VarMeta(Ez500_pre,Eso500_pre)
copy_VarMeta(Ez500_pre,Eso500_fv3)

copy_VarMeta(Ez500_pre,Ebg500_pre)
copy_VarMeta(Ez500_pre,Ebg500_fv3)

copy_VarMeta(Ez500_pre,Ero500_pre)
copy_VarMeta(Ez500_pre,Ero500_fv3)

copy_VarMeta(Ez500_pre,EroMSK_pre)
copy_VarMeta(Ez500_pre,EroMSK_fv3)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Make Plots;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Making pre- & post-fv3 error plots")

do f=1,nf-1
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_z500_error"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" HGT Error (1-11 June 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.5,0.5,11)
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "HGT Error [dam]"
  plot			= gsn_csm_contour_map(wks,Ez500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" HGT Error (12-30 June 2019)"
  plot2			= gsn_csm_contour_map(wks,Ez500_fv3(f,:,:),res)
  
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_So500_error"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" So Error (1-11 June 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.4,0.4,11)
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "So Error [m 100km~S~-1~N~]"
  plot			= gsn_csm_contour_map(wks,Eso500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" So Error (12-30 June 2019)"
  plot2			= gsn_csm_contour_map(wks,Eso500_fv3(f,:,:),res)
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_Bg500_error"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" Bg Error (1-11 June 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.5,0.5,11)
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "Bg Error [m 100km~S~-1~N~]"
  plot			= gsn_csm_contour_map(wks,Ebg500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" Bg Error (12-30 June 2019)"
  plot2			= gsn_csm_contour_map(wks,Ebg500_fv3(f,:,:),res)
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_RoMSK_error"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" Feature Mask Error (1-11 June 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.05,0.05,11)
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "MSK Error"
  plot			= gsn_csm_contour_map(wks,EroMSK_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" Feature Mask Error (12-30 June 2019)"
  plot2			= gsn_csm_contour_map(wks,EroMSK_fv3(f,:,:),res)
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_Ro500_error"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" Ro Error (1-11 June 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-100.0,100.0,11)
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "Ro Error [km]"
  plot			= gsn_csm_contour_map(wks,Ero500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" Ro Error (12-30 June 2019)"
  plot2			= gsn_csm_contour_map(wks,Ero500_fv3(f,:,:),res)
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  
end do
exit





print("plotting")
do y=0,ny-1
  if(mod(toint(YYYY(y)),4).eq.0)then
    TRUED(1) = 29
  else
    TRUED(1) = 28
  end if  
  do m=0,nm-1
    print(YYYY(y)+MM(m))
    do d=0,TRUED(m)-1;nd-1
      do h=0,nh-1
        do f=0,nf-1
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
          res@cnLevels		= fspan(5.0,40.0,8)
          res@cnFillColors	= (/0,11,10,9,8,7,6,5,4/)
          res@cnConstFEnableFill= True
          res@lbLabelBarOn	= True
	  res@lbTitleString	= "So [m 100km~S~-1~N~]"
          plot			= gsn_csm_contour_map(wks,so500(y,m,d,h,f,:,:),res)

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
	  
	  res3@cnLevelSelectionMode = "ExplicitLevels"
	  res3@cnLevels		= (/0.5/)
	  res3@cnFillOn		= True
	  res3@cnFillColors	= (/"transparent","red"/)
	  res3@cnFillOpacityF	= 0.2
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
          plres@gsLineThicknessF= 6.0

          ;plot the base polymarker and the radial circles

          do i=0,DS(y,m,d,h,f)-1
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
