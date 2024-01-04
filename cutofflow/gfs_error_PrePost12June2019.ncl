;debug_kasugaEA21.ncl
;;;;;;;;;;;;;
;Kevin Lupo
;klupo@ucar.edu
;plot forecast errors during june 2019 before and after the transition to fv3
;25 Mar 2022
;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Load;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"			; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"				; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"			; Standard load
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"			; Standard load
load "../ncl/shared/CLIMO_VARS.ncl"							; Date lists for climo plots

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
plevc	= 50000.


;Set initial time info & forecast hours
YYYYMMDDHH_0 	= "2018060100"
yyyymmddhh_0 	= cd_inv_calendar(2018,6,1,0,0,0,"hours since 1800-01-01 00:00:00",0)
YYYYMMDDHH_1 	= "2019053118"
yyyymmddhh_1 	= cd_inv_calendar(2019,5,31,18,0,0,"hours since 1800-01-01 00:00:00",0)
nt	     	= toint((yyyymmddhh_1-yyyymmddhh_0)/6)+1

YYYYMMDDHH_2 	= "2020060100"
yyyymmddhh_2 	= cd_inv_calendar(2020,6,1,0,0,0,"hours since 1800-01-01 00:00:00",0)
YYYYMMDDHH_3 	= "2021053118"
yyyymmddhh_3 	= cd_inv_calendar(2021,5,31,18,0,0,"hours since 1800-01-01 00:00:00",0)

fhour		= ispan(0,72,72)
nf		= dimsizes(fhour)
do f=0,nf-1
  fstr		= "f"+sprinti("%0.3i",fhour(f))
  do i=0,nt-1
    vind	= yyyymmddhh_0+(6*i)	
    iind	= vind-fhour(f)
    iind@units	= "hours since 1800-01-01 00:00:00"
    init	= cd_calendar(iind,0)
    iyyyy	= toint(init(0,0))
    imm		= toint(init(0,1))
    idd		= toint(init(0,2))
    ihh		= toint(init(0,3))
    DDIR	= str_concat((/PARENT,sprinti("%0.4i",iyyyy),"/",sprinti("%0.4i",iyyyy),sprinti("%0.2i",imm),sprinti("%0.2i",idd),"/"/))
    FILE	= str_concat((/DDIR,MODEL,".",RES,".",sprinti("%0.4i",iyyyy),sprinti("%0.2i",imm),sprinti("%0.2i",idd),sprinti("%0.2i",ihh),".",fstr,".",EXT/))    
    
    if(.not.fileexists(FILE))then
      print("File:"+FILE+" is missing...skipping...")
      continue
    end if
    print(FILE+"")
    fi		= addfile(FILE,"r")

    if(i.eq.0.and.f.eq.0)then
      print("Initializing Arrays")
      lon	= fi->lon_0({lon0:lon1})
      lat	= fi->lat_0({lat0:lat1})
      z500_pre	= new((/nf,nt,dimsizes(lat),dimsizes(lon)/),"float",lat@_FillValue)
      u500_pre	= z500_pre
      v500_pre	= z500_pre
    end if
    z500_pre(f,i,:,:) = fi->HGT_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    u500_pre(f,i,:,:) = fi->UGRD_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    v500_pre(f,i,:,:) = fi->VGRD_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    delete(fi)
  end do
end do


do f=0,nf-1
  fstr		= "f"+sprinti("%0.3i",fhour(f))
  do i=0,nt-1
    vind	= yyyymmddhh_2+(6*i)	
    iind	= vind-fhour(f)
    iind@units	= "hours since 1800-01-01 00:00:00"
    init	= cd_calendar(iind,0)
    iyyyy	= toint(init(0,0))
    imm		= toint(init(0,1))
    idd		= toint(init(0,2))
    ihh		= toint(init(0,3))
    DDIR	= str_concat((/PARENT,sprinti("%0.4i",iyyyy),"/",sprinti("%0.4i",iyyyy),sprinti("%0.2i",imm),sprinti("%0.2i",idd),"/"/))
    FILE	= str_concat((/DDIR,MODEL,".",RES,".",sprinti("%0.4i",iyyyy),sprinti("%0.2i",imm),sprinti("%0.2i",idd),sprinti("%0.2i",ihh),".",fstr,".",EXT/))    
    
    if(.not.fileexists(FILE))then
      print("File:"+FILE+" is missing...skipping...")
      continue
    end if
    print(FILE+"")
    fi		= addfile(FILE,"r")

    if(i.eq.0.and.f.eq.0)then
      print("Initializing Arrays")
      lon	= fi->lon_0({lon0:lon1})
      lat	= fi->lat_0({lat0:lat1})
      z500_fv3	= new((/nf,nt,dimsizes(lat),dimsizes(lon)/),"float",lat@_FillValue)
      u500_fv3	= z500_fv3
      v500_fv3	= z500_fv3
    end if
    z500_fv3(f,i,:,:) = fi->HGT_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    u500_fv3(f,i,:,:) = fi->UGRD_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    v500_fv3(f,i,:,:) = fi->VGRD_P0_L100_GLL0({plevc},{lat0:lat1},{lon0:lon1})
    delete(fi)
  end do
end do

z500_pre	= z500_pre/10.0
z500_fv3	= z500_fv3/10.0

Ez500_pre	= z500_pre-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),z500_pre(0,:,:,:),(/1,2,3/))
AEz500_pre	= dim_avg_n_Wrap(Ez500_pre,1)
Ez500_fv3	= z500_fv3-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),z500_fv3(0,:,:,:),(/1,2,3/))
AEz500_fv3	= dim_avg_n_Wrap(Ez500_fv3,1)

Eu500_pre	= u500_pre-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),u500_pre(0,:,:,:),(/1,2,3/))
AEu500_pre	= smth9_Wrap(dim_avg_n_Wrap(Eu500_pre,1),0.5,0.25,True)
Eu500_fv3	= u500_fv3-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),u500_fv3(0,:,:,:),(/1,2,3/))
AEu500_fv3	= smth9_Wrap(dim_avg_n_Wrap(Eu500_fv3,1),0.5,0.25,True)

Ev500_pre	= v500_pre-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),v500_pre(0,:,:,:),(/1,2,3/))
AEv500_pre	= smth9_Wrap(dim_avg_n_Wrap(Ev500_pre,1),0.5,0.25,True)
Ev500_fv3	= v500_fv3-conform_dims((/nf,nt,dimsizes(lat),dimsizes(lon)/),v500_fv3(0,:,:,:),(/1,2,3/))
AEv500_fv3	= smth9_Wrap(dim_avg_n_Wrap(Ev500_fv3,1),0.5,0.25,True)

AEz500_pre!1	= "lat"
AEz500_pre!2	= "lon"
AEz500_pre&lat	= lat
AEz500_pre&lon	= lon
copy_VarMeta(AEz500_pre,AEz500_fv3)
copy_VarMeta(AEz500_pre,AEu500_pre)
copy_VarMeta(AEz500_pre,AEu500_fv3)
copy_VarMeta(AEz500_pre,AEv500_pre)
copy_VarMeta(AEz500_pre,AEv500_fv3)

;copy_VarMeta(Ez500_pre,Eso500_pre)
;copy_VarMeta(Ez500_pre,Eso500_fv3)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Make Plots;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Making pre- & post-fv3 error plots")

do f=1,nf-1
  ;if(0.gt.1)then
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_z500_error_june2018may2019_june2020may2021"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" HGT Error (1 June 2018-31 May 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.2,0.2,11)*5.0
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "HGT Error [dam]"
  plot			= gsn_csm_contour_map(wks,AEz500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" HGT Error (1 June 2020-31 May 2021)"
  plot2			= gsn_csm_contour_map(wks,AEz500_fv3(f,:,:),res)
  
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  ;end if
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_u500_error_june2018may2019_june2020may2021"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" UWIND Error (1 June 2018-31 May 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.2,0.2,11)*5.0
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "UWIND Error [m/s]"
  plot			= gsn_csm_contour_map(wks,AEu500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" UWIND Error (1 June 2020-31 May 2021)"
  plot2			= gsn_csm_contour_map(wks,AEu500_fv3(f,:,:),res)
  
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
  
  load "../ncl/shared/MAPRES.ncl"
  wks_type	= poltype
  plot_type     = FDIR+"kasugaEA21_"+MODEL+"."+RES+"."+PLEV+"_"+"f"+sprinti("%0.3i",fhour(f))+"_v500_error_june2018may2019_june2020may2021"
  wks	     	= gsn_open_wks(wks_type,plot_type)
  gsn_define_colormap(wks,"amwg_blueyellowred")
  res	      	= POLRES
  pres	      	= POPANELRES
  res@gsnLeftString	= PLEV+" VWIND Error (1 June 2018-31 May 2019)"
  res@gsnRightString	= "f"+sprinti("%0.3i",fhour(f))
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= fspan(-0.2,0.2,11)*5.0
  res@cnFillColors	= (/2,4,5,6,8,0,0,11,13,14,15,17/)
  res@cnConstFEnableFill= True
  res@lbLabelBarOn	= True
  res@lbTitleString	= "VWIND Error [m/s]"
  plot			= gsn_csm_contour_map(wks,AEv500_pre(f,:,:),res)

  res@gsnLeftString	= PLEV+" VWIND Error (1 June 2020-31 May 2021)"
  plot2			= gsn_csm_contour_map(wks,AEv500_fv3(f,:,:),res)
  
  pres@lbLabelBarOn	= False
  gsn_panel(wks,(/plot,plot2/),(/1,2/),pres)
  system("convert "+plot_type+"."+wks_type+" -trim "+plot_type+"."+wks_type+"")
  delete([/res,pres/])
end do

end 
