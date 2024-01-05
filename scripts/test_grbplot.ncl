begin

f	= addfile("~/work_new/postdoc/kasugaEA21/version7/HGT_500mb/gridded/gfs.0p25.2018013112.f000.binary.grb","r")
var1	= f->PRES_GDS0_SFC
var2	= f->PRMSL_GDS0_SFC

printVarSummary(var1)
printVarSummary(var2)


 wks = gsn_open_wks("png","test_grbplot")          ; send graphics to PNG file


  res                   = True
  res@mpProjection      = "Mollweide"       ; choose projection
  res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpGridLatSpacingF = 30.               ; spacing for lat lines
  res@mpGridLonSpacingF = 30.               ; spacing for lon lines
  res@mpFillOn          = False
 
  res@cnFillOn          = True              ; color plot desired
  res@cnFillPalette     = "gui_default"     ; set color map
  res@cnLineLabelsOn    = False             ; turn off contour lines
  res@txFontHeightF     = 0.015 
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels		= 0.5
  res@cnFillColors	= (/"white","red"/)
  res@vpXF            = 0.1                 ; make plot bigger
  res@vpYF            = 0.9         
  res@vpWidthF        = 0.8
  res@vpHeightF       = 0.8

  res@lbLabelFontHeightF  = 0.015           ; label bar font height

  res@tiMainString       = "Example of Mollweide Projection"  ; add a title
  res@tiMainFontHeightF  = .018                               ; font height

  contour = gsn_csm_contour_map(wks,var2,res)  ; create the plot


end
