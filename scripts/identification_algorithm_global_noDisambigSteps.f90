!identification_algorithm_global_noDisambigSteps.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Kevin Lupo
!klupo@ucar.edu
!identify cutoff lows and preexisting troughs using kasugaEA21 algorithm
!
!29 Sep 2021
!	Computes slope function over number of radii for subset of northern hemisphere
!	Calculates distances using haversine function
!30 Sep 2021
!	Compute for NH and SH, omitting poles and tropics. rmax in KasugaEA21 does not require cross-pole gridpoint matching
!	Output slope function aslope(r,y,x) to netcdf for dev/debug
!	Maximize slope function wrt radius
!	Output So, Ro, BG
!	ID local maxima of So
!	Apply restrictions
!	ID gridpoints as features
!1 Oct 2021
!	Fixed bug with cyclic point
!	Output text file with feature info
!4 Oct 2021
!	Optimized search for features using pack() [appx 300x faster than searching all GPs]
!	Outputs [itime, fhour, So lat, lon, Ro, SR, zmin, zlat, zlon] 
!12 Oct 2021
!	Cleaned and added additional comments
! 	Corrected error that searched for local HEIGHT MAXIMA instead of MINIMA (this should correct some odd results comparing cutoffs and troughs)
! 	Modified criteria for nearby points - based on Slope Ratio now (keep smaller)
!14 Oct 2021
!	Added 9-point smoothing option for input height field
!	Nearby point criteria now based on So
!15 Oct 2021
!	Corrected error in linear interp of "lower" point at radius(r) from a point that incorrectly used "disti0" instead of "distj0"
!	Added "or" to distance check. If the distance between the points is less than EITHER optimal radius and one point is smaller than the other, wipe out the smaller So
!18 Oct 2021
!	Reduced Rmin to 100km
!19 Oct 2021
!	Set minimum So to 2.5 m/100km
!	Added smoothing option for So and SBG
!	Set Smoothing to 3-passes of a 9-point filter
!1 Nov 2021
!	Removed check to eliminate weak features within Ro of a greater slope feature. May help with tracking/assigning IDs
!	Fixed bug with finding height minima at cyclic point
!2 Nov 2021
!	Replaced "features within Ro" check with a criteria that determines if two cutoffs share a height minimum
!	.....Testing taking the larger Ro (Rationale being that the smaller feature is englufed by the larger)
!4 Nov 2021
!	Changed So cutoff to 5 to help tracking code (less likely that weak features will confuse later ID or tracking criteria)
!	If a trough is engulfed by the optimal radius of a cutoff, keep the cutoff
!	If a trough is engulfed by the optimal radius of a larger trough, keep the larger radius
!8 Nov 2021
! 	Wrapper script will now generate a dummy text file if model input is missing (should not have substantial effect on climo, but will have some tracking implications)
!	Fixed a bug that didn't output enough digits if So > 99.9 m/100km. This caused problems with tracking IO and prematurely terminated the output file for a given hour. Likely associated with TCs
!	For cleaner tracking, keep only features with So > 5.0 m/100km. This is consistent with the climatology of KasugaEA21
!11 Nov 2021
!	Changed output directory to ".../version5/..."
!	.dat files will include addition feature characteristics at So - 850-500-200 hPa HGT,TMP,RH,U,V
!17 Nov 2021
!	Not sure if this change is necessary, but adding an extra command line argument to increment the unit number of the output.dat file if running multiple instances of the wrapper script
!		...Should only need to make this change in this code, since tracking code cannot be run in parallel (and thus would not reference the same unit number concurrently)
!	Note that output erroneously duplicates SR into the BGo variable. Problem is noted and now accounted for in tracking code and ncl analysis
!7 Dec 2021
!	Modified submission code to batch submit ID algorithm to casper nodes
!8 Dec 2021
!	Feature data output to .dat file now includes area averages over Ro
!9 Dec 2021
!	Replaced RH output with qv mixing ratio
!10 Feb 2022
!	Added "debugonly" mode to only write netcdf files as output
!	...and "binary" for minimal output
!29 April 2022
!	Move to version 8
!	...In progress...
! 	Compute slope function in tropics - in progress...DONE 29 April 2022
!		- Simply removed if statement to avoid tropics
!	Compute slope function near poles - in progress...DONE 29 April 2022
!		- Maximum radius is limited to less than 0.5*circumference of earth at latitude to prevent longitude radial slope search from overlapping
!		- Cross polar radii are not tested
!		- Added array to store the maximum search radius at high latitudes (still need to eliminate features that have Ro=Rmax, even if Rmax is adjusted at latitudes)
!	Bypass out of domain tests for recording features...DONE 29 April 2022	
!	Remove non-KS21 disambiguation steps (shift this responsibility to tracking (consider tracking ONLY features with slope > 12 and plotting. Possible that original method was too preferentially selecting cutoffs over troughs)
!3 May 2022	
!	Final form 	- Smooth the height field 10x
!			- DO NOT smooth the resulting slopes	
!			- discard all features with a slope of less than 10 m/100km
!			- discard all features with a slope ratio greater than 2.25 (Necessary for noise)
!			- discard features with smaller So that share a height min with a larger So (Necessary for eroneously identified features)
!			- If features overlap by more than 90%, discard the feature with the smaller So (no preference for trough or cutoff)
!27 June 2022
!	Frustratingly, bug in the "Remove features at max-Ro search" step allowed features at the max alowed Ro to be retained (e.g., Ro=rmax). Moving to version 9
!		- Problem was that allowedR(:,:) was being set through "allowedR(y,:) = r", which uses an incorrect index ordering. Should be x,y (ie., allowedR(:,y) = r)
!		- This could yeild a variety of other errors. Needs to be rerun
!		- FIXED 27 June 2022
!	Since the code needs to be rerun to remove the ambiguous max Ro features & associated issues, taking this opportunity to add parameters for fixed-radius 500 hPa Z,qv,u,v, & T. 
!		- This is for matching later. Matching based on multi-level characteristics was redundant and confusing, and likely to cause issues for reviewers. Essentially, comparing Met vars over different radii would double-penalize differences in the radius (or artificially penalize trivial differences in the met fields).
!15 December 2022
!	Not changing anything substantial about Identification scheme (present results will not change)
!	Addding a Munoz EA20-like test to check for zonal wind reversal poleward of "cutoffs" 
!	If the mean zonal wind in area-averaged 2.75x2.75 (lat by lon) grid box poleward of "cutoff" So is less than zero (easterly/reversal), then point is a cutoff. Otherwise it is "closed"
!	Check avg(u(x0-5:x0+5,y-5:y+5)) starting at y=y0+6 ending at y=y0+35 (1.25 to 8.75 deg poleward of the cutoff).
! 	If any one of these 30 points has avg(u) less than zero - cutoff. Otherwise - closed.
! 	Moved to v10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
use netcdf

implicit none

real, parameter			:: rstep = 100000.0	![m] radii step (set to 100 km or 100000 m)
real, parameter			:: mres = 0.25		![deg] resolution of gfs data
real, parameter			:: SMP = 0.4		!Smoothing parameter p from NCL smth9
real, parameter			:: SMQ = 0.4		!Smoothing parameter q from NCL smth9
integer, parameter		:: NSM = 10		!Number of passes to smooth over
real, parameter			:: CEQUATOR = 40075000.0 ![m] circumference of earth at the equator
real				:: maxr, minr, srmax, CE, somin

integer				:: ncIDib, ncIDi, ncIDo, ncIDs, stat, istat, var_id, t_dimid, r_dimid, y_dimid, x_dimid, varid, incu, UNUM
integer			       	:: start(4), ct(4), dimids(4), ct_r(1), ct_x(1), ct_y(1), start_r(1), start_x(1), start_y(1), ct_t(1), start_t(1)
integer				:: start_nor(3), ct_nor(3), dimids_nor(3)
integer				:: nx, ny, nr, nt
integer				:: i, ii, j, k, x, y, r, a, b, np1, np2, np3, nz1, np, xm1, xp1
integer				:: i0, i1, j0, j1, xx0, xx1, yy0, yy1
integer				:: nrow
real				:: gpdistj, gpdisti, distj, disti, disti0, distj0, dist
real				:: xrslope, xlslope, yuslope, ylslope, left, right, upper, lower
real				:: d, haversine, hgtmin, c1c2_overlap, o, to_radian
real				:: Uavg
real				:: ST, FN, STsub, FNsub
real, allocatable		:: hgt_raw(:,:),hgt(:,:),aslopex(:,:), aslopexri(:,:),aslopexr(:,:),asmask(:,:),asmask2(:,:),asmask3(:,:),zmask(:,:),srat(:,:),BGx(:,:),BGxlat(:,:),BGxlon(:,:),cutflg(:,:),zval(:,:),zlat(:,:),zlon(:,:)
real, allocatable		:: z850(:,:),z200(:,:),t850(:,:),t500(:,:),t200(:,:),u850(:,:),u500(:,:),u200(:,:),v850(:,:),v500(:,:),v200(:,:),rh850(:,:),rh500(:,:),rh200(:,:)
real, allocatable		:: tc850(:,:),tc500(:,:),tc200(:,:),es850(:,:),es500(:,:),es200(:,:),mrs850(:,:),mrs500(:,:),mrs200(:,:),mr850(:,:),mr500(:,:),mr200(:,:)
real, allocatable		:: lat(:), lon(:), time(:), radius(:),lat2d(:,:),lon2d(:,:)
integer, allocatable		:: allowedR(:,:),allowedR_1D(:)
real, allocatable		:: aslope(:,:,:), BG(:,:,:), BGlat(:,:,:), BGlon(:,:,:)
real, allocatable		:: aslope_raw(:,:,:),BG_raw(:,:,:),BGlon_raw(:,:,:),BGlat_raw(:,:,:)
real, allocatable		:: features_p1(:,:),features_p2(:,:),features_p3(:,:),features_z(:,:),features(:,:)
real, allocatable		:: gridmask(:,:),gridmask_all(:,:)

!real, allocateable		:: hgt_crs(:,:)
character(len=80)		:: name

character(len=256)		:: infile
character(len=10)		:: maxr_c
character(len=10)		:: minr_c
character(len=5)		:: numr_c
character(len=5)		:: incu_c
character(len=5)		:: srmax_c
character(len=20)		:: varstr
character(len=5)		:: cyclic
character(len=256)		:: outfile,doutfile
character(len=350)		:: fstr
character(len=20)		:: itime
character(len=10)		:: fhour
character(len=10)		:: debug
character(len=256)		:: outdir
character(len=10)		:: smooth		!smoothing option - 9 point (smth9), none (none), area (area average)
character(len=10)		:: somin_c

!!!!!!Get some command line arguments

call get_command_argument(1,infile)
call get_command_argument(2,maxr_c)
read(maxr_c, *) maxr
call get_command_argument(3,minr_c)
read(minr_c, *) minr
call get_command_argument(4,numr_c)
read(numr_c, *) nr
call get_command_argument(5,srmax_c)
read(srmax_c, *) srmax
call get_command_argument(6,varstr)
call get_command_argument(7,cyclic)
call get_command_argument(8,debug)
call get_command_argument(9,incu_c)
read(incu_c, *) incu
call get_command_argument(10,smooth)
call get_command_argument(11,somin_c)
read(somin_c, *) somin
!!!!!Set the output unit number!!!!!
UNUM = 10+incu

!!!!!Name the output text file
outdir = "/glade/work/klupo/postdoc/kasugaEA21/version10/" // trim(varstr) // "/"
outfile = trim(outdir) // infile(31:54) // ".dat"

if(debug.eq."debug".or.debug.eq."debugonly")then
  doutfile = trim(outdir) // "gridded/" // infile(31:54) // "." // "debug" // ".nc"
endif
if(debug.eq."binary")then
  doutfile = trim(outdir) // "gridded/" // infile(31:54) // "." // "binary" // ".nc"
endif

itime = infile(40:49)
fhour = infile(51:54)


!!!!Set smoothing option!!!!!
!smooth = "area"	

call cpu_time(ST)

!!!!!!Open the netcdf file
call nc_check( nf90_open(trim(infile), nf90_nowrite, ncIDi), "kasugaEA21", "open " // trim(infile) )

!!!!!!Get nx and ny
call nc_check( nf90_inq_dimid(ncIDi, "longitude", var_id), "kasugaEA21", "inq_dimid longitude" )
call nc_check( nf90_inquire_dimension(ncIDi, var_id, name, nx), "kasugaEA21", "inquire_dimension longitude" )
call nc_check( nf90_inq_dimid(ncIDi, "latitude", var_id), "kasugaEA21", "inq_dimid latitude" )
call nc_check( nf90_inquire_dimension(ncIDi, var_id, name, ny), "kasugaEA21", "inquire_dimension latitude" )
call nc_check( nf90_inq_dimid(ncIDi, "time", var_id), "kasugaEA21", "inq_dimid time" )
call nc_check( nf90_inquire_dimension(ncIDi, var_id, name, nt), "kasugaEA21", "inquire_dimension time" )

!!!!!!Allocate arrays
allocate(aslope(nx,ny,nr),BG(nx,ny,nr),BGlat(nx,ny,nr),BGlon(nx,ny,nr))
allocate(aslope_raw(nx,ny,nr),BG_raw(nx,ny,nr),BGlat_raw(nx,ny,nr),BGlon_raw(nx,ny,nr))
allocate(hgt_raw(nx,ny),hgt(nx,ny),aslopex(nx,ny),aslopexri(nx,ny),aslopexr(nx,ny),asmask(nx,ny),asmask2(nx,ny),asmask3(nx,ny),zmask(nx,ny),srat(nx,ny),BGx(nx,ny),BGxlat(nx,ny),BGxlon(nx,ny))
allocate(z850(nx,ny),z200(nx,ny),t850(nx,ny),t500(nx,ny),t200(nx,ny),u850(nx,ny),u500(nx,ny),u200(nx,ny),v850(nx,ny),v500(nx,ny),v200(nx,ny),rh850(nx,ny),rh500(nx,ny),rh200(nx,ny))
allocate(tc850(nx,ny),tc500(nx,ny),tc200(nx,ny),es850(nx,ny),es500(nx,ny),es200(nx,ny),mrs850(nx,ny),mrs500(nx,ny),mrs200(nx,ny),mr850(nx,ny),mr500(nx,ny),mr200(nx,ny))
allocate(cutflg(nx,ny),zval(nx,ny),zlat(nx,ny),zlon(nx,ny))
allocate(lat(ny),lat2d(nx,ny))
allocate(lon(nx),lon2d(nx,ny))
allocate(time(nt))
allocate(gridmask(nx,ny),gridmask_all(nx,ny))
allocate(allowedR(nx,ny))
allocate(radius(nr))


!!!!!!Read lats and lons
call nc_check( nf90_inq_dimid(ncIDi, "latitude", var_id), "kasugaEA21", "inq_dimid latitude" )
call nc_check( nf90_get_var(ncIDi, var_id, lat), "kasugaEA21", "get_var lat" )
call nc_check( nf90_inq_dimid(ncIDi, "longitude", var_id), "kasugaEA21", "inq_dimid longitude" )
call nc_check( nf90_get_var(ncIDi, var_id, lon), "kasugaEA21", "get_var lon" )

call nc_check( nf90_inq_dimid(ncIDi, "time", var_id), "kasugaEA21", "inq_dimid time" )
call nc_check( nf90_get_var(ncIDi, var_id, time), "kasugaEA21", "get_var time" )

!!!!!!Make 2d lat lon arrays!!!!!
do x=1,nx
  lat2d(x,:) = lat
end do

do y=1,ny
  lon2d(:,y) = lon
end do

!!!!!!Read heights
call nc_check( nf90_inq_varid(ncIDi, trim(varstr), var_id), "kasugaEA21", "inq_varid " // trim(varstr) )
call nc_check( nf90_get_var(ncIDi, var_id, hgt), "kasugaEA21", "get_var hgt" ) 


!!!!!!Read auxilliary data
call nc_check( nf90_inq_varid(ncIDi, "HGT_850mb", var_id), "kasugaEA21", "inq_varid z850" )
call nc_check( nf90_get_var(ncIDi, var_id, z850), "kasugaEA21", "get_var z850" )

call nc_check( nf90_inq_varid(ncIDi, "HGT_200mb", var_id), "kasugaEA21", "inq_varid z200" )
call nc_check( nf90_get_var(ncIDi, var_id, z200), "kasugaEA21", "get_var z200" ) 

call nc_check( nf90_inq_varid(ncIDi, "TMP_850mb", var_id), "kasugaEA21", "inq_varid t850" )
call nc_check( nf90_get_var(ncIDi, var_id, t850), "kasugaEA21", "get_var t850" )

call nc_check( nf90_inq_varid(ncIDi, "TMP_500mb", var_id), "kasugaEA21", "inq_varid t500" )
call nc_check( nf90_get_var(ncIDi, var_id, t500), "kasugaEA21", "get_var t500" )

call nc_check( nf90_inq_varid(ncIDi, "TMP_200mb", var_id), "kasugaEA21", "inq_varid t200" )
call nc_check( nf90_get_var(ncIDi, var_id, t200), "kasugaEA21", "get_var t200" )

call nc_check( nf90_inq_varid(ncIDi, "UGRD_850mb", var_id), "kasugaEA21", "inq_varid u850" )
call nc_check( nf90_get_var(ncIDi, var_id, u850), "kasugaEA21", "get_var u850" )

call nc_check( nf90_inq_varid(ncIDi, "UGRD_500mb", var_id), "kasugaEA21", "inq_varid u500" )
call nc_check( nf90_get_var(ncIDi, var_id, u500), "kasugaEA21", "get_var u500" )

call nc_check( nf90_inq_varid(ncIDi, "UGRD_200mb", var_id), "kasugaEA21", "inq_varid u200" )
call nc_check( nf90_get_var(ncIDi, var_id, u200), "kasugaEA21", "get_var u200" ) 

call nc_check( nf90_inq_varid(ncIDi, "VGRD_850mb", var_id), "kasugaEA21", "inq_varid v850" )
call nc_check( nf90_get_var(ncIDi, var_id, v850), "kasugaEA21", "get_var v850" )

call nc_check( nf90_inq_varid(ncIDi, "VGRD_500mb", var_id), "kasugaEA21", "inq_varid v500" )
call nc_check( nf90_get_var(ncIDi, var_id, v500), "kasugaEA21", "get_var v500" )

call nc_check( nf90_inq_varid(ncIDi, "VGRD_200mb", var_id), "kasugaEA21", "inq_varid v200" )
call nc_check( nf90_get_var(ncIDi, var_id, v200), "kasugaEA21", "get_var v200" )

call nc_check( nf90_inq_varid(ncIDi, "RH_850mb", var_id), "kasugaEA21", "inq_varid rh850" )
call nc_check( nf90_get_var(ncIDi, var_id, rh850), "kasugaEA21", "get_var rh850" )

call nc_check( nf90_inq_varid(ncIDi, "RH_500mb", var_id), "kasugaEA21", "inq_varid rh500" )
call nc_check( nf90_get_var(ncIDi, var_id, rh500), "kasugaEA21", "get_var rh500" )

call nc_check( nf90_inq_varid(ncIDi, "RH_200mb", var_id), "kasugaEA21", "inq_varid rh200" )
call nc_check( nf90_get_var(ncIDi, var_id, rh200), "kasugaEA21", "get_var rh200" )

!!!!!!Close the netcdf file
call nc_check( nf90_close(ncIDi), "kasugaEA21", "close " // trim(infile) )
call cpu_time(FN)
print '("Time to Read = ",f6.3," seconds.")',FN-ST



!!!!!!Convert RH data to water vapor mixing ratio!!!!!!
!!!Convert T(K) to T(C)...
tc850 = t850 - 273.15
tc500 = t500 - 273.15
tc200 = t200 - 273.15

!!!Saturation vapor pressure... (Markowski and Richardson 2010, EQN 2.16, result is in hPa)
es850 = 6.112 * EXP( (17.67 * tc850)/(tc850 + 243.5) )
es500 = 6.112 * EXP( (17.67 * tc500)/(tc500 + 243.5) )
es200 = 6.112 * EXP( (17.67 * tc200)/(tc200 + 243.5) )

!!!Saturation mixing ratio... (weather.gov/media/epz/wxcalc/mixingRatio.pdf; req. es in hPa, returns g/kg)
mrs850 = 621.97 * (es850/(850.0-es850))
mrs500 = 621.97 * (es500/(500.0-es500))
mrs200 = 621.97 * (es200/(200.0-es200))

!!!Mixing ratio (from RH)
mr850 = (rh850 * mrs850) / 100.0
mr500 = (rh500 * mrs500) / 100.0
mr200 = (rh200 * mrs200) / 100.0


if(trim(smooth) == "smth9") then
call cpu_time(ST)
!!!!!!Smooth height field!!!!!!!!!!!9-point smoothing

hgt_raw = hgt

if(trim(cyclic) == 'no') then						!6 point smoothing at the east and west boundaries
do i=1,NSM  
  do y=2,ny-1								!Values at 90 and -90 are constant
    do x=2,nx-1								!Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
      hgt(x,y) = ((1.0-SMP-SMQ)*hgt_raw(x,y))+((SMP/4.0)*(hgt_raw(x,y+1)+hgt_raw(x+1,y)+hgt_raw(x,y-1)+hgt_raw(x-1,y)))+((SMQ/4.0)*(hgt_raw(x-1,y+1)+hgt_raw(x+1,y+1)+hgt_raw(x+1,y-1)+hgt_raw(x-1,y-1)))
    end do
    hgt(1,y)   = ((1.0-SMP-SMQ)*hgt_raw(1,y ))+((SMP/3.0)*(hgt_raw(1,y+1 )+hgt_raw(2,y   )+hgt_raw(1,y-1 )))+((SMQ/2.0)*(hgt_raw(2,y+1   )+hgt_raw(2,y-1   )))
    hgt(nx,y)  = ((1.0-SMP-SMQ)*hgt_raw(nx,y))+((SMP/3.0)*(hgt_raw(nx,y+1)+hgt_raw(nx,y-1)+hgt_raw(nx-1,y)))+((SMQ/2.0)*(hgt_raw(nx-1,y+1)+hgt_raw(nx-1,y-1)))
  end do
  hgt(:,1)  = hgt_raw(:,1)
  hgt(:,ny) = hgt_raw(:,ny)
  hgt_raw   = hgt
end do
end if

if(trim(cyclic) == 'yes') then						!9 point smoothing at the east and west boundaries by wrapping the cyclic point
do i=1,NSM
  do y=2,ny-1								!Values at 90 and -90 are constant
    do x=2,nx-1								!Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
      hgt(x,y) = ((1.0-SMP-SMQ)*hgt_raw(x,y))+((SMP/4.0)*(hgt_raw(x,y+1)+hgt_raw(x+1,y)+hgt_raw(x,y-1)+hgt_raw(x-1,y)))+((SMQ/4.0)*(hgt_raw(x-1,y+1)+hgt_raw(x+1,y+1)+hgt_raw(x+1,y-1)+hgt_raw(x-1,y-1)))
    end do
    hgt(1,y)  = ((1.0-SMP-SMQ)*hgt_raw(1,y))+((SMP/4.0)*(hgt_raw(1,y+1)+hgt_raw(2,y)+hgt_raw(1,y-1)+hgt_raw(nx,y)))+((SMQ/4.0)*(hgt_raw(nx,y+1)+hgt_raw(2,y+1)+hgt_raw(2,y-1)+hgt_raw(nx,y-1)))
    hgt(nx,y) = ((1.0-SMP-SMQ)*hgt_raw(nx,y))+((SMP/4.0)*(hgt_raw(nx,y+1)+hgt_raw(1,y)+hgt_raw(nx,y-1)+hgt_raw(nx-1,y)))+((SMQ/4.0)*(hgt_raw(nx-1,y+1)+hgt_raw(1,y+1)+hgt_raw(1,y-1)+hgt_raw(nx-1,y-1)))
  end do
  hgt_raw   = hgt
end do
end if
call cpu_time(FN)
print '("Time to Smooth HGT = ",f6.3," seconds.")',FN-ST
end if


if(trim(smooth) == "smth25") then
call cpu_time(ST)
!!!!!!Smooth height field!!!!!!!!!!!9-point smoothing

hgt_raw = hgt

if(trim(cyclic) == 'yes') then						!25 point smoothing at the east and west boundaries by wrapping the cyclic point
do i=1,1!NSM
  do y=3,ny-2								!Values at 90 and -90 are constant
    do x=3,nx-2								!Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
      hgt(x,y) = sum(hgt_raw(x-2:x+2,y-2:y+2))/25.0			!((1.0-SMP-SMQ)*hgt_raw(x,y))+((SMP/4.0)*(hgt_raw(x,y+1)+hgt_raw(x+1,y)+hgt_raw(x,y-1)+hgt_raw(x-1,y)))+((SMQ/4.0)*(hgt_raw(x-1,y+1)+hgt_raw(x+1,y+1)+hgt_raw(x+1,y-1)+hgt_raw(x-1,y-1)))
    end do
    
    hgt(2,y)	= (sum(hgt_raw(1:4,y-2:y+2))+sum(hgt_raw(nx,y-2:y+2)))/25.0
    hgt(1,y)	= (sum(hgt_raw(1:3,y-2:y+2))+sum(hgt_raw(nx-1:nx,y-2:y+2)))/25.0
    hgt(nx,y)	= (sum(hgt_raw(1:2,y-2:y+2))+sum(hgt_raw(nx-2:nx,y-2:y+2)))/25.0
    hgt(nx-1,y)	= (sum(hgt_raw(1,y-2:y+2))+sum(hgt_raw(nx-2:nx,y-2:y+2)))/25.0

  end do
  hgt_raw   = hgt
end do
end if
call cpu_time(FN)
print '("Time to Smooth HGT = ",f6.3," seconds.")',FN-ST
end if



if(trim(smooth) == "area") then						!Avoid. Slow
call cpu_time(ST)
!!!!!!Smooth height field!!!!!!!!!!!Area-averaging
hgt_raw = hgt
								
do y=1,ny								!Average all points within 150 km
  do x=1,nx
    gridmask = 0								
    do j=max(y-6,1),min(y+6,ny)
      do i=x-6,x+6
        if(i.lt.1)then
	  ii = ny-(abs(i))
	else if(i.gt.ny)then
	  ii = i-ny
	else 
	  ii = i
	end if
        
	!dist	= haversine(lat(y),lon(x),lat(j),lon(ii))
	!if(dist.le.150.0)then
	  gridmask(ii,j) = 1
	!else
	!  gridmask(i,j) = 0
	!end if
      end do
    end do
    hgt(x,y) = sum(hgt_raw,MASK=gridmask.eq.1)/sum(gridmask)
  end do
end do
call cpu_time(FN)
print '("Time to Smooth HGT = ",f6.3," seconds.")',FN-ST
end if





!!!!!!Calculate Slope Function!!!!!!
!Determine how many radii

call cpu_time(ST)

do r = 1, nr
  radius(r) = minr + ((r-1)*rstep)					! Set radius size
 			
  
  do y=1, ny
     CE = CEQUATOR*cos(to_radian(lat(y)))
     
     if(radius(r).ge.(0.5*CE)) cycle 					! At high latitudes, do not search radii that would overlap/circumnavigate the globe								
     if(radius(r).ge.(haversine(abs(lat(y)),10.0,90.0,10.0)*1000.0)) cycle		! At high latitudes, do not search radii that cross the pole
     allowedR(:,y) =  r							! Save the index of the max allowable radius at this latitude (used later to eliminate R=Rmax features)
    !if (abs(lat(y)).gt.89.75) cycle !	.lt.19.75.or.abs(lat(y)).gt.70.25) cycle		! Don't compute values poleward of 83 degrees (longitude slopes will overlap)
    
    									! Determine how many gridpoints to move to get to the specified radius 
    distj = haversine(lat(y),lon(1),lat(y+1),lon(1))*1000.0 		! Distance between gridpoints in the j/y direction. convert to meters
    gpdistj = (radius(r)/distj)						! Number of gridpoints in j/y to reach radius [m / (m/point) = points]					
    j0 = floor(gpdistj)							! Get integer gridpoints for linear interpolation
    j1 = ceiling(gpdistj)	
  
    
    									! Determine how many gridpoints to move to get to the specified radius 
    disti = haversine(lat(y),lon(1),lat(y),lon(2))*1000.0 		! Distance between gridpoints in the i/x direction. convert to meters
    gpdisti = (radius(r)/disti)  					! Number of gridpoints in i/x to reach radius [m / (m/point) = points]
    i0 = floor(gpdisti) 						! Get integer gridpoints for linear interpolation
    i1 = ceiling(gpdisti)
    
    do x=1, nx
      
      if (trim(cyclic) == 'no') then					! Data is not cyclic in longitude. Don't have to worry about cross polar points at the moment (maxr<20deg)
        if(((x+i0).ge.nx).or.((x-i0).le.1)) cycle			! Do not attempt to calculate slope if the test radius is outside of the domain bounds	
        								
									! Get HEIGHT at the radius distance to the EAST via linear interpolation between x+i0 and x+i1 [Point slope form y-y0 = m(x-x0); y=m(x-x0)+y0]
        disti0 = haversine(lat(y),lon(x),lat(y),lon(x+i0))*1000.0 	! x----------x+i0-x+i1
        xrslope = (hgt(x+i1,y)-hgt(x+i0,y))/disti			! |__disti0__|_disti_|
        right = (xrslope * (radius(r)-disti0)) + hgt(x+i0,y)		! |____radius____|____
    
    									! Get HEIGHT at the radius distance to the WEST via linear interpolation between x-i0 and x-i1 [Point slope form y-y0 = m(x-x0); y=m(x-x0)+y0]
        xlslope = (hgt(x-i1,y)-hgt(x-i0,y))/disti			! x-i1-x-i0----------x----------x+i0-x+i1
        left = (xlslope * (radius(r)-disti0)) + hgt(x-i0,y)		! |_disti_|__disti0__|
									! ____|____radius____|
      
        distj0 = haversine(lat(y),lon(x),lat(y+j0),lon(x))*1000.0	! Get HEIGHT at the radius distance to the NORTH via linear interpolation
        yuslope = (hgt(x,y+j1)-hgt(x,y+j0))/distj
        upper = (yuslope * (radius(r)-distj0)) + hgt(x,y+j0)
    
        ylslope = (hgt(x,y-j1)-hgt(x,y-j0))/distj			! Get HEIGHT at the radius distance to the SOUTH via linear interpolation
        lower = (ylslope * (radius(r)-distj0)) + hgt(x,y-j0)
      
        aslope(x,y,r) = (0.25)*(1/radius(r))*(right+left+upper+lower-(4.0*hgt(x,y)))*100000.0 	! Compute AS at radius (r) (Units of m/100km)
      									
									! Compute the background slope at radius (r)
        BGlat(x,y,r) = ((upper-lower)/(2*radius(r)))*100000.0 		! units of m/100km
        BGlon(x,y,r) = ((right-left)/(2*radius(r)))*100000.0 		! units of m/100km
        BG(x,y,r) = ((BGlat(x,y,r)**2 + BGlon(x,y,r)**2)**(0.5))	! *100000.0 		! units of m/100km
      endif
      
      
      
      if (trim(cyclic) == 'yes') then					! Data IS cyclic in longitude. Still don't have to worry about cross polar points at the moment (maxr<20deg). 
      									! Computation is similar to the non-cyclic case, but with adjusted gridpoints. See comments on above section.
      									! Testing in debug mode should show smooth/periodic longitudinal boundaries
	  
        if(((x+i0).ge.nx)) then						! Must consider case where the radius spans the cyclic longitude point TO THE RIGHT
	  
	  if((x+i0).gt.nx) then
	    xx0 = ((x+i0) - nx)						! if x+FLOOR is larger than nx, set xx0 to be HOW MUCH LARGER x+FLOOR is than nx
	    xx1 = ((x+i1) - nx)						! Similar for the ceiling
	  endif
	  if((x+i0).eq.nx) then
	    xx0 = nx							! if x+FLOOR is equal to nx, set xx0 to be nx
	    xx1 = 1							! set the ceiling point to be the first i index
	  endif
	  
	  disti0 = haversine(lat(y),lon(x),lat(y),lon(xx0))*1000.0 	! convert to meters
          xrslope = (hgt(xx1,y)-hgt(xx0,y))/disti
          right = (xrslope * (radius(r)-disti0)) + hgt(xx0,y)
    
	  xlslope = (hgt(x-i1,y)-hgt(x-i0,y))/disti
          left = (xlslope * (radius(r)-disti0)) + hgt(x-i0,y)
	  
	else if(((x-i0).le.1)) then					! Must consider case where the radius spans the cyclic longitude point TO THE LEFT

	  if((x-i0).lt.1) then					
	    xx0 = nx + (x-i0)						! if x-FLOOR is smaller than 1, set xx0 to be nx minus how much less than 1 (x-i0 is zero or negative)
	    xx1 = nx + (x-i1)						! similar for ceiling
	  endif
	  if((x-i0).eq.1) then
	    xx0 = 1							! if x-FLOOR equals 1, set xx0 to be 1
	    xx1 = nx							! set the ceiling point to be the final i index
	  endif
	  
	  disti0 = haversine(lat(y),lon(x),lat(y),lon(x+i0))*1000.0 	! convert to meters
          xrslope = (hgt(x+i1,y)-hgt(x+i0,y))/disti
          right = (xrslope * (radius(r)-disti0)) + hgt(x+i0,y)
	  
          xlslope = (hgt(xx1,y)-hgt(xx0,y))/disti
          left = (xlslope * (radius(r)-disti0)) + hgt(xx0,y)
	
	else
      
          disti0 = haversine(lat(y),lon(x),lat(y),lon(x+i0))*1000.0 	! convert to meters
          xrslope = (hgt(x+i1,y)-hgt(x+i0,y))/disti
          right = (xrslope * (radius(r)-disti0)) + hgt(x+i0,y)
     
          xlslope = (hgt(x-i1,y)-hgt(x-i0,y))/disti
          left = (xlslope * (radius(r)-disti0)) + hgt(x-i0,y)
        
	endif
	
	
	if((y-j0).le.1) then
	  
	  if((y-j0).lt.1) then
	    yy0 = 2+abs(y-j0)
	    yy1 = 2+abs(y-j1)
	  endif
	  if((y-j0).eq.1) then
	    yy0 = 1
	    yy1 = 2
	  endif
	  
	  if(lon(x).gt.180) then					! need to get point at x on other side of sphere
	    xx0 = x - (1/mres)*180
	  else if(lon(x).le.180) then
	    xx0 = x + (1/mres)*180
	  endif
          
	  distj0 = haversine(lat(y),lon(x),lat(y+j0),lon(x))*1000.0 	! convert to meters
          yuslope = (hgt(x,y+j1)-hgt(x,y+j0))/distj
          upper = (yuslope * (radius(r)-distj0)) + hgt(x,y+j0)
  	  
          ylslope = (hgt(xx0,yy1)-hgt(xx0,yy0))/distj
          lower = (ylslope * (radius(r)-distj0)) + hgt(xx0,yy0)
	  
	else if((y+j0).ge.ny) then
	  if((y+j0).gt.ny) then
	    yy0 = ny-(j0+y-ny)
	    yy1 = ny-(j1+y-ny)
	  endif
	  if((y+j0).eq.ny) then
	    yy0 = ny
	    yy1 = ny-1
	  endif
	  
	  if(lon(x).gt.180) then					! need to get point at x on other side of sphere
	    xx0 = x - (1/mres)*180
	  else if(lon(x).le.180) then
	    xx0 = x + (1/mres)*180
	  endif
          
	  distj0 = haversine(lat(y),lon(x),lat(y-j0),lon(x))*1000.0 	! convert to meters
          yuslope = (hgt(xx0,yy1)-hgt(xx0,yy0))/distj
          upper = (yuslope * (radius(r)-distj0)) + hgt(xx0,yy0)
  	  
          ylslope = (hgt(x,y-j1)-hgt(x,y-j0))/distj
          lower = (ylslope * (radius(r)-distj0)) + hgt(x,y-j0)
	
	else

	
	
	
        distj0 = haversine(lat(y),lon(x),lat(y+j0),lon(x))*1000.0 	! convert to meters
        yuslope = (hgt(x,y+j1)-hgt(x,y+j0))/distj
        upper = (yuslope * (radius(r)-distj0)) + hgt(x,y+j0)
    
        ylslope = (hgt(x,y-j1)-hgt(x,y-j0))/distj
        lower = (ylslope * (radius(r)-distj0)) + hgt(x,y-j0)
      end if
        aslope(x,y,r) = (0.25)*(1/radius(r))*(right+left+upper+lower-(4.0*hgt(x,y)))*100000.0 	! Units of m/100km
      
        BGlat(x,y,r) = ((upper-lower)/(2*radius(r)))*100000.0 		! units of m/100km
        BGlon(x,y,r) = ((right-left)/(2*radius(r)))*100000.0 		! units of m/100km
        BG(x,y,r) = ((BGlat(x,y,r)**2 + BGlon(x,y,r)**2)**(0.5))	! *100000.0 		! units of m/100km
      endif
            
end do ; end do ; end do
call cpu_time(FN)
print '("Time to compute AS & SBG= ",f6.3," seconds.")',FN-ST

!!!!!!!!!!!!!Smooth the slope function!!!!!!!!!
if(trim(smooth) == "smth9") then
call cpu_time(ST)
!!!!!!Smooth Slope field!!!!!!!!!!!9-point smoothing

aslope_raw = aslope
BG_raw     = BG
BGlat_raw  = BGlat
BGlon_raw  = BGlon

if(trim(cyclic) == 'no') then						!6 point smoothing at the east and west boundaries
do i=1,NSM  
  do r=1,nr
    do y=2,ny-1								
      if (abs(lat(y)).lt.20.or.abs(lat(y)).gt.70) cycle			! Don't compute values poleward of 70 degrees or between 20S and 20N
      do x=2,nx-1							! Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
        aslope(x,y,r) = ((1.0-SMP-SMQ)*aslope_raw(x,y,r))+((SMP/4.0)*(aslope_raw(x,y+1,r)+aslope_raw(x+1,y,r)+aslope_raw(x,y-1,r)+aslope_raw(x-1,y,r)))+((SMQ/4.0)*(aslope_raw(x-1,y+1,r)+aslope_raw(x+1,y+1,r)+aslope_raw(x+1,y-1,r)+aslope_raw(x-1,y-1,r)))
        
	BG(x,y,r) = ((1.0-SMP-SMQ)*BG_raw(x,y,r))+((SMP/4.0)*(BG_raw(x,y+1,r)+BG_raw(x+1,y,r)+BG_raw(x,y-1,r)+BG_raw(x-1,y,r)))+((SMQ/4.0)*(BG_raw(x-1,y+1,r)+BG_raw(x+1,y+1,r)+BG_raw(x+1,y-1,r)+BG_raw(x-1,y-1,r)))
        
	BGlat(x,y,r) = ((1.0-SMP-SMQ)*BGlat_raw(x,y,r))+((SMP/4.0)*(BGlat_raw(x,y+1,r)+BGlat_raw(x+1,y,r)+BGlat_raw(x,y-1,r)+BGlat_raw(x-1,y,r)))+((SMQ/4.0)*(BGlat_raw(x-1,y+1,r)+BGlat_raw(x+1,y+1,r)+BGlat_raw(x+1,y-1,r)+BGlat_raw(x-1,y-1,r)))

	BGlon(x,y,r) = ((1.0-SMP-SMQ)*BGlon_raw(x,y,r))+((SMP/4.0)*(BGlon_raw(x,y+1,r)+BGlon_raw(x+1,y,r)+BGlon_raw(x,y-1,r)+BGlon_raw(x-1,y,r)))+((SMQ/4.0)*(BGlon_raw(x-1,y+1,r)+BGlon_raw(x+1,y+1,r)+BGlon_raw(x+1,y-1,r)+BGlon_raw(x-1,y-1,r)))
      end do
      aslope(1,y,r)   = ((1.0-SMP-SMQ)*aslope_raw(1 ,y,r))+((SMP/3.0)*(aslope_raw(1 ,y+1,r)+aslope_raw(2 ,y  ,r)+aslope_raw(1 ,y-1,r)))+((SMQ/2.0)*(aslope_raw(2   ,y+1,r)+aslope_raw(2   ,y-1,r)))
      aslope(nx,y,r)  = ((1.0-SMP-SMQ)*aslope_raw(nx,y,r))+((SMP/3.0)*(aslope_raw(nx,y+1,r)+aslope_raw(nx,y-1,r)+aslope_raw(nx-1,y,r)))+((SMQ/2.0)*(aslope_raw(nx-1,y+1,r)+aslope_raw(nx-1,y-1,r)))
 
      BG(1,y,r)   = ((1.0-SMP-SMQ)*BG_raw(1 ,y,r))+((SMP/3.0)*(BG_raw(1 ,y+1,r)+BG_raw(2 ,y  ,r)+BG_raw(1 ,y-1,r)))+((SMQ/2.0)*(BG_raw(2   ,y+1,r)+BG_raw(2   ,y-1,r)))
      BG(nx,y,r)  = ((1.0-SMP-SMQ)*BG_raw(nx,y,r))+((SMP/3.0)*(BG_raw(nx,y+1,r)+BG_raw(nx,y-1,r)+BG_raw(nx-1,y,r)))+((SMQ/2.0)*(BG_raw(nx-1,y+1,r)+BG_raw(nx-1,y-1,r)))

      BGlat(1,y,r)   = ((1.0-SMP-SMQ)*BGlat_raw(1 ,y,r))+((SMP/3.0)*(BGlat_raw(1 ,y+1,r)+BGlat_raw(2 ,y  ,r)+BGlat_raw(1 ,y-1,r)))+((SMQ/2.0)*(BGlat_raw(2   ,y+1,r)+BGlat_raw(2   ,y-1,r)))
      BGlat(nx,y,r)  = ((1.0-SMP-SMQ)*BGlat_raw(nx,y,r))+((SMP/3.0)*(BGlat_raw(nx,y+1,r)+BGlat_raw(nx,y-1,r)+BGlat_raw(nx-1,y,r)))+((SMQ/2.0)*(BGlat_raw(nx-1,y+1,r)+BGlat_raw(nx-1,y-1,r)))

      BGlon(1,y,r)   = ((1.0-SMP-SMQ)*BGlon_raw(1 ,y,r))+((SMP/3.0)*(BGlon_raw(1 ,y+1,r)+BGlon_raw(2 ,y  ,r)+BGlon_raw(1 ,y-1,r)))+((SMQ/2.0)*(BGlon_raw(2   ,y+1,r)+BGlon_raw(2   ,y-1,r)))
      BGlon(nx,y,r)  = ((1.0-SMP-SMQ)*BGlon_raw(nx,y,r))+((SMP/3.0)*(BGlon_raw(nx,y+1,r)+BGlon_raw(nx,y-1,r)+BGlon_raw(nx-1,y,r)))+((SMQ/2.0)*(BGlon_raw(nx-1,y+1,r)+BGlon_raw(nx-1,y-1,r)))
    end do ; end do
  aslope_raw = aslope
  BG_raw     = BG
  BGlat_raw  = BGlat
  BGlon_raw  = BGlon

end do
end if

if(trim(cyclic) == 'yes') then						!9 point smoothing at the east and west boundaries by wrapping the cyclic point
do i=1,0!NSM
  do r=1,nr
    do y=2,ny-1								
      if (abs(lat(y)).lt.20.or.abs(lat(y)).gt.70) cycle			! Don't compute values poleward of 70 degrees or between 20S and 20N
      do x=2,nx-1							! Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
        aslope(x,y,r) = ((1.0-SMP-SMQ)*aslope_raw(x,y,r))+((SMP/4.0)*(aslope_raw(x,y+1,r)+aslope_raw(x+1,y,r)+aslope_raw(x,y-1,r)+aslope_raw(x-1,y,r)))+((SMQ/4.0)*(aslope_raw(x-1,y+1,r)+aslope_raw(x+1,y+1,r)+aslope_raw(x+1,y-1,r)+aslope_raw(x-1,y-1,r)))
        
	BG(x,y,r) = ((1.0-SMP-SMQ)*BG_raw(x,y,r))+((SMP/4.0)*(BG_raw(x,y+1,r)+BG_raw(x+1,y,r)+BG_raw(x,y-1,r)+BG_raw(x-1,y,r)))+((SMQ/4.0)*(BG_raw(x-1,y+1,r)+BG_raw(x+1,y+1,r)+BG_raw(x+1,y-1,r)+BG_raw(x-1,y-1,r)))
        
	BGlat(x,y,r) = ((1.0-SMP-SMQ)*BGlat_raw(x,y,r))+((SMP/4.0)*(BGlat_raw(x,y+1,r)+BGlat_raw(x+1,y,r)+BGlat_raw(x,y-1,r)+BGlat_raw(x-1,y,r)))+((SMQ/4.0)*(BGlat_raw(x-1,y+1,r)+BGlat_raw(x+1,y+1,r)+BGlat_raw(x+1,y-1,r)+BGlat_raw(x-1,y-1,r)))

	BGlon(x,y,r) = ((1.0-SMP-SMQ)*BGlon_raw(x,y,r))+((SMP/4.0)*(BGlon_raw(x,y+1,r)+BGlon_raw(x+1,y,r)+BGlon_raw(x,y-1,r)+BGlon_raw(x-1,y,r)))+((SMQ/4.0)*(BGlon_raw(x-1,y+1,r)+BGlon_raw(x+1,y+1,r)+BGlon_raw(x+1,y-1,r)+BGlon_raw(x-1,y-1,r)))
      end do
      
      aslope(1,y,r)   = ((1.0-SMP-SMQ)*aslope_raw(1 ,y,r))+((SMP/4.0)*(aslope_raw(1 ,y+1,r)+aslope_raw(2,y,r)+aslope_raw(1 ,y-1,r)+aslope_raw(nx  ,y,r)))+((SMQ/4.0)*(aslope_raw(nx  ,y+1,r)+aslope_raw(2,y+1,r)+aslope_raw(2,y-1,r)+aslope_raw(nx  ,y-1,r)))
      aslope(nx,y,r)  = ((1.0-SMP-SMQ)*aslope_raw(nx,y,r))+((SMP/4.0)*(aslope_raw(nx,y+1,r)+aslope_raw(1,y,r)+aslope_raw(nx,y-1,r)+aslope_raw(nx-1,y,r)))+((SMQ/4.0)*(aslope_raw(nx-1,y+1,r)+aslope_raw(1,y+1,r)+aslope_raw(1,y-1,r)+aslope_raw(nx-1,y-1,r)))
      
      BG(1,y,r)   = ((1.0-SMP-SMQ)*BG_raw(1 ,y,r))+((SMP/4.0)*(BG_raw(1 ,y+1,r)+BG_raw(2,y,r)+BG_raw(1 ,y-1,r)+BG_raw(nx  ,y,r)))+((SMQ/4.0)*(BG_raw(nx  ,y+1,r)+BG_raw(2,y+1,r)+BG_raw(2,y-1,r)+BG_raw(nx  ,y-1,r)))
      BG(nx,y,r)  = ((1.0-SMP-SMQ)*BG_raw(nx,y,r))+((SMP/4.0)*(BG_raw(nx,y+1,r)+BG_raw(1,y,r)+BG_raw(nx,y-1,r)+BG_raw(nx-1,y,r)))+((SMQ/4.0)*(BG_raw(nx-1,y+1,r)+BG_raw(1,y+1,r)+BG_raw(1,y-1,r)+BG_raw(nx-1,y-1,r)))
    
      BGlat(1,y,r)   = ((1.0-SMP-SMQ)*BGlat_raw(1 ,y,r))+((SMP/4.0)*(BGlat_raw(1 ,y+1,r)+BGlat_raw(2,y,r)+BGlat_raw(1 ,y-1,r)+BGlat_raw(nx  ,y,r)))+((SMQ/4.0)*(BGlat_raw(nx  ,y+1,r)+BGlat_raw(2,y+1,r)+BGlat_raw(2,y-1,r)+BGlat_raw(nx  ,y-1,r)))
      BGlat(nx,y,r)  = ((1.0-SMP-SMQ)*BGlat_raw(nx,y,r))+((SMP/4.0)*(BGlat_raw(nx,y+1,r)+BGlat_raw(1,y,r)+BGlat_raw(nx,y-1,r)+BGlat_raw(nx-1,y,r)))+((SMQ/4.0)*(BGlat_raw(nx-1,y+1,r)+BGlat_raw(1,y+1,r)+BGlat_raw(1,y-1,r)+BGlat_raw(nx-1,y-1,r)))

      BGlon(1,y,r)   = ((1.0-SMP-SMQ)*BGlon_raw(1 ,y,r))+((SMP/4.0)*(BGlon_raw(1 ,y+1,r)+BGlon_raw(2,y,r)+BGlon_raw(1 ,y-1,r)+BGlon_raw(nx  ,y,r)))+((SMQ/4.0)*(BGlon_raw(nx  ,y+1,r)+BGlon_raw(2,y+1,r)+BGlon_raw(2,y-1,r)+BGlon_raw(nx  ,y-1,r)))
      BGlon(nx,y,r)  = ((1.0-SMP-SMQ)*BGlon_raw(nx,y,r))+((SMP/4.0)*(BGlon_raw(nx,y+1,r)+BGlon_raw(1,y,r)+BGlon_raw(nx,y-1,r)+BGlon_raw(nx-1,y,r)))+((SMQ/4.0)*(BGlon_raw(nx-1,y+1,r)+BGlon_raw(1,y+1,r)+BGlon_raw(1,y-1,r)+BGlon_raw(nx-1,y-1,r)))
    end do ; end do
  aslope_raw = aslope
  BG_raw     = BG
  BGlat_raw  = BGlat
  BGlon_raw  = BGlon
  end do
end if
call cpu_time(FN)
print '("Time to Smooth Slopes = ",f6.3," seconds.")',FN-ST
end if






if(trim(smooth) == "smth25") then
call cpu_time(ST)
!!!!!!Smooth height field!!!!!!!!!!!9-point smoothing

aslope_raw = aslope
BG_raw     = BG
BGlat_raw  = BGlat
BGlon_raw  = BGlon

if(trim(cyclic) == 'yes') then						!25 point smoothing at the east and west boundaries by wrapping the cyclic point
do i=1,1!NSM
  do r=1,nr
  do y=3,ny-2								!Values at 90 and -90 are constant
    do x=3,nx-2								!Weighted avereage of CENTER(W=(1-P-Q)), UP/RIGHT/DOWN/LEFT(W=P), and UL/UR/LR/LL(W=Q)
      aslope(x,y,r) = sum(aslope_raw(x-2:x+2,y-2:y+2,r))/25.0		      !((1.0-SMP-SMQ)*hgt_raw(x,y))+((SMP/4.0)*(hgt_raw(x,y+1)+hgt_raw(x+1,y)+hgt_raw(x,y-1)+hgt_raw(x-1,y)))+((SMQ/4.0)*(hgt_raw(x-1,y+1)+hgt_raw(x+1,y+1)+hgt_raw(x+1,y-1)+hgt_raw(x-1,y-1)))
      BG(x,y,r) = sum(BG_raw(x-2:x+2,y-2:y+2,r))/25.0		      
      BGlat(x,y,r) = sum(BGlat_raw(x-2:x+2,y-2:y+2,r))/25.0	      
      BGlon(x,y,r) = sum(BGlon_raw(x-2:x+2,y-2:y+2,r))/25.0	      
    end do
    
    aslope(2,y,r)    = (sum(aslope_raw(1:4,y-2:y+2,r))+sum(aslope_raw(nx,y-2:y+2,r)))/25.0
    BG(1,y,r)    = (sum(BG_raw(1:3,y-2:y+2,r))+sum(BG_raw(nx-1:nx,y-2:y+2,r)))/25.0
    BGlat(nx,y,r)   = (sum(BGlat_raw(1:2,y-2:y+2,r))+sum(BGlat_raw(nx-2:nx,y-2:y+2,r)))/25.0
    BGlon(nx-1,y,r) = (sum(BGlon_raw(1,y-2:y+2,r))+sum(BGlon_raw(nx-2:nx,y-2:y+2,r)))/25.0

  end do; end do
  aslope_raw = aslope
  BG_raw     = BG
  BGlat_raw  = BGlat
  BGlon_raw  = BGlon
end do
end if
call cpu_time(FN)
print '("Time to Smooth HGT = ",f6.3," seconds.")',FN-ST
end if






if(trim(smooth) == "area") then				! avoid, slow
call cpu_time(ST)
!!!!!!Smooth height field!!!!!!!!!!!Area-averaging
aslope_raw = aslope
BG_raw     = BG
BGlat_raw  = BGlat
BGlon_raw  = BGlon

do y=1,ny								!Average all points within 150 km
  do x=1,nx
    gridmask = 0								
    do j=max(y-6,1),min(y+6,ny)
      do i=x-6,x+6
        if(i.lt.1)then
	  ii = ny-(abs(i))
	else if(i.gt.ny)then
	  ii = i-ny
	else 
	  ii = i
	end if
        
	!dist	= haversine(lat(y),lon(x),lat(j),lon(ii))
	!if(dist.le.150.0)then
	  gridmask(ii,j) = 1
	!else
	!  gridmask(i,j) = 0
	!end if
      end do
    end do
    do r=1,nr
      aslope(x,y,r) = sum(aslope_raw(:,:,r),MASK=gridmask.eq.1)/sum(gridmask)
      BG(x,y,r)     = sum(BG_raw(:,:,r)    ,MASK=gridmask.eq.1)/sum(gridmask)
      BGlat(x,y,r)  = sum(BGlat_raw(:,:,r) ,MASK=gridmask.eq.1)/sum(gridmask)
      BGlon(x,y,r)  = sum(BGlon_raw(:,:,r) ,MASK=gridmask.eq.1)/sum(gridmask)
    end do
  end do
end do
call cpu_time(FN)
print '("Time to Smooth Slope = ",f6.3," seconds.")',FN-ST
end if







call cpu_time(ST)   
!!!!!!Maximize Slope Function!!!!!!
aslopex = maxval(aslope,dim=3)						! Maximize the slope function over the radius dimension




!!!!!!Record optimal radius!!!!!!
aslopexri = maxloc(aslope,dim=3)					! Find the index of the optimal radius and record the optimal radius, BG slope at the optimal radius
do y=1,ny
  do x=1,nx
    aslopexr(x,y) = radius(aslopexri(x,y))
    BGx(x,y) = BG(x,y,aslopexri(x,y))
    BGxlon(x,y) = BGlon(x,y,aslopexri(x,y))
    BGxlat(x,y) = BGlat(x,y,aslopexri(x,y))
end do ; end do



!!!!!!Record Slope ratio!!!!!
srat = BGx/aslopex							! Compute the slope ratio (Background/Optimal Slope)



!!!!!!Determine if So is local max & hgt is local min!!!!!
do y=2,ny-1								! Find local MAXIMA of the optimal slope function and MINIMA of height 
  do x=2,nx-1
    ! Find over "normal" points (no cyclic wrapping)
    if((aslopex(x,y).gt.aslopex(x-1,y+1)).and.(aslopex(x,y).gt.aslopex(x,y+1)).and.(aslopex(x,y).gt.aslopex(x+1,y+1)).and.(aslopex(x,y).gt.aslopex(x+1,y)).and.(aslopex(x,y).gt.aslopex(x+1,y-1)).and.(aslopex(x,y).gt.aslopex(x,y-1)).and.(aslopex(x,y).gt.aslopex(x-1,y-1)).and.(aslopex(x,y).gt.aslopex(x-1,y))) then
      asmask(x,y) = 1    
    end if
					
    if((hgt(x,y).lt.hgt(x-1,y+1)).and.(hgt(x,y).lt.hgt(x,y+1)).and.(hgt(x,y).lt.hgt(x+1,y+1)).and.(hgt(x,y).lt.hgt(x+1,y)).and.(hgt(x,y).lt.hgt(x+1,y-1)).and.(hgt(x,y).lt.hgt(x,y-1)).and.(hgt(x,y).lt.hgt(x-1,y-1)).and.(hgt(x,y).lt.hgt(x-1,y))) then
      zmask(x,y) = 1
    end if
  
  end do
  
  ! Wrap cyclic point from x = 1
  x = 1
  xm1 = nx
  if((aslopex(x,y).gt.aslopex(xm1,y+1)).and.(aslopex(x,y).gt.aslopex(x,y+1)).and.(aslopex(x,y).gt.aslopex(x+1,y+1)).and.(aslopex(x,y).gt.aslopex(x+1,y)).and.(aslopex(x,y).gt.aslopex(x+1,y-1)).and.(aslopex(x,y).gt.aslopex(x,y-1)).and.(aslopex(x,y).gt.aslopex(xm1,y-1)).and.(aslopex(x,y).gt.aslopex(xm1,y))) then
    asmask(x,y) = 1    
  end if
					
  if((hgt(x,y).lt.hgt(xm1,y+1)).and.(hgt(x,y).lt.hgt(x,y+1)).and.(hgt(x,y).lt.hgt(x+1,y+1)).and.(hgt(x,y).lt.hgt(x+1,y)).and.(hgt(x,y).lt.hgt(x+1,y-1)).and.(hgt(x,y).lt.hgt(x,y-1)).and.(hgt(x,y).lt.hgt(xm1,y-1)).and.(hgt(x,y).lt.hgt(xm1,y))) then
    zmask(x,y) = 1
  end if
  
  ! Wrap cyclic point from x = nx
  x = nx
  xp1 = 1
  if((aslopex(x,y).gt.aslopex(x-1,y+1)).and.(aslopex(x,y).gt.aslopex(x,y+1)).and.(aslopex(x,y).gt.aslopex(xp1,y+1)).and.(aslopex(x,y).gt.aslopex(xp1,y)).and.(aslopex(x,y).gt.aslopex(xp1,y-1)).and.(aslopex(x,y).gt.aslopex(x,y-1)).and.(aslopex(x,y).gt.aslopex(x-1,y-1)).and.(aslopex(x,y).gt.aslopex(x-1,y))) then
    asmask(x,y) = 1    
  end if
					
  if((hgt(x,y).lt.hgt(x-1,y+1)).and.(hgt(x,y).lt.hgt(x,y+1)).and.(hgt(x,y).lt.hgt(xp1,y+1)).and.(hgt(x,y).lt.hgt(xp1,y)).and.(hgt(x,y).lt.hgt(xp1,y-1)).and.(hgt(x,y).lt.hgt(x,y-1)).and.(hgt(x,y).lt.hgt(x-1,y-1)).and.(hgt(x,y).lt.hgt(x-1,y))) then
    zmask(x,y) = 1
  end if    
end do 

call cpu_time(FN)
print '("Time to compute So, Ro, SR, min(HGT)= ",f6.3," seconds.")',FN-ST




call cpu_time(ST)
!!!!!!Allocate arrays to hold features (should be more efficient than looking at every GP)!!!!!!
np1 = int(sum(asmask))
nz1 = int(sum(zmask))
allocate(features_p1(np1,10),allowedR_1D(np1))
allocate(features_z(nz1,3))


!!!!!!Assign feature data to arrays!!!!!!
allowedR_1D	 = pack(allowedR,asmask.eq.1)
features_p1(:,1) = pack(aslopex,asmask.eq.1)				! Take advantage of fortran pack()...Analagous to a combination of NCL's where() and ndtooned() procedures. Make 1D arrays of features where local maxima are identified and save to first pass (p1) list variable
features_p1(:,2) = pack(lat2d,asmask.eq.1)				! @AS lat				
features_p1(:,3) = pack(lon2d,asmask.eq.1)				! @AS lon
features_p1(:,4) = pack(aslopexr,asmask.eq.1)				! @AS Ro
features_p1(:,5) = pack(aslopexri,asmask.eq.1)				! @AS Ro index
features_p1(:,6) = pack(srat,asmask.eq.1)				! @AS SRAT
features_p1(:,7) = pack(BGx,asmask.eq.1)				! @AS BG
features_p1(:,8) = pack(BGxlat,asmask.eq.1)				! @AS BG y-comp
features_p1(:,9) = pack(BGxlon,asmask.eq.1)				! @AS BG x-comp
features_p1(:,10) = 1							! @AS flag (for next step)
!features_p1(:,11) = pack(z850,asmask.eq.1)
!features_p1(:,12) = pack(hgt_raw,asmask.eq.1)
!features_p1(:,13) = pack(z200,asmask.eq.1)
!features_p1(:,14) = pack(t850,asmask.eq.1)
!features_p1(:,15) = pack(t500,asmask.eq.1)
!features_p1(:,16) = pack(t200,asmask.eq.1)
!features_p1(:,17) = pack(u850,asmask.eq.1)
!features_p1(:,18) = pack(u500,asmask.eq.1)
!features_p1(:,19) = pack(u200,asmask.eq.1)
!features_p1(:,20) = pack(v850,asmask.eq.1)
!features_p1(:,21) = pack(v500,asmask.eq.1)
!features_p1(:,22) = pack(v200,asmask.eq.1)
!features_p1(:,23) = pack(rh850,asmask.eq.1)
!features_p1(:,24) = pack(rh500,asmask.eq.1)
!features_p1(:,25) = pack(rh200,asmask.eq.1)


b = 1									! Start a counter
									! Check each feature. reset flag to 0 based on the below criteria
do a = 1, np1
  !if(abs(features_p1(a,2)).le.20.or.abs(features_p1(a,2)).ge.70) then	! Fail test - out of domain
  !  features_p1(a,10) = 0
  !  cycle
  !end if
  if(abs(features_p1(a,2)).eq.90) then					! Fail test - pole
    features_p1(a,10) = 0
    cycle
  end if
  if(features_p1(a,1).lt.somin) then					! Fail test - So is less that or equal to zero (or a larger user defined threshold value)
    features_p1(a,10) = 0
    cycle
  end if
  if(features_p1(a,5).eq.1.or.features_p1(a,5).eq.allowedR_1D(a)) then	! Fail test - Ro is rmin or rmax or max allowable rmax near poles
  !if(features_p1(a,5).eq.1.or.features_p1(a,5).eq.nr) then		! Fail test - Ro is rmin or rmax
    features_p1(a,10) = 0
    cycle
  end if
  if(features_p1(a,6).ge.srmax) then					! Fail test - SR is too large
    features_p1(a,10) = 0	
    cycle
  end if
  b=b+1									! All tests passed - Incremement the counter
end do

!!!!!!Allocate arrays to hold second pass (distance threshold) features (should be more efficient than looking at every GP)!!!!!!
np2 = b-1
allocate(features_p2(np2,10))

features_p2(:,1) = pack(features_p1(:,1),features_p1(:,10).eq.1)
features_p2(:,2) = pack(features_p1(:,2),features_p1(:,10).eq.1)
features_p2(:,3) = pack(features_p1(:,3),features_p1(:,10).eq.1)
features_p2(:,4) = pack(features_p1(:,4),features_p1(:,10).eq.1)
features_p2(:,5) = pack(features_p1(:,5),features_p1(:,10).eq.1)
features_p2(:,6) = pack(features_p1(:,6),features_p1(:,10).eq.1)
features_p2(:,7) = pack(features_p1(:,7),features_p1(:,10).eq.1)
features_p2(:,8) = pack(features_p1(:,8),features_p1(:,10).eq.1)
features_p2(:,9) = pack(features_p1(:,9),features_p1(:,10).eq.1)
features_p2(:,10) = 1
!features_p2(:,11) = pack(features_p1(:,11),features_p1(:,10).eq.1)
!features_p2(:,12) = pack(features_p1(:,12),features_p1(:,10).eq.1)
!features_p2(:,13) = pack(features_p1(:,13),features_p1(:,10).eq.1)
!features_p2(:,14) = pack(features_p1(:,14),features_p1(:,10).eq.1)
!features_p2(:,15) = pack(features_p1(:,15),features_p1(:,10).eq.1)
!features_p2(:,16) = pack(features_p1(:,16),features_p1(:,10).eq.1)
!features_p2(:,17) = pack(features_p1(:,17),features_p1(:,10).eq.1)
!features_p2(:,18) = pack(features_p1(:,18),features_p1(:,10).eq.1)
!features_p2(:,19) = pack(features_p1(:,19),features_p1(:,10).eq.1)
!features_p2(:,20) = pack(features_p1(:,20),features_p1(:,10).eq.1)
!features_p2(:,21) = pack(features_p1(:,21),features_p1(:,10).eq.1)
!features_p2(:,22) = pack(features_p1(:,22),features_p1(:,10).eq.1)
!features_p2(:,23) = pack(features_p1(:,23),features_p1(:,10).eq.1)
!features_p2(:,24) = pack(features_p1(:,24),features_p1(:,10).eq.1)
!features_p2(:,25) = pack(features_p1(:,25),features_p1(:,10).eq.1)


if(0.gt.1)then		! Skip this block...this is the original "if features is within Ro of another feature"
do a=1, np2							       ! For features that are too close together, pick the steeper slope
  if(features_p2(a,10).eq.0)cycle				       ! Skip features that have already been de-flagged
  do i=1, np2
    if(i.eq.a) cycle
    if(features_p2(i,10).eq.0)cycle				       ! Skip features that have already been de-flagged
    d = haversine(features_p2(a,2),features_p2(a,3),features_p2(i,2),features_p2(i,3))*1000.0
    if((d.lt.features_p2(a,4).or.d.lt.features_p2(i,4)).and.features_p2(i,1).lt.features_p2(a,1)) then
      features_p2(i,10) = 0
      !exit ; cycle
    end if
  end do
end do


do a=1, np2							       ! For features that are close together, eliminate the smaller Ro if a feature with a larger Ro is within it
  if(features_p2(a,10).eq.0)cycle				       ! Skip features that have already been de-flagged
  do i=1, np2
    if(i.eq.a) cycle
    if(features_p2(i,10).eq.0)cycle				       ! Skip features that have already been de-flagged
    x = minloc((/features_p2(a,4),features_p2(i,4)/),1)		       ! Select the feature with smaller radius for comparison to d
    if(x.eq.1) y = a
    if(x.eq.2) y = i
    d = haversine(features_p2(a,2),features_p2(a,3),features_p2(i,2),features_p2(i,3))*1000.0    
    if(d.lt.features_p2(y,4)) then 					
      features_p2(y,10) = 0
      !exit ; cycle
    end if
  end do
end do
end if


!!!!!!Allocate arrays to hold third pass (local height minima) features (should be more efficient than looking at every GP)!!!!!!
!First need to gather local minima
features_z(:,1) = pack(hgt,zmask.eq.1)
features_z(:,2) = pack(lat2d,zmask.eq.1)
features_z(:,3) = pack(lon2d,zmask.eq.1)

np3 = int(sum(features_p2(:,10)))
allocate(features_p3(np3,12))

features_p3(:,1) = pack(features_p2(:,1),features_p2(:,10).eq.1)
features_p3(:,2) = pack(features_p2(:,2),features_p2(:,10).eq.1)
features_p3(:,3) = pack(features_p2(:,3),features_p2(:,10).eq.1)
features_p3(:,4) = pack(features_p2(:,4),features_p2(:,10).eq.1)
features_p3(:,5) = pack(features_p2(:,6),features_p2(:,10).eq.1)
features_p3(:,6) = pack(features_p2(:,7),features_p2(:,10).eq.1)
features_p3(:,7) = pack(features_p2(:,8),features_p2(:,10).eq.1)
features_p3(:,8) = pack(features_p2(:,9),features_p2(:,10).eq.1)
features_p3(:,9:11) = -9999.9
features_p3(:,12)= 1
!features_p3(:,13) = pack(features_p2(:,11),features_p2(:,10).eq.1)
!features_p3(:,14) = pack(features_p2(:,12),features_p2(:,10).eq.1)
!features_p3(:,15) = pack(features_p2(:,13),features_p2(:,10).eq.1)
!features_p3(:,16) = pack(features_p2(:,14),features_p2(:,10).eq.1)
!features_p3(:,17) = pack(features_p2(:,15),features_p2(:,10).eq.1)
!features_p3(:,18) = pack(features_p2(:,16),features_p2(:,10).eq.1)
!features_p3(:,19) = pack(features_p2(:,17),features_p2(:,10).eq.1)
!features_p3(:,20) = pack(features_p2(:,18),features_p2(:,10).eq.1)
!features_p3(:,21) = pack(features_p2(:,19),features_p2(:,10).eq.1)
!features_p3(:,22) = pack(features_p2(:,20),features_p2(:,10).eq.1)
!features_p3(:,23) = pack(features_p2(:,21),features_p2(:,10).eq.1)
!features_p3(:,24) = pack(features_p2(:,22),features_p2(:,10).eq.1)
!features_p3(:,25) = pack(features_p2(:,23),features_p2(:,10).eq.1)
!features_p3(:,26) = pack(features_p2(:,24),features_p2(:,10).eq.1)
!features_p3(:,27) = pack(features_p2(:,25),features_p2(:,10).eq.1)

do a=1, np3
  hgtmin = 100000.0
  do i=1, nz1
    d = haversine(features_p3(a,2),features_p3(a,3),features_z(i,2),features_z(i,3))*1000.0
    if(d.lt.features_p3(a,4).and.features_z(i,1).lt.hgtmin) then
       features_p3(a,9:11) = features_z(i,:)
       hgtmin = features_z(i,1)
    end if
  end do
end do

do a=1,np3
  i = findloc(lon,features_p3(a,3),1)
  j = findloc(lat,features_p3(a,2),1)  
  asmask2(i,j) = 1
end do


! Check for shared height minima! Check for shared height minima
do a=1, np3
  if(features_p3(a,9).eq.-9999.9)cycle
  if(features_p3(a,12).eq.0)cycle
  do i=1, np3
    if(i.eq.a) cycle
    if(features_p3(i,9).eq.-9999.9)cycle
    if(features_p3(i,12).eq.0)cycle
    x = minloc((/features_p3(a,1),features_p3(i,1)/),1)		! Select the feature with smaller So for comparison to d
    !x = minloc((/features_p3(a,4),features_p3(i,4)/),1)			! Select the feature with smaller Ro for comparison to d
    if(x.eq.1) y = a
    if(x.eq.2) y = i
    
    d = haversine(features_p3(a,10),features_p3(a,11),features_p3(i,10),features_p3(i,11))
    if(d.lt.10.0) then
       features_p3(y,12) = 0						! Fail feature (a or i) if it shares a Zmin with a stronger So
       !features_p3(y,12) = 0						! Fail feature (a or i) if it shares a Zmin with a Larger Ro
    end if
  end do
end do


if(0.gt.1) then

! Check for height minima within optimal radius of a different cutoff
do a=1, np3
  if(features_p3(a,9).eq.-9999.9)cycle
  do i=1, np3
    if(i.eq.a) cycle
    if(features_p3(i,9).eq.-9999.9)cycle
    
    x = minloc((/features_p3(a,1),features_p3(i,1)/),1)		       ! Select the feature with larger So 
    if(x.eq.1) y = a
    if(x.eq.2) y = i
    
    d = haversine(features_p3(a,2),features_p3(a,3),features_p3(i,10),features_p3(i,11))*1000.0
    if(d.lt.features_p3(a,4)) then
       features_p3(y,12) = 0						! Fail smaller Ro feature (a or i) if i's Zmin is within Ro of a 
    end if
  end do
end do
end if


! Check for engulfed features within larger features
do a=1, np3
  !if(features_p3(a,9).ne.-9999.9)cycle
  if(features_p3(a,12).eq.0)cycle
  do i=1, np3
    if(i.eq.a) cycle
    if(features_p3(i,12).eq.0)cycle
    !if(features_p3(i,9).ne.-9999.9)cycle
    x = minloc((/features_p3(a,1),features_p3(i,1)/),1)		       ! Select the feature with smaller So for comparison to d
   ! x = minloc((/features_p3(a,4),features_p3(i,4)/),1)		       ! Select the feature with smaller Ro for comparison to d
    if(x.eq.1) y = a
    if(x.eq.2) y = i
    
    d = haversine(features_p3(a,2),features_p3(a,3),features_p3(i,2),features_p3(i,3))*1000.0 	! to meters
    o = c1c2_overlap(d,features_p3(a,4),features_p3(i,4))
    if(o.ge.0.9) then
       features_p3(y,12) = 0						! Fail feature (a or i) if it has the smaller So of two 90% overlapping features
    end if
  end do
end do

if(0.gt.1)then
! Check for engulfed troughs within larger cutoffs
do a=1, np3
  if(features_p3(a,9).eq.-9999.9)cycle
  if(features_p3(a,12).eq.0)cycle
  do i=1, np3
    if(i.eq.a) cycle
    if(features_p3(i,12).eq.0)cycle
    if(features_p3(i,9).ne.-9999.9)cycle
    x = minloc((/features_p3(a,4),features_p3(i,4)/),1)		       ! Select the feature with smaller So for comparison to d
    if(x.eq.1) cycle						       ! Omit the case if the cutoff is smaller
    if(x.eq.2) y = i
    
    d = haversine(features_p3(a,2),features_p3(a,3),features_p3(i,2),features_p3(i,3))*1000.0 	! to meters
    o = c1c2_overlap(d,features_p3(a,4),features_p3(i,4))
    if(o.gt.0.80) then
       features_p3(y,12) = 0						! Fail the trough if it is smaller than the cutoff and 80% engulfed by the cutoff
    end if
  end do
end do
end if

!!!!!!Allocate arrays to hold final features (should be more efficient than looking at every GP)!!!!!!
np = int(sum(features_p3(:,12)))
allocate(features(np,32))

features(:,1) = pack(features_p3(:,1),features_p3(:,12).eq.1)
features(:,2) = pack(features_p3(:,2),features_p3(:,12).eq.1)
features(:,3) = pack(features_p3(:,3),features_p3(:,12).eq.1)
features(:,4) = pack(features_p3(:,4),features_p3(:,12).eq.1)
features(:,5) = pack(features_p3(:,5),features_p3(:,12).eq.1)
features(:,6) = pack(features_p3(:,6),features_p3(:,12).eq.1)		 
features(:,7) = pack(features_p3(:,7),features_p3(:,12).eq.1)
features(:,8) = pack(features_p3(:,8),features_p3(:,12).eq.1)
features(:,9) = pack(features_p3(:,9),features_p3(:,12).eq.1)
features(:,10) = pack(features_p3(:,10),features_p3(:,12).eq.1)
features(:,11) = pack(features_p3(:,11),features_p3(:,12).eq.1)

call cpu_time(STsub)
!!!!!!Assign area-average feature data to arrays!!!!!
gridmask_all = 0
do a = 1, np
  gridmask = 0
  i = findloc(lon,features(a,3),1)
  j = findloc(lat,features(a,2),1)    
      								! Determine how many gridpoints to move to get to the specified radius 
  distj = haversine(lat(j),lon(i),lat(j-1),lon(i))*1000.0 	! Distance between gridpoints in the j/y direction. convert to meters
  gpdistj = features(a,4)/distj 				! Number of gridpoints in j/y to reach radius [m / (m/point) = points]
  j1 = ceiling(gpdistj)
  
  do y = j-j1, j+j1
    do x = 1, nx
      dist = haversine(features(a,2),features(a,3),lat(y),lon(x))*1000.0
      if(dist.le.features(a,4)) then
        gridmask(x,y) = 1
	gridmask_all(x,y) = 1
      endif
    end do
  end do
  
  features(a,12) = sum(z850,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,13) = sum(hgt_raw,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,14) = sum(z200,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,15) = sum(t850,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,16) = sum(t500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,17) = sum(t200,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,18) = sum(u850,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,19) = sum(u500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,20) = sum(u200,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,21) = sum(v850,MASK=gridmask.eq.1)/sum(gridmask) 
  features(a,22) = sum(v500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,23) = sum(v200,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,24) = sum(mr850,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,25) = sum(mr500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,26) = sum(mr200,MASK=gridmask.eq.1)/sum(gridmask)
end do
call cpu_time(FNsub)  
print '("Subtime to do feature Ro area averages= ",f6.3," seconds.")',FNsub-STsub


call cpu_time(STsub)
!!!!!!Assign FIXED RADIUS area-average feature data to arrays!!!!!
gridmask_all = 0
do a = 1, np
  gridmask = 0
  i = findloc(lon,features(a,3),1)
  j = findloc(lat,features(a,2),1)    
      								! Determine how many gridpoints to move to get to the specified radius 
  distj = haversine(lat(j),lon(i),lat(j-1),lon(i))*1000.0 	! Distance between gridpoints in the j/y direction. convert to meters
  gpdistj = 600000.0/distj 					! Number of gridpoints in j/y to reach 600 km [m / (m/point) = points]
  j1 = ceiling(gpdistj)
  
  do y = j-j1, j+j1
    do x = 1, nx
      dist = haversine(features(a,2),features(a,3),lat(y),lon(x))*1000.0
      if(dist.le.600000.0) then
        gridmask(x,y) = 1
	gridmask_all(x,y) = 1
      endif
    end do
  end do
  
  features(a,27) = sum(hgt_raw,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,28) = sum(t500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,29) = sum(u500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,30) = sum(v500,MASK=gridmask.eq.1)/sum(gridmask)
  features(a,31) = sum(mr500,MASK=gridmask.eq.1)/sum(gridmask)
end do
call cpu_time(FNsub)  
print '("Subtime to do feature 600 km area averages= ",f6.3," seconds.")',FNsub-STsub



call cpu_time(STsub)
!!!!!!Assign Trough, Closed, and Cutoff Tags to features (adapted from Munoz et al 2020)!!!!!
gridmask_all = 0
do a = 1, np
  
  if(features(a,2).gt.70.0.or.features(a,2).lt.-70.0.or.abs(features(a,2)).lt.20.0) then
    features(a,32) = 0	! Unclassified (polar or equator)
    cycle
  endif
  if(features(a,9).eq.-9999.9) then
    features(a,32) = 1	! Trough (no z min)
    cycle
  endif  
  features(a,32)   = 2	! Closed low (not unclassified and has a zmin)
  
  gridmask = 0
  i0 = findloc(lon,features(a,3),1)
  j0 = findloc(lat,features(a,2),1)    
  
  nrow = 0
  do j = 6, 35
    x = i0
    if(features(a,2).ge.20) then
      y = j0+j ! NHEM poleward
    end if
    if(features(a,2).le.-20 )then
      y = j0-j ! SHEM poleward
    end if
    
    gridmask(x-5:x+5,y-5:y+5) = 1
    Uavg = sum(u500,MASK=gridmask.eq.1)/sum(gridmask)
    
    if(Uavg.lt.-10) then
      nrow = nrow + 1
    else
      nrow = 0
    end if
    
    if(nrow.eq.10)then
      features(a,32) = 3    ! Cutoff low (closed low that has a poleward wind reversal for 2.5 degrees)
      exit		    ! criteria is met for 10 points
    endif
  
  end do
end do
call cpu_time(FNsub)  
print '("Subtime to add feature type tags= ",f6.3," seconds.")',FNsub-STsub

!features(:,12) = pack(features_p3(:,13),features_p3(:,12).eq.1)
!features(:,13) = pack(features_p3(:,14),features_p3(:,12).eq.1)
!features(:,14) = pack(features_p3(:,15),features_p3(:,12).eq.1)
!features(:,15) = pack(features_p3(:,16),features_p3(:,12).eq.1)
!features(:,16) = pack(features_p3(:,17),features_p3(:,12).eq.1)
!features(:,17) = pack(features_p3(:,18),features_p3(:,12).eq.1)
!features(:,18) = pack(features_p3(:,19),features_p3(:,12).eq.1)
!features(:,19) = pack(features_p3(:,20),features_p3(:,12).eq.1)
!features(:,20) = pack(features_p3(:,21),features_p3(:,12).eq.1)
!features(:,21) = pack(features_p3(:,22),features_p3(:,12).eq.1)
!features(:,22) = pack(features_p3(:,23),features_p3(:,12).eq.1)
!features(:,23) = pack(features_p3(:,24),features_p3(:,12).eq.1)
!features(:,24) = pack(features_p3(:,25),features_p3(:,12).eq.1)
!features(:,25) = pack(features_p3(:,26),features_p3(:,12).eq.1)
!features(:,26) = pack(features_p3(:,27),features_p3(:,12).eq.1)

do a=1,np
  i = findloc(lon,features(a,3),1)
  j = findloc(lat,features(a,2),1)  
  asmask3(i,j) = 1
end do

if(debug.ne."debugonly".and.debug.ne."binaryonly")then
!!!!Write to text file!!!!
  open(UNUM, file=trim(outfile), status="unknown")

  fstr = "(A7, 1x, A7, 1x, A11, 1x, A7, 1x, A7, 1x, A7, 1x, A4, 1x, A15, 1x, A16, 1x, A16, 1x, A7, 1x, A7, 1x, A7, 1x, A22, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A11, 1x, A17, 1x, A17, 1x, A17, 1x, A17, 1x, A17)"
  write(UNUM,fstr) "ITIME","FHOUR","So(m/100km)","LAT(N)","LON(E)","Ro(km)","SR","BGo(m/100km)","BGo-lat(m/100km)","BGo-lon(m/100km)","ZMIN(m)","ZLAT(N)","ZLON(E)","U(0)/T(1)/Cl(2)/Cu(3)","Z850(m)","Z500(m)","Z200(m)","T850(K)","T500(K)","T200(K)","U850(m/s)","U500(m/s)","U200(m/s)","V850(m/s)","V500(m/s)","V200(m/s)","MR850(g/kg)","MR500(g/kg)","MR200(g/kg)","600kmZ500(m)","600kmT500(K)","600kmU500(m/s)","600kmV500(m/s)","600kmMR500(g/kg)"

  fstr = "(A11, 1x, A4, 1x, F6.2, 1x, F6.2, 1x, F6.2, 1x, F10.2, 1x, F6.2, 1x, F6.2, 1x, F6.2, 1x, F6.2, 1x, F8.2, 1x, F8.2, 1x, F8.2, 1x, F3.1, 1x, F8.2, 1x, F8.2, 1x, F8.2, 1x, F5.1,1x,F5.1,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.4,1x,F7.4,1x,F7.4,1x,F8.2,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.4)"
  do b=1,np  
    write(UNUM,fstr) itime, fhour, features(b,1), features(b,2), features(b,3), features(b,4)/1000.0, features(b,5), features(b,6), features(b,7), features(b,8), features(b,9), features(b,10), features(b,11), features(b,32), features(b,12), features(b,13), features(b,14), features(b,15), features(b,16), features(b,17), features(b,18), features(b,19), features(b,20), features(b,21), features(b,22), features(b,23), features(b,24), features(b,25), features(b,26), features(b,27), features(b,28), features(b,29), features(b,30), features(b,31)
  end do 
  close(UNUM)
  call cpu_time(FN)
  print '("Time to make list and output = ",f6.3," seconds.")',FN-ST
end if

  

if(debug.eq."debug".or.debug.eq."debugonly")then ! Enter debug mode
call cpu_time(ST)

!!!!For debugging!!!!    
! Create the output file
    write(*,*) "Creating:  " // trim(doutfile)
    !Define grid dimensions
    call nc_check( nf90_create( trim(doutfile), nf90_clobber, ncIDo ), "kasugaEA21", "open kasugaEA21_slope" )
    call nc_check( nf90_def_dim(ncIDo, "time", nf90_unlimited, t_dimid), "kasugaEA21", "def dim time kasugaEA21_slope" )
    
    call nc_check( nf90_def_dim(ncIDo, "latitude", ny, y_dimid), "kasugaEA21", "def dim south_north kasugaEA21_slope" )

    call nc_check( nf90_def_dim(ncIDo, "radius", nr, r_dimid), "kasugaEA21", "def dim bottom_top kasugaEA21_slope" )
    
    call nc_check( nf90_def_dim(ncIDo, "longitude", nx, x_dimid), "kasugaEA21", "def dim west_east kasugaEA21_slope" )

    dimids = (/ x_dimid, y_dimid, r_dimid, t_dimid /)
    ct = (/ nx, ny, nr, nt /)
    start =  (/ 1, 1, 1, 1 /) 
   
    dimids_nor = (/ x_dimid, y_dimid, t_dimid /)
    ct_nor = (/ nx, ny, nt /)
    start_nor =  (/ 1, 1, 1 /) 

    ct_t = (/ nt /)
    start_t = (/ 1 /)
        
    ct_x = (/ nx /)
    start_x = (/ 1 /)
    
    ct_y = (/ ny /)
    start_y = (/ 1 /)
    
    ct_r = (/ nr /)
    start_r = (/ 1 /)    

    !Initiailize vars
    call nc_check( nf90_def_var(ncIDo, "HGT_500mb", nf90_real, dimids_nor, varid), "kasugaEA21", "def var hgt kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m"), "kasugaEA21", "def units hgt kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "geopotential height"), "kasugaEA21", "def units hgt kasugaEA21_slope" )
    
    call nc_check( nf90_def_var(ncIDo, "ASLOPE", nf90_real, dimids, varid), "kasugaEA21", "def var aslope kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units aslope kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "average slope function"), "kasugaEA21", "def units aslope kasugaEA21_slope" )
    
    call nc_check( nf90_def_var(ncIDo, "BGSLOPE", nf90_real, dimids, varid), "kasugaEA21", "def var bgslope kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units bgslope kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "background slope function"), "kasugaEA21", "def units bgslope kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "So", nf90_real, dimids_nor, varid), "kasugaEA21", "def var So kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units So kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "optimal slope function"), "kasugaEA21", "def units So kasugaEA21_slope" )
    
    call nc_check( nf90_def_var(ncIDo, "BGo", nf90_real, dimids_nor, varid), "kasugaEA21", "def var BGo kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units BGo kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "background slope function at optimal radius"), "kasugaEA21", "def units BGo kasugaEA21_slope" )
    
    call nc_check( nf90_def_var(ncIDo, "BGo-lat", nf90_real, dimids_nor, varid), "kasugaEA21", "def var BGol-at kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units BGo-lat kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "y-comp of background slope function at optimal radius"), "kasugaEA21", "def units BGo-lat kasugaEA21_slope" )

    call nc_check( nf90_def_var(ncIDo, "BGo-lon", nf90_real, dimids_nor, varid), "kasugaEA21", "def var BGol-lon kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m 100km-1"), "kasugaEA21", "def units BGo-lon kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "x-comp of background slope function at optimal radius"), "kasugaEA21", "def units BGo-lon kasugaEA21_slope" )
   
    call nc_check( nf90_def_var(ncIDo, "Ro", nf90_real, dimids_nor, varid), "kasugaEA21", "def var Ro kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m"), "kasugaEA21", "def units Ro kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "optimal radius"), "kasugaEA21", "def units Ro kasugaEA21_slope" )
    
    call nc_check( nf90_def_var(ncIDo, "SR", nf90_real, dimids_nor, varid), "kasugaEA21", "def var Ro kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", ""), "kasugaEA21", "def units Ro kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "slope ratio"), "kasugaEA21", "def units SR kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "ASMASK", nf90_real, dimids_nor, varid), "kasugaEA21", "def var asmask kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", " "), "kasugaEA21", "def units asmask kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "mask for So maxima"), "kasugaEA21", "def units asmask kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "ASMASK2", nf90_real, dimids_nor, varid), "kasugaEA21", "def var asmask2 kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", " "), "kasugaEA21", "def units asmask2 kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "mask for So maxima after distance checks"), "kasugaEA21", "def units asmask2 kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "ASMASK3", nf90_real, dimids_nor, varid), "kasugaEA21", "def var asmask23kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", " "), "kasugaEA21", "def units asmask3 kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "mask for So maxima after shared zmin checks"), "kasugaEA21", "def units asmask3 kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "MASKRo", nf90_real, dimids_nor, varid), "kasugaEA21", "def var MASKRo kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "binary"), "kasugaEA21", "def units MASKRo kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "binary mask for radius of features after all checks"), "kasugaEA21", "def units MASKRo kasugaEA21_slope" )       
    
    call nc_check( nf90_def_var(ncIDo, "time", nf90_real, t_dimid, varid), "kasugaEA21", "def var time kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "seconds since 1970-01-01 00:00:00.0 0:00"), "kasugaEA21", "def units time kasugaEA21_slope" )  

    call nc_check( nf90_def_var(ncIDo, "latitude", nf90_real, y_dimid, varid), "kasugaEA21", "def var latitude kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "degrees_north"), "kasugaEA21", "def units latitude kasugaEA21_slope" )  
    
    call nc_check( nf90_def_var(ncIDo, "radius", nf90_real, r_dimid, varid), "kasugaEA21", "def var radius kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "m"), "kasugaEA21", "def units radius kasugaEA21_slope" )  
    
    call nc_check( nf90_def_var(ncIDo, "longitude", nf90_real, x_dimid, varid), "kasugaEA21", "def var longitude kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "degrees_east"), "kasugaEA21", "def units longitude kasugaEA21_slope" )  
   
    call nc_check( nf90_enddef(ncIDo), "kasugaEA21", "end define mode kasugaEA21_slope" ) 


    !Write Vars


    call nc_check( nf90_inq_varid(ncIDo, "HGT_500mb", varid), "kasugaEA21", "inq_varid hgt")
    call nc_check( nf90_put_var(ncIDo, varid, hgt, start_nor, ct_nor), "kasugaEA21", "put var hgt kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "ASLOPE", varid), "kasugaEA21", "inq_varid aslope")
    call nc_check( nf90_put_var(ncIDo, varid, aslope, start, ct), "kasugaEA21", "put var aslope kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "BGSLOPE", varid), "kasugaEA21", "inq_varid bgslope")
    call nc_check( nf90_put_var(ncIDo, varid, BG, start, ct), "kasugaEA21", "put var bgslope kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "So", varid), "kasugaEA21", "inq_varid So")
    call nc_check( nf90_put_var(ncIDo, varid, aslopex, start_nor, ct_nor), "kasugaEA21", "put var So kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "BGo", varid), "kasugaEA21", "inq_varid BGo")
    call nc_check( nf90_put_var(ncIDo, varid, BGx, start_nor, ct_nor), "kasugaEA21", "put var BGo kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "BGo-lat", varid), "kasugaEA21", "inq_varid BGo-lat")
    call nc_check( nf90_put_var(ncIDo, varid, BGxlat, start_nor, ct_nor), "kasugaEA21", "put var BGo-lat kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "BGo-lon", varid), "kasugaEA21", "inq_varid BGo-lon")
    call nc_check( nf90_put_var(ncIDo, varid, BGxlon, start_nor, ct_nor), "kasugaEA21", "put var BGo-lon kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "Ro", varid), "kasugaEA21", "inq_varid Ro")
    call nc_check( nf90_put_var(ncIDo, varid, aslopexr, start_nor, ct_nor), "kasugaEA21", "put var Ro kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "SR", varid), "kasugaEA21", "inq_varid SR")
    call nc_check( nf90_put_var(ncIDo, varid, srat, start_nor, ct_nor), "kasugaEA21", "put var SR kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "ASMASK", varid), "kasugaEA21", "inq_varid asmask")
    call nc_check( nf90_put_var(ncIDo, varid, asmask, start_nor, ct_nor), "kasugaEA21", "put var asmask kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "ASMASK2", varid), "kasugaEA21", "inq_varid asmask2")
    call nc_check( nf90_put_var(ncIDo, varid, asmask2, start_nor, ct_nor), "kasugaEA21", "put var asmask2 kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "ASMASK3", varid), "kasugaEA21", "inq_varid asmask3")
    call nc_check( nf90_put_var(ncIDo, varid, asmask3, start_nor, ct_nor), "kasugaEA21", "put var asmask3 kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "MASKRo", varid), "kasugaEA21", "inq_varid MASKRo")
    call nc_check( nf90_put_var(ncIDo, varid, gridmask_all, start_nor, ct_nor), "kasugaEA21", "put var MASKRo kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "time", varid), "kasugaEA21", "inq_varid lat")
    call nc_check( nf90_put_var(ncIDo, varid, time, start_t, ct_t), "kasugaEA21", "put var time kasugaEA21_slope" )
    call nc_check( nf90_inq_varid(ncIDo, "latitude", varid), "kasugaEA21", "inq_varid lat")
    call nc_check( nf90_put_var(ncIDo, varid, lat, start_y, ct_y), "kasugaEA21", "put var lat kasugaEA21_slope" )
    call nc_check( nf90_inq_varid(ncIDo, "radius", varid), "kasugaEA21", "inq_varid radius")
    call nc_check( nf90_put_var(ncIDo, varid, radius, start_r, ct_r), "kasugaEA21", "put var radius kasugaEA21_slope" )
    call nc_check( nf90_inq_varid(ncIDo, "longitude", varid), "kasugaEA21", "inq_varid lon")
    call nc_check( nf90_put_var(ncIDo, varid, lon, start_x, ct_x), "kasugaEA21", "put var lon kasugaEA21_slope" )
    
    ! Close the output file
    call nc_check( nf90_close(ncIDo), "kasugaEA21", "close " // trim(doutfile) )
  
  
call cpu_time(FN)
print '("Time for netcdf= ",f6.3," seconds.")',FN-ST
end if ! End debug mode


if(debug.eq."binary")then ! Enter debug mode
call cpu_time(ST)

!!!!For debugging!!!!    
! Create the output file
    write(*,*) "Creating:  " // trim(doutfile)
    !Define grid dimensions
    call nc_check( nf90_create( trim(doutfile), nf90_clobber, ncIDo ), "kasugaEA21", "open kasugaEA21_slope" )
    call nc_check( nf90_def_dim(ncIDo, "time", nf90_unlimited, t_dimid), "kasugaEA21", "def dim time kasugaEA21_slope" )
    
    call nc_check( nf90_def_dim(ncIDo, "latitude", ny, y_dimid), "kasugaEA21", "def dim south_north kasugaEA21_slope" )
    
    call nc_check( nf90_def_dim(ncIDo, "longitude", nx, x_dimid), "kasugaEA21", "def dim west_east kasugaEA21_slope" )
   
    dimids_nor = (/ x_dimid, y_dimid, t_dimid /)
    ct_nor = (/ nx, ny, nt /)
    start_nor =  (/ 1, 1, 1 /) 
        
    ct_x = (/ nx /)
    start_x = (/ 1 /)
    
    ct_y = (/ ny /)
    start_y = (/ 1 /)
    
     ct_t = (/ nt /)
     start_t = (/ 1 /)
      

    !Initiailize vars
    call nc_check( nf90_def_var(ncIDo, "ASMASK3", nf90_real, dimids_nor, varid), "kasugaEA21", "def var asmask23kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", " "), "kasugaEA21", "def units asmask3 kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "mask for So maxima after shared zmin checks"), "kasugaEA21", "def units asmask3 kasugaEA21_slope" ) 
    
    call nc_check( nf90_def_var(ncIDo, "MASKRo", nf90_real, dimids_nor, varid), "kasugaEA21", "def var MASKRo kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "binary"), "kasugaEA21", "def units MASKRo kasugaEA21_slope" )  
    call nc_check( nf90_put_att(ncIDo, varid, "description ", "binary mask for radius of features after all checks"), "kasugaEA21", "def units MASKRo kasugaEA21_slope" )       
    
    
    call nc_check( nf90_def_var(ncIDo, "time", nf90_real, t_dimid, varid), "kasugaEA21", "def var time kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "seconds since 1970-01-01 00:00:00.0 0:00"), "kasugaEA21", "def units time kasugaEA21_slope" )  
    call nc_check( nf90_def_var(ncIDo, "latitude", nf90_real, y_dimid, varid), "kasugaEA21", "def var latitude kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "degrees_north"), "kasugaEA21", "def units latitude kasugaEA21_slope" )  
    call nc_check( nf90_def_var(ncIDo, "longitude", nf90_real, x_dimid, varid), "kasugaEA21", "def var longitude kasugaEA21_slope" )
    call nc_check( nf90_put_att(ncIDo, varid, "units ", "degrees_east"), "kasugaEA21", "def units longitude kasugaEA21_slope" )  
   
    call nc_check( nf90_enddef(ncIDo), "kasugaEA21", "end define mode kasugaEA21_slope" ) 


    !Write Vars
    call nc_check( nf90_inq_varid(ncIDo, "ASMASK3", varid), "kasugaEA21", "inq_varid asmask3")
    call nc_check( nf90_put_var(ncIDo, varid, asmask3, start_nor, ct_nor), "kasugaEA21", "put var asmask3 kasugaEA21_slope" )
    
    call nc_check( nf90_inq_varid(ncIDo, "MASKRo", varid), "kasugaEA21", "inq_varid MASKRo")
    call nc_check( nf90_put_var(ncIDo, varid, gridmask_all, start_nor, ct_nor), "kasugaEA21", "put var MASKRo kasugaEA21_slope" )
    
    
    call nc_check( nf90_inq_varid(ncIDo, "time", varid), "kasugaEA21", "inq_varid lat")
    call nc_check( nf90_put_var(ncIDo, varid, time, start_t, ct_t), "kasugaEA21", "put var time kasugaEA21_slope" )
    call nc_check( nf90_inq_varid(ncIDo, "latitude", varid), "kasugaEA21", "inq_varid lat")
    call nc_check( nf90_put_var(ncIDo, varid, lat, start_y, ct_y), "kasugaEA21", "put var lat kasugaEA21_slope" )
    call nc_check( nf90_inq_varid(ncIDo, "longitude", varid), "kasugaEA21", "inq_varid lon")
    call nc_check( nf90_put_var(ncIDo, varid, lon, start_x, ct_x), "kasugaEA21", "put var lon kasugaEA21_slope" )  

    ! Close the output file
    call nc_check( nf90_close(ncIDo), "kasugaEA21", "close " // trim(doutfile) )
  
  
call cpu_time(FN)
print '("Time for netcdf= ",f6.3," seconds.")',FN-ST
end if ! End binary mode

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!https://rosettacode.org/wiki/Haversine_formula#Fortran

function to_radian(degree) result(rad)
          ! degrees to radians
          real,intent(in) :: degree
          real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
          real :: rad
 
          rad = degree*deg_to_rad
end function to_radian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
          ! great circle distance -- adapted from Matlab 
          real,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real :: a,c,dist,dlat,dlon,lat1,lat2
          real,parameter :: radius = 6372.8 
 
          dlat = to_radian(deglat2-deglat1)
          dlon = to_radian(deglon2-deglon1)
          lat1 = to_radian(deglat1)
          lat2 = to_radian(deglat2)
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          c = 2*asin(sqrt(a))
          dist = radius*c
end function haversine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!https://diego.assencio.com/?index=8d6ca3d82151bad815f78addf9b5c1c6
!Compute the percent overlap between two circles (relative to the smaller circle)
!Value of 1.0 means the smaller circle is englufed completely by the larger

function c1c2_overlap(distance,r1,r2) result (overlap)
          real,intent(in) 	:: distance,r1,r2
          real 			:: tempvar,dist1,dist2,area1,area2,radius1,radius2
	  real			:: toofar, within, maxover
 	  real, parameter	:: PI = 4.0*atan(1.0)
	  
	  radius1 = r1
	  radius2 = r2
          if(radius2.gt.radius1)then
	    !print *,radius1
	    !print *,radius2
	    tempvar = radius2
	    radius2 = radius1
	    radius1 = tempvar
	    !print *,radius1
	    !print *,radius2
	    !stop
	  end if
	
	  toofar = radius1+radius2
	  within = abs(radius1-radius2)
	  maxover = PI*(min(radius2,radius1)**2.0)
	  if(distance.ge.toofar)then
	    overlap = 0		!At most overlap at a single point
	  else if(distance.le.within)then
	    overlap = maxover/maxover
          else
	  
	    dist1 = ((radius1**2.0)-(radius2**2.0)+(distance**2.0))/(2.0*distance)
	    dist2 = distance-dist1
	  
	    area1 = ((radius1**2.0)*acos(dist1/radius1))-(dist1*(((radius1**2.0)-(dist1**2))**(0.5)))
 	    area2 = ((radius2**2.0)*acos(dist2/radius2))-(dist2*(((radius2**2.0)-(dist2**2))**(0.5)))

	    overlap = (area1+area2)/maxover
	  end if
end function c1c2_overlap



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nc_check(istatus, subr_name, context)

use netcdf

implicit none

integer,          intent (in)          :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context

character(len=129) :: error_msg

! if no error, nothing to do here.  we are done.
if( istatus == nf90_noerr) return

! something wrong.  construct an error string
if (present(context) ) then
   error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))
else
   error_msg = nf90_strerror(istatus)
endif

write(6,*) error_msg
stop

end subroutine nc_check




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SLOPE2D - Compute the 2D slope function following Kasuga et al (2021)
!
!INPUT:
!  radii	list of radii over which to calculate slope	[m]
!  lat		latitudes 					[deg N]
!  lon		longitudes					[deg E]
!  hgt		geopotential height				[m]
!  
!  ny		size of grid in y direction			[]
!  nx		size of grid in x direction			[]
!  nr		number of radii to check			[]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine SLOPE2D(radii, lat, lon, hgt, ny, nx, nr)

!implicit none

!integer, intent(in)	:: ny, nx, nr
!real, intent(in)	:: hgt(ny,nx)
!real, intent(in)	:: lat(ny), lon(nx), radii(nr)
!real, intent(out)	:: slope(ny,nx)

!integer			:: r, j

!do r = 1, nr
  
!  do j = 0, 20
!end
