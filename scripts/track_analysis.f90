!track_kasugaEA21
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Kevin Lupo
!klupo@ucar.edu
!track cutoff lows and preexisting troughs identified using the KS21 algorithm
!
!21 Oct 2021
!	Simple check if feature(s) at time(f) is within Ro(s or t) of feature(t) at time(f-1)...
!	...If within distance, assign the (f) feature the same ID as the (f-1) feature
!	...Works better than expected, some ambiguity, but a good starting point
!22 Oct 2021
!	In cases where multiple features are within Ro, take the closest
!27 Oct 2021
!	Defined a penalty function to determine if a feature is likely found at subsequent times
!	Test predicting the location of a cutoff at the next time step by estimating...
!	...the background geostrophic wind from BGo (ug=-g/f*dz/dy; vg=g/f*dz/dx)
!	Geostrophic wind test can work, but tabled for now
!	Changed criteria to be based on circle overlap
!29 Oct 2021
!	Moving to version2...track only cutoffs that have a So max of at least 15 m/100km or troughs with So ge 10 m/100km
!4 Nov 2021 
!	Many changes - back to 10 & 5 So mins for cutoffs and troughs
!		- Penalty terms for change in So, BG, Ro, dx, dy
!5 Nov 2021
!	Moving to version3...more strict elimination criteria (e.g., https://journals.ametsoc.org/view/journals/clim/14/18/1520-0442_2001_014_3863_votnac_2.0.co_2.xml)
!8 Nov 2021 (relevant changes to general algorithm code)
! 	Algorithm wrapper script will now generate a dummy text file if model input is missing (should not have substantial effect on climo, but will have some tracking implications)
!	Fixed a bug that didn't output enough digits if So > 99.9 m/100km. This caused problems with tracking IO and prematurely terminated the output file for a given hour. Likely associated with TCs
!	For cleaner tracking, keep only features with So > 5.0 m/100km. This is consistent with the climatology of KasugaEA21!	
!11 Nov 2021
! 	To save time in analysis code, a feature's DX, DY, and DISTance travelled between ID times is saved to the output file
!	"...", feature's DT between ID times is also saved (relevant if a time is skipped due to missing files, feature not ID'd...)
!13 Nov 2021
!	Adapted to read more information from .dat file (850,500,200 Z, T, U, V, RH)
!17 Nov 2021
!	Added temporary fix to deal with bug in algorithm output where BGo = SR
!10 Dec 2021
!	BGo bug is fixed
!	Data now in version6
!14 Dec 2021
!	Added feature duration to output
!15 Dec 2021
!	Added FTARGET track quality metric to output (PTY-OVR)
!	Added placeholder columns for forecast error (ERROR-Y, ERROR-X, ERROR-TOT)
!16 Dec 2021
!	For analysis purposes, added a maxval(pack()) step to find the maximum duration of a feature. This will help to only analyze features that existed beyond a given threshold, without requiring the plotting/analysis code to find these maxima itself
!28 Dec 2021
!	Bug identified that erroneously does not remove duplicate IDs from certain times (for example, ID14900 in 2016022906.f000). Not only does this cause errors in the analysis tracks, but also forecast tracks. Unclear how pervasive this error is in forecasts or exactly what is causing it...
!30 Dec 2021
!	Bug was in the loop that checked for duplicate IDs. ID(n,s) could be reset to a negative value by setting ID(n,t)=-1*ID(n,s) when t=s, which cased problems if even more IDs matched ID(n,s). Added "TESTID" integer variable to preserve the original value of ID(n,s)
!14 Jan 2022
! 	Commented out running averages. They are unused.
!9 Feb 2022
!	.dat files are now split between northern and southern hemispheres for efficiency. Tracking is now done seperately for each hemisphere. Some changes necessary to directory and file references in wrapper .csh script and .f90 code
!10 Feb 2022
!	Added flags to indicate times of max So and min Z500
!28 April 2022
!	Output redirected to version 8
!	Increased "direction change" threshold to 1000km x and y
!	Revised penalty terms (see notes in section)
!3 May 2022
!	Removed hemisphere separation. If any exist, features are now permitted to cross the equator and be tracked poleward of 70Lat
!18 May 2022
!	Added a check that uses the hstry_dt variable to determine if more than 1 of the last 4 instances of an ID were missed. In this case, tracking of that ID is discontinued to prevent ambiguous/erroneously long tracks
!27 June 2022
!	Move to version9 following corrections to ID code related to max radius limit
!	Expanded input datastream to include fixed radius 500hPa Z,T,U,V,MR
!28 June 2022
! 	Testing configuration that uses the penalty term as the "Target" only (no combined penalty-overlap...why should a potentially highly-penalized feature be considered simply because it's radius is smaller?)
!		- Commit to this configuration. It's simpler, and has no impact on the three 2020 cutoffs (actually improves the starting point of the 4 April cutoff over Alaska by removing "sharp" motions)
!	Read user selected penalty norms and penalty max from wrapper script
!	Test configuration that reduces pmax to 1.0
!29 June 2022
!	Wrapper script uses a variety of user-selectable inputs to tune the tracking scheme
!22 Sep 2022
!	Fixed displacement error placeholders and missing values to allow for errors > 10000km
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
use netcdf

implicit none
real, parameter			:: GRAV = 9.81 ! m/s
real, parameter			:: ROT = 7.2921*(10.0**-5) ! rad/s
real, parameter			:: DT = 6.0*3600.0 ! s

integer, parameter		:: npass = 4
integer				:: ipass

real, allocatable		:: COR(:,:),Plat(:,:),Plon(:,:),vg(:,:),ug(:,:)
real				:: pdla, pdlo
real				:: apdist, apxrat, apyrat
integer				:: nx, ny, nf, ns, nt, nsmax, s, f, t, tt, q, nadlalo
integer				:: i, ii, j, k, x, y, r, a, b, np1, np2, nz1, np, g, c, n
integer				:: i0, i1, j0, j1, xx0, xx1
integer				:: istat,istat2
real				:: d, haversine, to_radian, c1c2_overlap, hgtmin, dmin, dthresh, p, pmin, o, smin, sm
real				:: ST, FN
real				:: dlat, dlon, PTYy, PTYx
integer, allocatable		:: ID(:,:),nid(:,:)
integer				:: TESTID
real, allocatable		:: d2prev(:,:),PTY(:,:), o2prev(:,:), o2next(:,:),acd2prev(:,:),ad2prev(:,:),acdlat2prev(:,:),adlat(:,:),acdlon2prev(:,:),adlon(:,:),dlon2prev(:,:),dlat2prev(:,:)
real, allocatable		:: So(:,:), Slat(:,:), Slon(:,:), Ro(:,:), SR(:,:), BGo(:,:), BGoy(:,:), BGox(:,:), Zmin(:,:), Zlat(:,:), Zlon(:,:)
real, allocatable		:: z850(:,:), z500(:,:), z200(:,:), t850(:,:), t500(:,:), t200(:,:), u850(:,:), u500(:,:), u200(:,:), v850(:,:), v500(:,:), v200(:,:), mr850(:,:), mr500(:,:), mr200(:,:)
real, allocatable		:: FIXEDz500(:,:), FIXEDt500(:,:), FIXEDu500(:,:),FIXEDv500(:,:),FIXEDmr500(:,:)

real, allocatable		:: lat(:), lon(:), radius(:),lat2d(:,:),lon2d(:,:)
integer, allocatable		:: isFMaxSo(:,:),isFMinZ(:,:),idxSo(:),idxZ(:)



!real, allocatable		:: runav_So(:,:),runav_BG(:,:),runav_Ro(:,:)
real, allocatable		:: hstry_So(:,:,:),hstry_BG(:,:,:),hstry_Ro(:,:,:)
!real, allocatable		:: runav_la(:,:),runav_lo(:,:),runav_dx(:,:),runav_dy(:,:)
real, allocatable		:: hstry_la(:,:,:),hstry_lo(:,:,:),hstry_dx(:,:,:),hstry_dy(:,:,:)
real, allocatable		:: saveDX(:,:),saveDY(:,:),saveDIST(:,:),saveDT(:,:)
real, allocatable		:: hstry_dt(:,:,:),saveDUR(:,:)
integer				:: EXCD
real				:: P_SLOPE,P_BACKG,P_RADII,P_DX,P_DY

real,allocatable		:: FERRX(:,:),FERRY(:,:),FERR(:,:)

real					:: CSOMIN,TSOMIN,OMIN,PMAX,DMAX,EMAX,OXMAX,OYMAX,SNORM,BNORM,RNORM
real					:: OVERLAP,DISTANCE
real					:: PDX,PDY,DX,DY,EDIST,OYRAT,OXRAT
real					:: SPTY,BPTY,RPTY,PSUM
real					:: TARGET
real, allocatable			:: FTARGET(:,:),saveFTARG(:,:)
real, allocatable			:: maxID(:,:)
integer					:: prior

real				:: user_pmax,user_normso,user_normbgo,user_normro,user_emax,user_xoppmax,user_yoppmax,user_dmax



character(len=80)		:: name
character(len=270)		:: list
character(len=300)		:: line, line2
character(len=256)		:: infile
character(len=20)		:: varstr
character(len=10)		:: hemis
character(len=5)		:: hext
character(len=10)		:: YYYY
character(len=10)		:: MM
character(len=5)		:: DD
character(len=5)		:: hh
character(len=20),allocatable   :: fhour(:,:)
character(len=256),allocatable	:: outfile(:)
character(len=400)		:: fstr
character(len=20),allocatable	:: itime(:,:)
character(len=270)		:: outdir
character(len=20)		:: didinc
character(len=10)		:: st_pmax,st_normso,st_normbgo,st_normro,st_emax,st_xoppmax,st_yoppmax,st_dmax

!!!!!!Get some command line arguments
call get_command_argument(1,list)
call get_command_argument(2,varstr)
call get_command_argument(3,st_pmax)
call get_command_argument(4,st_normso)
call get_command_argument(5,st_normbgo)
call get_command_argument(6,st_normro)
call get_command_argument(7,st_emax)
call get_command_argument(8,st_xoppmax)
call get_command_argument(9,st_yoppmax)
call get_command_argument(10,st_dmax)

read(st_pmax    , * ) user_pmax
read(st_normso  , * ) user_normso
read(st_normbgo , * ) user_normbgo
read(st_normro  , * ) user_normro
read(st_emax    , * ) user_emax   
read(st_xoppmax , * ) user_xoppmax
read(st_yoppmax , * ) user_yoppmax
read(st_dmax    , * ) user_dmax   

!call get_command_argument(3,hemis)
!call get_command_argument(4,hext)


!!!!!Name the output text file
outdir = "/glade/work/klupo/postdoc/kasugaEA21/version9/" // trim(varstr) // "/"



	
!!!!!Get dimensions for variable lists!!!!!
call cpu_time(ST)
nf = 0
nsmax = 0
open(unit = 10, file = trim(list), status='old', access='sequential', form='formatted', action='read')
do
  read(10,'(a256)',iostat=istat) infile     
  if(istat /= 0)exit
  
  ns=0
  open(unit = 11, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do
    if(ns.eq.0)then
      read(11,*,iostat=istat2) line
    else
      read(11,*,iostat=istat2) line
    end if    
    if(istat2 /=0)exit
    ns=ns+1
    nsmax = max(ns,nsmax)
  end do
  close(11)
  
  nf=nf+1
end do
close(10)
call cpu_time(FN)
print '("Time to get max file dimensions = ",f6.3," seconds.")',FN-ST


allocate(itime(nf,nsmax),fhour(nf,nsmax),ID(nf,nsmax),So(nf,nsmax),Slat(nf,nsmax),Slon(nf,nsmax),Ro(nf,nsmax),SR(nf,nsmax),BGo(nf,nsmax),BGoy(nf,nsmax),BGox(nf,nsmax),Zmin(nf,nsmax),Zlat(nf,nsmax),Zlon(nf,nsmax),d2prev(nf,nsmax),o2prev(nf,nsmax),o2next(nf,nsmax),PTY(nf,nsmax))
allocate(z850(nf,nsmax), z500(nf,nsmax), z200(nf,nsmax), t850(nf,nsmax), t500(nf,nsmax), t200(nf,nsmax), u850(nf,nsmax), u500(nf,nsmax), u200(nf,nsmax), v850(nf,nsmax), v500(nf,nsmax), v200(nf,nsmax), mr850(nf,nsmax), mr500(nf,nsmax), mr200(nf,nsmax))
allocate(FIXEDz500(nf,nsmax), FIXEDt500(nf,nsmax), FIXEDu500(nf,nsmax),FIXEDv500(nf,nsmax),FIXEDmr500(nf,nsmax))
allocate(saveDX(nf,nsmax),saveDY(nf,nsmax),saveDIST(nf,nsmax),saveDT(nf,nsmax))
allocate(FTARGET(nf,nsmax))
allocate(COR(nf,nsmax),Plat(nf,nsmax),Plon(nf,nsmax),vg(nf,nsmax),ug(nf,nsmax))
allocate(outfile(nf))
allocate(nid(nf,nsmax))
!allocate(runav_So(nf,nsmax),runav_BG(nf,nsmax),runav_Ro(nf,nsmax))
allocate(hstry_So(nf,nsmax,300),hstry_BG(nf,nsmax,300),hstry_Ro(nf,nsmax,300))
!allocate(runav_la(nf,nsmax),runav_lo(nf,nsmax),runav_dx(nf,nsmax),runav_dy(nf,nsmax))
allocate(hstry_la(nf,nsmax,300),hstry_lo(nf,nsmax,300),hstry_dx(nf,nsmax,300),hstry_dy(nf,nsmax,300))
allocate(hstry_dt(nf,nsmax,300),saveDUR(nf,nsmax))
allocate(saveFTARG(nf,nsmax))
allocate(FERRX(nf,nsmax),FERRY(nf,nsmax),FERR(nf,nsmax))
allocate(maxID(nf,nsmax))
allocate(isFMaxSo(nf,nsmax),isFMinZ(nf,nsmax),idxSo(2),idxZ(2))

!!!!!Set dummy values for lists!!!!!
itime = "-9999"
fhour = "-9999"
ID = -1
So = -9999.9
Slat = -9999.9
Slon = -9999.9
Ro = -9999.9
SR = -9999.9
BGo = -9999.9
BGoy = -9999.9
BGox = -9999.9
Zmin = -9999.9
Zlat = -9999.9
Zlon = -9999.9
saveDX = -9999.9
saveDY = -9999.9
saveDIST = -9999.9
saveDT = 0.0
saveDUR = 0.0
saveFTARG = 9999.9
FTARGET = 9999.9
hstry_dt = 0.0
COR = -9999.9
ug = -9999.9
vg = -9999.9
nid = 0
FERRY = -99999.9
FERRX = -99999.9
FERR = -99999.9
maxID = -9999.9
isFMaxSo = 0
isFMinZ = 0


call cpu_time(ST)
!!!!!Read Vars from .dat files!!!!!
open(unit = 10, file = trim(list), status='old', access='sequential', form='formatted', action='read')
do f=1,nf
  read(10,'(a256)',iostat=istat) infile     
  if(istat /= 0)exit
  
  open(unit = 11, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do s=0,nsmax
    if(s.eq.0)then
      read(11,*,iostat=istat2) line
    else
      read(11,*,iostat=istat2) itime(f,s), fhour(f,s), So(f,s), Slat(f,s), Slon(f,s), Ro(f,s), SR(f,s), BGo(f,s), BGoy(f,s), BGox(f,s), Zmin(f,s), Zlat(f,s), Zlon(f,s), z850(f,s), z500(f,s), z200(f,s), t850(f,s), t500(f,s), t200(f,s), u850(f,s), u500(f,s), u200(f,s), v850(f,s), v500(f,s), v200(f,s), mr850(f,s), mr500(f,s), mr200(f,s),FIXEDz500(f,s), FIXEDt500(f,s), FIXEDu500(f,s),FIXEDv500(f,s),FIXEDmr500(f,s)
    end if    
    if(istat2 /=0)exit
  end do
  close(11)
  outfile(f) = trim(infile(1:37)) // "version9" // trim(infile(46:80)) // ".track" 
  
end do
close(10)
call cpu_time(FN)
print '("Time to read data = ",f6.3," seconds.")',FN-ST


!!!!!Estimate the predicted points!!!!!
!do f=1,nf
!  do s=1,nsmax
!    if(itime(f,s).eq."-9999") exit
!    COR(f,s) = 2.0*ROT*to_radian(Slat(f,s))
!    vg(f,s) = (GRAV/COR(f,s))*(BGox(f,s)*(1/100000.0))
!    ug(f,s) = (-1.0)*(GRAV/COR(f,s))*(BGoy(f,s)*(1/100000.0))
!    pdx = ug(f,s)*DT
!    pdy = vg(f,s)*DT
    
!    d = haversine(Slat(f,s),Slon(f,s),Slat(f,s)+1.0,Slon(f,s))*1000.0
!    pdla = pdx/d
!    d = haversine(Slat(f,s),0.0,Slat(f,s),1.0)*1000.0
!    pdlo = pdy/d
    
!    Plat(f,s) = Slat(f,s)+pdla
!    Plon(f,s) = Slon(f,s)+pdlo
!    if(Plon(f,s).ge.360.0)then
!      Plon(f,s) = Plon(f,s)-360.0
!    end if
!    if(Plon(f,s).lt.0.0)then
!      Plon(f,s) = 360.0-Plon(f,s)
!    end if
    
    
!  end do
!end do


!!!!!Find cutoffs and troughs meeting minimum So criteria (C So ge 10, T So ge 5)!!!!!
!!!!!Match current time (n) features (s) to previous time (n-1) features (t)!!!!!
!!!!!If n=1, all qualifying features assigned a unique integer ID (i), i incremented by 1!!!!!
!!!!!!!!!!FAIL TEST - So is too small (ONLY IF IT IS THE FIRST TIME A FEATURE IS IDENTIFIED)!!!!!
!!!!!If n>1, features (s) at time (n) are compared to features (t) at time (n-1) to determine which (if any) n-1 feature matches time (n) feature.
!!!!!If a match is made, ID (i) is carried forward to feature (n,s).
!!!!!!!!!!FAIL TEST - Optimal radii of features (n,s) and (n-1,t) do not overlap
!!!!!!!!!!FAIL TEST - Distance between features exceeds 1500 km
!!!!!!!!!!If number of occurances of ID (i) is at least 2, use position history to estimate the expected position of the feature at time (n)
!!!!!!!!!!!!!!!FAIL TEST - Position of feature (n,s) is more than 1000km different than the expected position
!!!!!!!!!!!!!!!FAIL TEST - delta-x between feature (n,s) and feature (n-1,t) is greater than 300 km in the opposite direction of the expected position
!!!!!!!!!!!!!!!FAIL TEST - delta-y between feature (n,s) and feature (n-1,t) is greater than 600 km in the opposite direction of the expected position
!!!!!!!!!!PENALTY - Large change in So between (n-1,t) and (n,s) ---- PTY = abs(So(n,s)-So(n-1,t))/10.0
!!!!!!!!!!PENALTY - Large change in BGo between (n-1,t) and (n,s) ---- PTY = abs(BGo(n,s)-BGo(n-1,t))/10.0
!!!!!!!!!!PENALTY - Large change in Ro between (n-1,t) and (n,s) ---- PTY = abs(Ro(n,s)-Ro(n-1,t))/900.0
!!!!!!!!!!FAIL TEST - Total penalty term is greater than 1.5
!!!!!!!!!!If multiple features (n-1,t) pass all tests - Keep the one with the smallest value of PENALTY-%OVERLAP
!!!!!!If the same (n-1,t) feature is a best match to multiple (n,s) features, keep the match with feature with the longest history
!!!!!!!!!!If the same history, keep the match with the smallest value of PENALTY-%OVERLAP
!!!!!!!!!!Reassign the "unmatched" feature (n,s) an ID of -1*ID(n-1,t) so that it may be tested on subsequent pass to determine if it is an acceptable match for a different point 
!!!!!!If a feature (n,s) has no match from (n-1,t), repeat the above examination for features at (n-2,t) that were not propagated to time (n-1)



!!!!!Set some failure thresholds
CSOMIN = 5.0		! Minimum So for a cutoff to be identified
TSOMIN = 5.0		! Minimum So for a trough to be identified
OMIN   = 0.0		! Minimum overlap percentage for matching features
PMAX   = user_pmax	! 1.5		! Maximum allowable penalty points before a feature is disqualified
DMAX   = user_dmax	! 1500		! Maximum distance between matches allowed
EMAX   = user_emax	! 1500		! Maximum "error" distance between predicted and candidate points
OXMAX  = user_xoppmax	! 1000.0		! Maximum distance a matching point can be in the opposite x-direction of the predicted point (revised from 300km to 1000km since this seems to have caused problems with the september 2020 colorado cutoff track
OYMAX  = user_yoppmax	! 1000.0		! ".."y-direction"" More lenient to allow for features moving around larger scale troughs to move to the north (revised from 600, same as above)
SNORM  = user_normso	! 15.0		! Normalization for the So penalty 	! Motivated by Sep 2020 case, changed to 15 from 10)
BNORM  = user_normbgo	! 40.0		! Normalization for the BGo penalty	! motivated by Sep 2020 case, changed to 40 from 15)
RNORM  = user_normro	! 1000.0		! Normalization for the Ro penalty	! motivated by the Sep 2020 case, changed to 1000 from 900

call cpu_time(ST)
i=1
do n=1,nf
  !print *,n
  do ipass=1,npass
    do s=1,nsmax   	
      if(itime(n,s).eq."-9999") exit									! Break the (s) loop since it has reached the end of the feature list
      if(ID(n,s).gt.0) cycle										! Skip this (s) because it already has an ID from a previous pass
      
      if(n.gt.1) then
        prior = n-1											! Start matching points
        TARGET = 9999.9
        do t=1,nsmax
	  if(itime(prior,t).eq."-9999") exit								! Break the (t) loop since it has reached the end of the feature list
          if(ID(prior,t).lt.0) cycle									! Skip this (t) because it has no ID
	  if(-1*ID(prior,t).eq.ID(n,s)) cycle								! Skip this (t) because (n-1,t) better matches a different (n,s), but is also the the best match to this (n,s) 
	  
	  if(count(hstry_dt(prior,t,(max(1,nid(prior,t)-3)):nid(prior,t)).eq.12).gt.1) cycle		! No longer tracking a feature if at least 2 of the last 4 times were skipped
	  
	  DISTANCE 	= haversine(Slat(n,s),Slon(n,s),Slat(prior,t),Slon(prior,t))
	  DX 		= haversine(Slat(prior,t),Slon(n,s),Slat(prior,t),Slon(prior,t))
	  DY 		= haversine(Slat(n,s),Slon(prior,t),Slat(prior,t),Slon(prior,t))	  
													! Get the signs of the test point motion
	  if((Slat(n,s)-hstry_la(prior,t,nid(prior,t))).lt.0) DY = -1.0*DY	  
	  if(((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).lt.0).and.((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.-100.0)) DX = -1.0*DX
	  if(((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.0).and.((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.100.0))  DX = -1.0*DX
	  
	  
	  
	  if(DISTANCE.gt.DMAX) cycle									! Skip this (t) because it is too far from (n,s)
	  OVERLAP 	= c1c2_overlap(DISTANCE,Ro(n,s),Ro(prior,t))
	  
	  
	  
	  if(OVERLAP.le.OMIN) cycle									! Skip this (t) because it does not overlap with (n,s)
	  
	  if(nid(prior,t).gt.1) then									! Check DX and DY against predicted motion
	    PDX 	= abs(hstry_dx(prior,t,nid(prior,t)))
	    PDY 	= abs(hstry_dy(prior,t,nid(prior,t)))
	    
	    												! Get the signs of the predicted motion
	    if((hstry_la(prior,t,nid(prior,t))-hstry_la(prior,t,nid(prior,t)-1)).lt.0) PDY = -1.0*PDY	  
	    if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).lt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.-100.0)) PDX = -1.0*PDX
	    if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.100.0))  PDX = -1.0*PDX
	  												
													
	    
	    EDIST 	= ((PDX-DX)**2+(PDY-DY)**2)**(0.5)
	    OXRAT 	= DX/PDX
	    OYRAT 	= DY/PDY
	    
	    !if(n.eq.29.and.s.eq.8.and.t.eq.7) then
	    !  print *,DISTANCE
	    !  print *,OVERLAP
	    !  print *,PSUM
	    !  print *,abs(PDX-DX)
	    !  print *,abs(PDY-DY)
	    !  print *,EDIST
	    !end if
	        
	    if(OXRAT.lt.0.0.and.abs(PDX-DX).gt.OXMAX) cycle							! Skip this (t) because it moved too far in the opposite of the expected x-direction
	    if(OYRAT.lt.0.0.and.abs(PDY-DY).gt.OYMAX) cycle							! Skip this (t) because it moved too far in the opposite of the expected y-direction
	    if(EDIST.gt.EMAX) cycle									! Skip this (t) because it is too far from its expected position
	  end if
	  
	  SPTY 		= abs(So(n,s)-So(prior,t))/SNORM						! Penalize this (t) based on the difference between it and So(n,s)
	  BPTY 		= abs(BGo(n,s)-BGo(prior,t))/BNORM						! Penalize this (t) based on the difference between it and BGo(n,s)
	  RPTY 		= abs(Ro(n,s)-Ro(prior,t))/RNORM						! Penalize this (t) based on the difference between it and Ro(n,s)	  
	  PSUM 		= (1.0*SPTY)+(1.0*BPTY)+(1.0*RPTY)						! Some cases (e.g., Sep 2020 colorado) demonstrate that the backgroudn term has to much influence. Weights added here
	  
	  
	  !if((PSUM-OVERLAP).gt.PMAX) cycle									! Skip this (t) because it accumulated too many penalty points
	  if(PSUM.gt.PMAX) cycle  
	  
	  !if((PSUM-OVERLAP).lt.TARGET) then								! Compare the combined penalty & overlap term to the previous "TARGET" combined term to beat
	  if(PSUM.lt.TARGET) then
	    ID(n,s)	= ID(prior,t)									! Assign a new history and ID to feature (n,s)
	    nid(n,s) 	= nid(prior,t) + 1  								! Assign characteristic history
	    hstry_So(n,s,:) 		= hstry_So(prior,t,:)
            hstry_BG(n,s,:) 		= hstry_BG(prior,t,:)
            hstry_Ro(n,s,:) 		= hstry_Ro(prior,t,:)
	    hstry_So(n,s,nid(n,s)) 	= So(n,s)
            hstry_BG(n,s,nid(n,s)) 	= BGo(n,s)
            hstry_Ro(n,s,nid(n,s)) 	= Ro(n,s)
	  												! Assign position history
	    hstry_la(n,s,:) 		= hstry_la(prior,t,:)
            hstry_lo(n,s,:) 		= hstry_lo(prior,t,:)
            hstry_dx(n,s,:) 		= hstry_dx(prior,t,:)
	    hstry_dy(n,s,:) 		= hstry_dy(prior,t,:)	  
	    hstry_la(n,s,nid(n,s)) 	= Slat(n,s)
            hstry_lo(n,s,nid(n,s)) 	= Slon(n,s)
            hstry_dx(n,s,nid(n,s)) 	= DX
	    hstry_dy(n,s,nid(n,s)) 	= DY
	  												! Compute running averages of histories
	    !runav_So(n,s) 		= sum(hstry_So(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 
	    !runav_BG(n,s) 		= sum(hstry_BG(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 
	    !runav_Ro(n,s) 		= sum(hstry_Ro(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 

            !runav_la(n,s) 		= sum(hstry_la(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
            !runav_lo(n,s) 		= sum(hstry_lo(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
            !runav_dx(n,s) 		= sum(hstry_dx(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
            !runav_dy(n,s) 		= sum(hstry_dy(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
	    
	    saveDX(n,s) 		= DX
	    saveDY(n,s)			= DY
	    saveDIST(n,s)		= DISTANCE
	    saveDT(n,s)			= 6.0
	    hstry_dt(n,s,:)		= hstry_dt(prior,t,:)
	    hstry_dt(n,s,nid(n,s)) 	= saveDT(n,s)
	    saveDUR(n,s)		= sum(hstry_dt(n,s,1:nid(n,s)))
	    TARGET			= PSUM!-OVERLAP
	    FTARGET(n,s)		= TARGET
	    saveFTARG(n,s)		= TARGET
	  end if !TARGET if()
	end do ! t do()
	
	
	
	if(n.gt.2.and.ID(n,s).lt.0) then								! Check criteria at n-2 if an ID hasn't been propagated to n from n-1
	  prior = n-2
	  do t=1,nsmax
	    if(itime(prior,t).eq."-9999") cycle								! Skip this (n-2,t) because we have reached the end of the feature list at n-2
	    if(ID(prior,t).lt.0) cycle									! Skip this (n-2,t) because it does not have an ID
	    if(-1*ID(prior,t).eq.ID(n,s)) cycle								! Skip this (n-2,t) because it better matches a different (n,s), but is also the best match to this (n,s)
	    if(any(ID(n-1,:).eq.ID(prior,t))) cycle							! Skip this (n-2,t) because it was already propagated to (n-1)
	    
	    if(count(hstry_dt(prior,t,(max(1,nid(prior,t)-3)):nid(prior,t)).eq.12).gt.1) cycle		! No longer tracking a feature if at least 2 of the last 4 times were skipped

	    
	    DISTANCE 	= haversine(Slat(n,s),Slon(n,s),Slat(prior,t),Slon(prior,t))
	    DX 		= haversine(Slat(prior,t),Slon(n,s),Slat(prior,t),Slon(prior,t))
	    DY 		= haversine(Slat(n,s),Slon(prior,t),Slat(prior,t),Slon(prior,t))
													! Get the signs of the test point motion
	    if((Slat(n,s)-hstry_la(prior,t,nid(prior,t))).lt.0) DY = -1.0*DY	  
	    if(((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).lt.0).and.((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.-100.0)) DX = -1.0*DX
	    if(((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.0).and.((Slon(n,s)-hstry_lo(prior,t,nid(prior,t))).gt.100.0))  DX = -1.0*DX

	    if(DISTANCE.gt.DMAX) cycle									! Skip this (t) because it is too far from (n,s)
	  
	    OVERLAP 	= c1c2_overlap(DISTANCE,Ro(n,s),Ro(prior,t))
	    if(OVERLAP.le.OMIN) cycle									! Skip this (t) because it does not overlap with (n,s)
	  
	    if(nid(prior,t).gt.1) then									! Check DX and DY against predicted motion
	      PDX 	= abs(hstry_dx(prior,t,nid(prior,t)))
	      PDY 	= abs(hstry_dy(prior,t,nid(prior,t)))
	    
	    												! Get the signs of the predicted motion
	      if((hstry_la(prior,t,nid(prior,t))-hstry_la(prior,t,nid(prior,t)-1)).lt.0) PDY = -1.0*PDY	  
	      if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).lt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.-100.0)) PDX = -1.0*PDX
	      if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.100.0))  PDX = -1.0*PDX
	  												
	    
	      EDIST 	= ((PDX-DX)**2+(PDY-DY)**2)**(0.5)
	      OXRAT 	= DX/PDX
	      OYRAT 	= DY/PDY
	    	    
	      if(OXRAT.lt.0.0.and.abs(PDX-DX).gt.OXMAX) cycle						! Skip this (t) because it moved too far in the opposite of the expected x-direction
	      if(OYRAT.lt.0.0.and.abs(PDY-DY).gt.OYMAX) cycle						! Skip this (t) because it moved too far in the opposite of the expected y-direction
	      if(EDIST.gt.EMAX) cycle									! Skip this (t) because it is too far from its expected position
	    end if
	  
	    SPTY 		= abs(So(n,s)-So(prior,t))/SNORM					! Penalize this (t) based on the difference between it and So(n,s)
	    BPTY 		= abs(BGo(n,s)-BGo(prior,t))/BNORM					! Penalize this (t) based on the difference between it and BGo(n,s)
	    RPTY 		= abs(Ro(n,s)-Ro(prior,t))/RNORM					! Penalize this (t) based on the difference between it and Ro(n,s)	  
	    PSUM 		= (1.0*SPTY)+(1.0*BPTY)+(1.0*RPTY)						! Some cases (e.g., Sep 2020 colorado) demonstrate that the backgroudn term has to much influence. Weights added here
	    !if((PSUM-OVERLAP).gt.PMAX) cycle									! Skip this (t) because it accumulated too many penalty points
	    if(PSUM.gt.PMAX) cycle
	    !if((PSUM-OVERLAP).lt.TARGET) then								! Compare the combined penalty & overlap term to the previous "TARGET" combined term to beat
	    if(PSUM.lt.TARGET) then
	      ID(n,s)	= ID(prior,t)									! Assign a new history and ID to feature (n,s)
	      nid(n,s) 	= nid(prior,t) + 1
	      	    											! Assign characteristic history
	      hstry_So(n,s,:) 		= hstry_So(prior,t,:)
              hstry_BG(n,s,:) 		= hstry_BG(prior,t,:)
              hstry_Ro(n,s,:) 		= hstry_Ro(prior,t,:)
	      hstry_So(n,s,nid(n,s)) 	= So(n,s)
              hstry_BG(n,s,nid(n,s)) 	= BGo(n,s)
              hstry_Ro(n,s,nid(n,s)) 	= Ro(n,s)
	  												! Assign position history
	      hstry_la(n,s,:) 		= hstry_la(prior,t,:)
              hstry_lo(n,s,:) 		= hstry_lo(prior,t,:)
              hstry_dx(n,s,:) 		= hstry_dx(prior,t,:)
  	      hstry_dy(n,s,:) 		= hstry_dy(prior,t,:)	  
	      hstry_la(n,s,nid(n,s)) 	= Slat(n,s)
              hstry_lo(n,s,nid(n,s)) 	= Slon(n,s)
              hstry_dx(n,s,nid(n,s)) 	= DX
	      hstry_dy(n,s,nid(n,s)) 	= DY
	  												! Compute running averages of histories
	      !runav_So(n,s) 		= sum(hstry_So(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 
	      !runav_BG(n,s) 		= sum(hstry_BG(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 
	      !runav_Ro(n,s) 		= sum(hstry_Ro(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/min(nid(n,s),4) 

              !runav_la(n,s) 		= sum(hstry_la(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
              !runav_lo(n,s) 		= sum(hstry_lo(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
              !runav_dx(n,s) 		= sum(hstry_dx(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
              !runav_dy(n,s) 		= sum(hstry_dy(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
	    
	      saveDX(n,s) 		= DX
	      saveDY(n,s)		= DY
	      saveDIST(n,s)		= DISTANCE
	      saveDT(n,s)		= 12.0
	      hstry_dt(n,s,:)		= hstry_dt(prior,t,:)
	      hstry_dt(n,s,nid(n,s)) 	= saveDT(n,s)
	      saveDUR(n,s)		= sum(hstry_dt(n,s,1:nid(n,s)))
	      TARGET			= PSUM!-OVERLAP
	      FTARGET(n,s)		= TARGET
	      saveFTARG(n,s)		= TARGET
	    end if ! TARGET if()
	  end do ! t do()
	end if ! n gt 2 if()
	    
        if(count(ID(n,:).eq.ID(n,s)).gt.1.and.ID(n,s).gt.0) then					! Check if multiple (n,s) features matched the same (n-1,t) feature
          TESTID		= ID(n,s)
	  a = minloc(FTARGET(n,:),1,MASK = ID(n,:).eq.ID(n,s))
	  
	  do t=1,nsmax
	    if(ID(n,t).eq.ID(n,s).and.t.ne.a) then
	      ID(n,t) 		= -1*TESTID
	      nid(n,t) 		= 0
	      FTARGET(n,t) 	= 9999.9
	      d2prev(n,t) 	= 0.0
	      o2prev(n,t) 	= 0.0
	      hstry_So(n,t,:) 	= 0.0
	      hstry_BG(n,t,:) 	= 0.0
	      hstry_Ro(n,t,:) 	= 0.0
	      hstry_dx(n,t,:) 	= 0.0
	      hstry_dy(n,t,:) 	= 0.0
	      hstry_la(n,t,:) 	= 0.0
	      hstry_lo(n,t,:) 	= 0.0
	      !runav_So(n,t) 	= 0.0
	      !runav_BG(n,t) 	= 0.0
	      !runav_Ro(n,t) 	= 0.0
	      !runav_dx(n,t) 	= 0.0
	      !runav_dy(n,t) 	= 0.0
	      !runav_la(n,t) 	= 0.0
	      !runav_lo(n,t) 	= 0.0
	      saveDX(n,t) 	= -9999.9
	      saveDY(n,t)	= -9999.9
	      saveDIST(n,t)	= -9999.9
	      saveDT(n,t)	= 0.0
	      hstry_dt(n,t,:) 	= 0.0
	      saveDUR(n,t)	= 0.0
	      saveFTARG(n,t)	= 9999.9
	      
	    end if
          end do
	  if(ID(n,s).eq.-14900) then
	    print *,"IDs post fix"
	    print *,ID(n,:)
	  end if
        end if
      end if ! n gt 1 if()
      
      if(ID(n,s).lt.0.and.((Zmin(n,s).gt.0.0.and.So(n,s).ge.CSOMIN).or.(Zmin(n,s).lt.0.0.and.So(n,s).ge.TSOMIN))) then		! If no ID has been assigned, this is a new feature. Give it a new ID and increment integer (i)
      																! ....and initialize history arrays
        ID(n,s) 		= i
        i = i+1
        nid(n,s) 		= 1
	
	hstry_So(n,s,nid(n,s)) 	= So(n,s)
        hstry_BG(n,s,nid(n,s)) 	= BGo(n,s)
        hstry_Ro(n,s,nid(n,s)) 	= Ro(n,s)
      
        hstry_la(n,s,nid(n,s)) 	= Slat(n,s)
        hstry_lo(n,s,nid(n,s)) 	= Slon(n,s)
        hstry_dx(n,s,nid(n,s)) 	= 0.0
        hstry_dy(n,s,nid(n,s)) 	= 0.0
     
       ! runav_So(n,s) 		= sum(hstry_So(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
       ! runav_BG(n,s) 		= sum(hstry_BG(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
       ! runav_Ro(n,s) 		= sum(hstry_Ro(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))

        !runav_la(n,s) 		= sum(hstry_la(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
        !runav_lo(n,s) 		= sum(hstry_lo(n,s,(max(1,nid(n,s)-3)):nid(n,s)))/float(min(nid(n,s),4))
        !runav_dx(n,s) 		= sum(hstry_dx(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
        !runav_dy(n,s) 		= sum(hstry_dy(n,s,(max(2,nid(n,s)-2)):nid(n,s)))/float(min(nid(n,s),3))
      	
      end if  
      
      
      !if(n.gt.1.and.nid(n,s).eq.1) then 								! It is possible that the first identified instance of a feature may be preceded by a non-qualifying feature. Need to check prior times
      	! Since the ID algorithm currently uses a So minimum of 5.0 m/100km, it may not be necessary to do this step if the CSOMIN and TSOMIN are both set to 5.0
	! ...With this lower-limit, all features will qualify. Does place a heavier burden on the tracking scheme...
	! ...but should be much faster without having to re-check previous times.
      !end if			
      
      
      
      
      
             
    end do ! s do()
  end do ! npass do()
end do ! n do()
call cpu_time(FN)
print '("Time to track features = ",f6.3," seconds.")',FN-ST

call cpu_time(ST)
do n=1,nf
  do s=1,nsmax
    if(ID(n,s).gt.0) then
      maxID(n,s) = maxval(pack(saveDUR,ID.eq.ID(n,s)))
    end if
  end do
end do
call cpu_time(FN)
print '("Time to find max duration = ",f8.3," seconds.")',FN-ST

call cpu_time(ST)
do i=1,maxval(ID)
  idxSo = maxloc(So, mask = ID.eq.i)
  isFMaxSo(idxSo(1),idxSo(2)) = 1
  idxZ  = minloc(z500, mask = ID.eq.i)
  isFMinZ(idxZ(1),idxZ(2)) = 1
  !print *,i
  !print *,idxSo
end do
call cpu_time(FN)
print '("Time to apply maxSo and minZ flags = ",f6.3," seconds.")',FN-ST  


!!!!Write to text file!!!!
do f=1,nf
  open(10, file=trim(outfile(f)), status="unknown")
  
  fstr = "(A7,1x,A7,1x,A7,1x,A11,1x,A7,1x,A7,1x,A8,1x,A7,1x,A4,1x,A15,1x,A16,1x,A16,1x,A7,1x,A7,1x,A7,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A17,1x,A17,1x,A17,1x,A17,1x,A17,1x,A6,1x,A6,1x,A8,1x,A6,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)"
  write(10,fstr) "ITIME","FHOUR","ID","So(m/100km)","LAT(N)","LON(E)","SoFlag","Ro(km)","SR","BGo(m/100km)","BGo-lat(m/100km)","BGo-lon(m/100km)","ZMIN(m)","ZLAT(N)","ZLON(E)","ZFlag","Z850(m)","Z500(m)","Z200(m)","T850(K)","T500(K)","T200(K)","U850(m/s)","U500(m/s)","U200(m/s)","V850(m/s)","V500(m/s)","V200(m/s)","MR850(g/kg)","MR500(g/kg)","MR200(g/kg)","600kmZ500(m)","600kmT500(K)","600kmU500(m/s)","600kmV500(m/s)","600kmMR500(g/kg)","DY(km)","DX(km)","DIST(km)","DT(h)","DUR(h)","MAXDUR(h)","PTY-OVR", "FERRY(km)","FERRX(km)","FERR(km)"

  fstr = "(A11,1x,A4,1x,I9,1x,F6.2,1x,F6.2,1x,F6.2,1x,I3,1x,F10.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,I3,1x,F8.2,1x,F8.2,1x,F8.2,1x,F5.1,1x,F5.1,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.4,1x,F7.4,1x,F7.4,1x,F8.2,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.4,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F10.2,1x,F10.2,1x,F10.2)"
  do s=1,nsmax
    if(itime(f,s).eq."-9999") exit  
    write(10,fstr) itime(f,s), fhour(f,s), ID(f,s), So(f,s), Slat(f,s), Slon(f,s), isFMaxSo(f,s), Ro(f,s), SR(f,s), BGo(f,s), BGoy(f,s), BGox(f,s), Zmin(f,s), Zlat(f,s), Zlon(f,s), isFMinZ(f,s), z850(f,s), z500(f,s), z200(f,s), t850(f,s), t500(f,s), t200(f,s), u850(f,s), u500(f,s), u200(f,s), v850(f,s), v500(f,s), v200(f,s), mr850(f,s), mr500(f,s), mr200(f,s),FIXEDz500(f,s), FIXEDt500(f,s), FIXEDu500(f,s),FIXEDv500(f,s),FIXEDmr500(f,s), saveDY(f,s), saveDX(f,s), saveDIST(f,s), saveDT(f,s), saveDUR(f,s),maxID(f,s),saveFTARG(f,s),FERRY(f,s),FERRX(f,s),FERR(f,s)
  end do 
  close(10)
end do


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
