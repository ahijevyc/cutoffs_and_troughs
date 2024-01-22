!trackfx_kasugaEA21_version3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Kevin Lupo
!klupo@ucar.edu
!track cutoff lows and preexisting troughs identified using the KS21 algorithm through different forecast hours
!
!21 Oct 2021
!	Simple check if feature(s) at time(f) is within Ro(s or t) of feature(t) at time(f-1)...
!	...If within distance, assign the (f) feature the same ID as the (f-1) feature
!	...Works better than expected, some ambiguity, but a good starting point
!22 Oct 2021
!	In cases where multiple features are within Ro, take the closest
!-----------!
!
!26 Oct 2021
!	Begin reading forecast data (testing on f006-f012 in 2018)
!8 Nov 2021
!	Move to version 2 with many changes from analysis tracking code
!9 Nov 2021
!	Fixed memory issue due to wrapper script continuously appending to "ffile" list rather than refreshing the file. Much faster now.
!11 Nov 2021
! 	To save time in analysis code, a feature's DX, DY, and DISTance travelled between ID times is saved to the output file
!	"...", feature's DT between ID times is also saved (relevant if a time is skipped due to missing files, feature not ID'd...)
!13 Nov 2021
!	Adapted to read more information from .dat file (850,500,200 Z, T, U, V, RH)
!17 Nov 2021
!	Added temporary fix to manually calculate BGo from SR and So. Algorithm was erroneously outputting SR as the value of BGo
!10 Dec 2021
!	Fixed BGo bug
!	Data now is in /version6
!10 Dec 2021	
!	trackfx version 3
!	Developing capability for features to be tracked in forecasts initiallized prior to the feature existing in the analysis
!11-13 Dec 2021
!	Code -can- ID features from prior initialization times, however, these tracks are messy, difficult to interpret, and difficult to plot
!14 Dec 2021
!	Added limit to prevent features from being ID'd at fhours valid beyond the verification lifetime of a feature. This led to some messy plots at longer fhours
!15 Dec 2021
!	Added duration of feature to output. This incorporates the verification & forecast history
!	Added FTARGET track quality metric to output (PTY-OVR)
!	Added Distance errors (X, Y, Total) to output
!16 Dec 2021
!	Added output column that lists the maximum duration of the feature ID in the verification dataset
!28 Dec 2021
!	(From track) Bug identified that erroneously does not remove duplicate IDs from certain times (for example, ID14900 in 2016022906.f000). Not only does this cause errors in the analysis tracks, but also forecast tracks. Unclear how pervasive this error is in forecasts or exactly what is causing it...
!30 Dec 2021
!	(From track) Bug was in the loop that checked for duplicate IDs. ID(n,s) could be reset to a negative value by setting ID(n,t)=-1*ID(n,s) when t=s, which cased problems if even more IDs matched ID(n,s). Added "TESTID" integer variable to preserve the original value of ID(n,s)
!	Added "prior itime tracking" capability back to code since negative ID bug has been resolved
!8 Jan 2022
!	When counting number of fhours present at times during features' lifespan, it appears that "prior time tracked" features have a strange drop at feature life hour 6, but rebound at feature hour 12 and later. Unsurprisingly, counts from these "prior" tracks are less than counts from tracks initialized after the feature exists (e.g., there are more instances of "f048" at life hour 48 (48 hour fx init at life hour 0) than there are insances of "f048" at life hour 42 (48 hour fx init at life hour -6)
!	Error and distance information for prior tracks is incorrect...this might be associated with the above problem
!8 Jan 2022
! 	Several changes...
!	---Added an extra step to set any negative ID values to -1 (missing flag) since setting a feature between 1:s-1 to a negative value at ipass=npass would result in a negative ID being passed to the output file
!	---Originally, the history variables for "prior" features were being set the history of the verification data, which does not make sense since "prior" features by definition do not originate from a verification time. Set the first history step to be from the feature, itself.
!	---For "prior" features...Increased the maximum allowable penalty to PMAX instead of 0.5
!	---For "prior" features...Increased the maxumum allowable "DISTANCE" to 750.0 km from 350.0 km (relative to the verification point at feature life time 0). This should be okay, since the verification point acts as the predicted point.
!9 Jan 2022
!	Commented out all "runav_" variables. They are unused.
!	Corrected a problem with searching for IDs at previous verification times that caused the code to run incredibly slowly after two years. When checking for (a) an analysis history or (b) if the "prior" search was testing against a first occurance of a verification ID, the code would search ALL prior times (1:n). This is unneccessary, as, if less than 2 of ID(n-3:n) equal ID(n,t) or ID(n+f,t), then there is no analysis history AND it is the first insance of the ID. Changed count() to only search ID(max(1,n-5):n) or ID(max(1,n-5):n+f) for speed.
!10 Jan 2022
!	Corrected a new error that passed a missing DISTANCE into the overlap function in the "prior" forecast step
!	Since the changes on 9 Jan increased speed, returned the npass value to 4 passes.
!14 Jan 2022
!	Following conversation with Craig on 13 Jan...
!	Increased "prior" maximum predicted error to 1000km for consistency with other forecasts (was 750km)
!	Removed minimum overlap step from "prior"
!	Redirected output to version7
!9 Feb 2022
!	.dat files are now split between northern and southern hemispheres for efficiency. Tracking is now done seperately for each hemisphere. Some changes necessary to directory and file references in wrapper .csh script and .f90 code
!13-14 Feb 2022
!	Added flags to indicate times of max So and min Z500 (refer to the ANALYSIS TIME)
!	Disabled "prior" tracks
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Major Revisions. See track_forecast.f90.old for previous details!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!27 April 2022
!	Tracking concept revisited - 
!	Since forecasts are INDEPENDENT, can track forecast hours at different itimes independently
!	Rough outline of implementation plan - 
!	1. Track features from analysis time to first fx hour
!	2. Any features that do not correspond to a tracked f000 ID are then MATCHED to an f000 ID at the CURRENT VALID TIME
!	3. All feature IDs (including those that were matched) are tracked to the next forecast hour
!	4. Repeat matching of untracked IDs at CURRENT VALID TIME
!	5-n. Repeat steps through fhour 240
!	-Conceptually, this is a compromise between the matching scheme and the tracking scheme. It should be more efficient than the original tracking scheme, but also retain useful information (e.g., feature motion)
!3 May 2022
!	Features that are matched are compared to the previous valid time to attempt to find the effective feature velocity 
!	Added a flag for features tracked, matched, or neither
!4 May 2022
!	Once a track is "lost" through forecast space, no attempt is made to restart/match this track; could lead to features being "tracked" that aren't necessarily coherent in forecast space.
!	Added a limiter to prevent forecasts initialized near the end of the analysis period from attempting to match/track beyond the final analysis date
!	Fixed a bug that was causeing subsequently tracked, matched features to not retain the "current" duration of the feature they matched to
!4 May 2022
!	Prepared to batch submit all months separately. Should take ~5 minutes on casper to run all forecasts
!5 May 2022
!	Added auxilliary output file that holds the lat/lon of corresponding verification features in order to compute location-based position (and other) error characteristics. 
!	....For now, skip the auxilliary stream and simply add the corresponding valid lat and lon to the final 2 columns of the output file
!	...Adding valid lat, valid lon, valid So, valid Ro, Valid Zmin
!18 May 2022
!	Forbid matching if abs(Slat-Vlat).gt.30 -> cannot match a feature across the pole (can cause enormous errors if matched feature is tracked & is obviously not representative of the verification feature)
!	Decrease normalization factor for matching penalty to 750 km - Ensure the match is close to the valid time feature
!	Adjust weighting of kasuga, met, and distance penalty terms (1.0 met, 0.8 kasuga, 1.2 distance) - penalize the distance most harshly, kasuga least, and met middle
!	Some definitions - KPTY is the average of the normalized So, Ro, BGx, and BGy penalties (4 terms)
!			 - METPTY is the average of the normalized penalties for all met terms at each level (15 terms)
!			 - DPTY is the normalized distance penalty (1 term)
!	Use the history_dt variable to test the most recent 4 instances of an ID. if 2 or more of the most recent 4 history times are "12 hours", then 2 or more times were missed. Discontinue tracking this feature.
!30 June 2022
!	From analysis tracking, user inputs are now used to configure the "tracking" component of the scheme
!	Updated code to read fixed-radius variables from input files
! 	Updated matching scheme to use fixed-radius parameters with user-selectable normalization factors (currently 2*std from the 2020 test period over the globe)
!22 Sep 2022
!	Unexpected position errors can exceed 10000km 0_o Need to rerun tracking to prevent output fields of ******** in .track data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
use netcdf

implicit none

integer				:: nx, ny, nf, ns, nt, nsmax, nsfmax, s, f, t, tt, ss, n
integer				:: i, ii, j, k, x, y, r, a, b, np1, np2, nz1, np
integer				:: i0, i1, j0, j1, xx0, xx1
integer				:: istat,istat2,istat3
real				:: d, hgtmin, dmin,dthresh
real				:: haversine, c1c2_overlap
real				:: ST, FN
integer, allocatable		:: ID(:,:), nid(:,:), TrackMatch(:,:)
real, allocatable		:: d2prev(:,:),o2prev(:,:)
real, allocatable		:: So(:,:), Slat(:,:), Slon(:,:), Ro(:,:), SR(:,:), BGo(:,:), BGoy(:,:), BGox(:,:), Zmin(:,:), Zlat(:,:), Zlon(:,:)
real, allocatable		:: z850(:,:), z500(:,:), z200(:,:), t850(:,:), t500(:,:), t200(:,:), u850(:,:), u500(:,:), u200(:,:), v850(:,:), v500(:,:), v200(:,:), mr850(:,:), mr500(:,:), mr200(:,:)
real, allocatable		:: FIXEDz500(:,:), FIXEDt500(:,:), FIXEDu500(:,:),FIXEDv500(:,:),FIXEDmr500(:,:)

integer, allocatable		:: isFMaxSo(:,:),isFMinZ(:,:)

integer, allocatable		:: vID(:,:)
real, allocatable		:: vSo(:,:), vSlat(:,:), vSlon(:,:), vRo(:,:), vSR(:,:), vBGo(:,:), vBGoy(:,:), vBGox(:,:), vZmin(:,:), vZlat(:,:), vZlon(:,:)
real, allocatable		:: vz850(:,:), vz500(:,:), vz200(:,:), vt850(:,:), vt500(:,:), vt200(:,:), vu850(:,:), vu500(:,:), vu200(:,:), vv850(:,:), vv500(:,:), vv200(:,:), vmr850(:,:), vmr500(:,:), vmr200(:,:)
real, allocatable		:: vFIXEDz500(:,:), vFIXEDt500(:,:), vFIXEDu500(:,:),vFIXEDv500(:,:),vFIXEDmr500(:,:)

integer, allocatable		:: visFMaxSo(:,:),visFMinZ(:,:)
integer				:: holdSoFlag, holdZFlag

!real, allocatable		:: runav_So(:,:),runav_BG(:,:),runav_Ro(:,:),runav_la(:,:),runav_lo(:,:),runav_dx(:,:),runav_dy(:,:)
real, allocatable		:: hstry_So(:,:,:),hstry_BG(:,:,:),hstry_Ro(:,:,:), hstry_la(:,:,:),hstry_lo(:,:,:),hstry_dx(:,:,:),hstry_dy(:,:,:)

real, allocatable		:: hstry_dt(:,:,:)
real, allocatable		:: saveDUR(:,:),FERRX(:,:),FERRY(:,:),FERR(:,:),saveFTARG(:,:)

integer, parameter		:: npass = 4
integer				:: ipass
real				:: CSOMIN,TSOMIN,OMIN,PMAX,DMAX,EMAX,OXMAX,OYMAX,SNORM,BNORM,RNORM,BNORMXY
real				:: OVERLAP,DISTANCE
real				:: PDX,PDY,DX,DY,EDIST,OYRAT,OXRAT,pSlat,pSlon
real				:: SPTY,BPTY,RPTY,DPTY,PSUM
real				:: TARGET
real, allocatable		:: FTARGET(:,:)
integer				:: prior,validind,TESTID
real, allocatable		:: saveDX(:,:),saveDY(:,:),saveDIST(:,:),saveDT(:,:),saveVLat(:,:),saveVLon(:,:),saveVSo(:,:),saveVRo(:,:),saveVZmin(:,:)
real, allocatable		:: vsaveDX(:,:),vsaveDY(:,:),vsaveDIST(:,:),vsaveDT(:,:)
real, allocatable		:: vsaveDUR(:,:),vFERRX(:,:),vFERRY(:,:),vFERR(:,:),vsaveFTARG(:,:)
real				:: ERROR_DIST,ERROR_DY,ERROR_DX,VAL_LAT,VAL_LON,VAL_SO,VAL_RO,VAL_ZMIN
real, allocatable		:: vmaxID(:,:),maxID(:,:)

real				:: Z850NORM, Z500NORM, Z200NORM, Z850PTY, Z500PTY, Z200PTY
real				:: T850NORM, T500NORM, T200NORM, T850PTY, T500PTY, T200PTY
real				:: U850NORM, U500NORM, U200NORM, U850PTY, U500PTY, U200PTY
real				:: V850NORM, V500NORM, V200NORM, V850PTY, V500PTY, V200PTY
real				:: Q850NORM, Q500NORM, Q200NORM, Q850PTY, Q500PTY, Q200PTY
real				:: ZPTY, TPTY, UPTY, VPTY, QPTY, METPTY, METMAX
real				:: KPTY, KMAX, TOTPTY
real				:: BNORMXY_m,SNORM_m,PMAX_m,DMAX_m,RNORM_m

real				:: user_pmax,user_normso,user_normbgo,user_normro,user_emax,user_xoppmax,user_yoppmax,user_dmax,user_Mpmax,user_Mnorm_so,user_Mnorm_bgo,user_Mnorm_ro,user_Mnorm_dist,user_Mnorm_z500,user_Mnorm_t500,user_Mnorm_u500,user_Mnorm_v500,user_Mnorm_mr500


integer				:: nvalid, VER
character(len=10)		:: nvstr
character(len=80)		:: name
character(len=400)		:: ffilelist
character(len=400)		:: line, line2, line3
character(len=400)		:: vfile, listfile, infile, ffile
character(len=400)		:: iodir
character(len=10)		:: hemis
character(len=5)		:: hext
character(len=20),allocatable   :: fhour(:,:),vfhour(:,:)
character(len=400),allocatable	:: outfile(:)				!,auxoutfile(:)
character(len=450)		:: fstr
character(len=20),allocatable	:: itime(:,:),vitime(:,:)
character(len=400)		:: outdir
character(len=10)		:: st_pmax,st_normso,st_normbgo,st_normro,st_emax,st_xoppmax,st_yoppmax,st_dmax,st_Mpmax,st_Mnorm_so,st_Mnorm_bgo,st_Mnorm_ro,st_Mnorm_dist,st_Mnorm_z500,st_Mnorm_t500,st_Mnorm_u500,st_Mnorm_v500,st_Mnorm_mr500

!!!!!!Get some command line arguments
call get_command_argument(1,vfile)
call get_command_argument(2,ffile)
call get_command_argument(3,nvstr)
call get_command_argument(4,iodir)
call get_command_argument(5,st_pmax)
call get_command_argument(6,st_normso)
call get_command_argument(7,st_normbgo)
call get_command_argument(8,st_normro)
call get_command_argument(9,st_emax)
call get_command_argument(10,st_xoppmax)
call get_command_argument(11,st_yoppmax)
call get_command_argument(12,st_dmax)
call get_command_argument(13,st_Mpmax)
call get_command_argument(14,st_Mnorm_so)
call get_command_argument(15,st_Mnorm_bgo)
call get_command_argument(16,st_Mnorm_ro)
call get_command_argument(17,st_Mnorm_dist)
call get_command_argument(18,st_Mnorm_z500)
call get_command_argument(19,st_Mnorm_t500)
call get_command_argument(20,st_Mnorm_u500)
call get_command_argument(21,st_Mnorm_v500)
call get_command_argument(22,st_Mnorm_mr500)


read(nvstr , * ) nvalid
read(st_pmax    , * ) user_pmax
read(st_normso  , * ) user_normso
read(st_normbgo , * ) user_normbgo
read(st_normro  , * ) user_normro
read(st_emax    , * ) user_emax   
read(st_xoppmax , * ) user_xoppmax
read(st_yoppmax , * ) user_yoppmax
read(st_dmax    , * ) user_dmax   

read(st_Mpmax	     , * )user_Mpmax	
read(st_Mnorm_so    , * ) user_Mnorm_so	
read(st_Mnorm_bgo   , * ) user_Mnorm_bgo  
read(st_Mnorm_ro    , * ) user_Mnorm_ro	
read(st_Mnorm_dist  , * ) user_Mnorm_dist 
read(st_Mnorm_z500  , * ) user_Mnorm_z500 
read(st_Mnorm_t500  , * ) user_Mnorm_t500 
read(st_Mnorm_u500  , * ) user_Mnorm_u500 
read(st_Mnorm_v500  , * ) user_Mnorm_v500 
read(st_Mnorm_mr500 , * ) user_Mnorm_mr500


!!!!!Name the output text file
outdir = iodir
write(*,*) 'outdir=',trim(outdir)

!!!!!Get dimensions for forecast lists!!!!!
call cpu_time(ST)
nf    = 0
nsfmax = 0
open(unit = 10, file = trim(ffile), status='old', access='sequential', form='formatted', action='read')
do
  read(10,'(a256)',iostat=istat) infile   
  if(istat /= 0)exit
  ns=0
  open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do
    if(s.eq.0)then
      read(12,*,iostat=istat3) line2
    else
      read(12,*,iostat=istat3) line3
    end if    
    
    if(istat3 /=0)exit
      
    ns=ns+1
    nsfmax = max(ns,nsfmax)
  end do    
  close(12)  
  nf=nf+1
end do
close(10)
write(*,'(A,I4)')'nf=',nf

allocate(itime(nf,nsfmax),fhour(nf,nsfmax),ID(nf,nsfmax),TrackMatch(nf,nsfmax),So(nf,nsfmax),Slat(nf,nsfmax),Slon(nf,nsfmax),Ro(nf,nsfmax),SR(nf,nsfmax),BGo(nf,nsfmax),BGoy(nf,nsfmax),BGox(nf,nsfmax),Zmin(nf,nsfmax),Zlat(nf,nsfmax),Zlon(nf,nsfmax),d2prev(nf,nsfmax),o2prev(nf,nsfmax))
allocate(outfile(nf))
!allocate(auxoutfile(nf))
allocate(z850(nf,nsfmax), z500(nf,nsfmax), z200(nf,nsfmax), t850(nf,nsfmax), t500(nf,nsfmax), t200(nf,nsfmax), u850(nf,nsfmax), u500(nf,nsfmax), u200(nf,nsfmax), v850(nf,nsfmax), v500(nf,nsfmax), v200(nf,nsfmax), mr850(nf,nsfmax), mr500(nf,nsfmax), mr200(nf,nsfmax))
allocate(FIXEDz500(nf,nsfmax), FIXEDt500(nf,nsfmax), FIXEDu500(nf,nsfmax),FIXEDv500(nf,nsfmax),FIXEDmr500(nf,nsfmax))

allocate(nid(nf,nsfmax))
!allocate(runav_So(nf,nsfmax),runav_BG(nf,nsfmax),runav_Ro(nf,nsfmax))
allocate(FTARGET(nf,nsfmax))
allocate(saveDX(nf,nsfmax),saveDY(nf,nsfmax),saveDIST(nf,nsfmax),saveDT(nf,nsfmax),saveVLat(nf,nsfmax),saveVLon(nf,nsfmax),saveVSo(nf,nsfmax),saveVRo(nf,nsfmax),saveVZmin(nf,nsfmax))
allocate(hstry_So(nf,nsfmax,300),hstry_BG(nf,nsfmax,300),hstry_Ro(nf,nsfmax,300),hstry_dt(nf,nsfmax,300))
!allocate(runav_la(nf,nsfmax),runav_lo(nf,nsfmax),runav_dx(nf,nsfmax),runav_dy(nf,nsfmax))
allocate(hstry_la(nf,nsfmax,300),hstry_lo(nf,nsfmax,300),hstry_dx(nf,nsfmax,300),hstry_dy(nf,nsfmax,300))
allocate(isFMaxSo(nf,nsfmax),isFMinZ(nf,nsfmax))

allocate(saveDUR(nf,nsfmax))
allocate(saveFTARG(nf,nsfmax))
allocate(FERRX(nf,nsfmax),FERRY(nf,nsfmax),FERR(nf,nsfmax))
allocate(maxID(nf,nsfmax))

call cpu_time(FN)
print '("Time to determine fx array sizes = ",f6.3," seconds.")',FN-ST


!!!!!Get dimensions for verification lists!!!!!
call cpu_time(ST)
nsmax = 0
nt=0
open(unit = 11, file = trim(vfile), status='old', access='sequential',form='formatted', action='read')
do
  read(11,'(a256)',iostat=istat) infile
  if(istat /=0)exit
    
  ns=0
  open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do
    if(s.eq.0)then
      read(12,*,iostat=istat2) line2
    else
      read(12,*,iostat=istat2) line3
    end if        
    if(istat2 /=0)exit
    ns=ns+1
    nsmax = max(ns,nsmax)
  end do
    
  close(12)
  nt=nt+1
end do
close(11)  
write(*,*)'ns,nt=',ns,nt

allocate(vitime(nt,nsmax),vfhour(nt,nsmax),vID(nt,nsmax),vSo(nt,nsmax),vSlat(nt,nsmax),vSlon(nt,nsmax),vRo(nt,nsmax),vSR(nt,nsmax),vBGo(nt,nsmax),vBGoy(nt,nsmax),vBGox(nt,nsmax),vZmin(nt,nsmax),vZlat(nt,nsmax),vZlon(nt,nsmax))
allocate(vsaveDX(nt,nsmax),vsaveDY(nt,nsmax),vsaveDIST(nt,nsmax),vsaveDT(nt,nsmax))
allocate(vz850(nt,nsmax), vz500(nt,nsmax), vz200(nt,nsmax), vt850(nt,nsmax), vt500(nt,nsmax), vt200(nt,nsmax), vu850(nt,nsmax), vu500(nt,nsmax), vu200(nt,nsmax), vv850(nt,nsmax), vv500(nt,nsmax), vv200(nt,nsmax), vmr850(nt,nsmax), vmr500(nt,nsmax), vmr200(nt,nsmax))
allocate(vFIXEDz500(nt,nsmax), vFIXEDt500(nt,nsmax), vFIXEDu500(nt,nsmax),vFIXEDv500(nt,nsmax),vFIXEDmr500(nt,nsmax))

allocate(visFMaxSo(nt,nsfmax),visFMinZ(nt,nsfmax))
allocate(vsaveDUR(nt,nsmax))
allocate(vsaveFTARG(nt,nsmax))
allocate(vFERRX(nt,nsmax),vFERRY(nt,nsmax),vFERR(nt,nsmax))
allocate(vmaxID(nt,nsmax))
call cpu_time(FN)
print '("Time to determine analysis array sizes = ",f8.3," seconds.")',FN-ST







!!!!!Set dummy values for lists!!!!!
itime = "-9999"
fhour = "-9999"
ID = -1
TrackMatch = -1
isFMaxSo = -1	! Note this and the isFMinZ flags refer to the VERIFICATION feature
isFMinZ = -1
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
FERRY = -99999.9
FERRX = -99999.9
FERR = -99999.9
maxID = -9999.9
saveVLat = -9999.9
saveVLon = -9999.9
saveVSo = -9999.9
saveVRo = -9999.9
saveVZmin = -9999.9

vitime = "-9999"
vfhour = "-9999"
vID = -1
vSo = -9999.9
vSlat = -9999.9
vSlon = -9999.9
vRo = -9999.9
vSR = -9999.9
vBGo = -9999.9
vBGoy = -9999.9
vBGox = -9999.9
vZmin = -9999.9
vZlat = -9999.9
vZlon = -9999.9
vsaveDX = -9999.9
vsaveDY = -9999.9
vsaveDIST = -9999.9
vsaveDT = 0.0
vsaveDUR = 0.0
vsaveFTARG = 9999.9
!vFTARGET = 9999.9
vFERRY = -99999.9
vFERRX = -99999.9
vFERR = -99999.9




!!!!!Read Verif vars from .track files!!!!!
call cpu_time(ST)
open(unit = 11, file = trim(vfile), status='old', access='sequential', form='formatted', action='read')
do t=1,nt
  read(11,'(a256)',iostat=istat) infile     
  if(istat /= 0)then
      write(*,*)'read ',trim(infile), 'istat=',istat
      exit
  endif
  
  open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do s=0,nsmax
    if(s.eq.0)then
      read(12,*,iostat=istat2) line
    else
      read(12,*,iostat=istat2) vitime(t,s), vfhour(t,s), vID(t,s), vSo(t,s), vSlat(t,s), vSlon(t,s), visFMaxSo(t,s), vRo(t,s), vSR(t,s), vBGo(t,s), vBGoy(t,s), vBGox(t,s), vZmin(t,s), vZlat(t,s), vZlon(t,s), visFMinZ(t,s), vz850(t,s), vz500(t,s), vz200(t,s), vt850(t,s), vt500(t,s), vt200(t,s), vu850(t,s), vu500(t,s), vu200(t,s), vv850(t,s), vv500(t,s), vv200(t,s), vmr850(t,s), vmr500(t,s), vmr200(t,s), vFIXEDz500(t,s), vFIXEDt500(t,s), vFIXEDu500(t,s), vFIXEDv500(t,s), vFIXEDmr500(t,s), vsaveDY(t,s), vsaveDX(t,s), vsaveDIST(t,s), vsaveDT(t,s), vsaveDUR(t,s),vmaxID(t,s), vsaveFTARG(t,s),vFERRY(t,s),vFERRX(t,s),vFERR(t,s)
    end if    
    if(istat2 /=0)exit
  end do
  close(12)
end do
close(11)



call cpu_time(FN)
print '("Time to gather analysis data = ",f6.3," seconds.")',FN-ST

!!!!!Read Forecast vars from .dat files!!!!!
call cpu_time(ST)
open(unit = 10, file = trim(ffile), status='old', access='sequential', form='formatted', action='read')
do f=1,nf
  read(10,'(a256)',iostat=istat) infile      
  if(istat /= 0)then
      write(*,*)'read ',trim(infile), 'istat=',istat
      exit
  endif

  open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  write(*,'(A,I5)')'nsfmax=',nsfmax
  do s=0,nsfmax
    if(s.eq.0)then
      read(12,*,iostat=istat3) line
    else
      read(12,*,iostat=istat3) itime(f,s), fhour(f,s), So(f,s), Slat(f,s), Slon(f,s), Ro(f,s), SR(f,s), BGo(f,s), BGoy(f,s), BGox(f,s), Zmin(f,s), Zlat(f,s), Zlon(f,s), z850(f,s), z500(f,s), z200(f,s), t850(f,s), t500(f,s), t200(f,s), u850(f,s), u500(f,s), u200(f,s), v850(f,s), v500(f,s), v200(f,s), mr850(f,s), mr500(f,s), mr200(f,s), FIXEDz500(f,s), FIXEDt500(f,s), FIXEDu500(f,s), FIXEDv500(f,s), FIXEDmr500(f,s)
    end if      
    if(istat3 /= 0)then
        write(*,*)'read ',trim(infile), ' s=',s,' istat3=',istat3
        exit
    endif
  end do  
  close(12)
  outfile(f) = trim(infile(1:len_trim(infile)-3)) // "track"
  !auxoutfile(f) = trim(outdir) // trim(infile(57:80) // ".track.error" 
end do

close(10)
call cpu_time(FN)
print '("Time gather forecast data = ",f8.3," seconds.")',FN-ST









!!!!!Find forecast points that are within threshold distance of analysis points!!!!!
!!!!!Do not create "new features" in the forecast, only propagate identified features from the analysis forward through the forecast
!!!!!Match current forecast hour (f) features (s) to previous time (f-1 or "v") features (t)!!!!!
!!!!!If f=1, compare forecast features to the analysis time!!!!!
!!!!!If f>1, compare forecast features to the previous forecast time!!!!!
!!!!!If a match is made, ID is carried forward to feature (f,n,s).
!!!!!!!!!!FAIL TEST - Optimal radii of features (f,n,s) and (f-1,n,t) do not overlap
!!!!!!!!!!FAIL TEST - Distance between features exceeds 1500 km
!!!!!!!!!!If number of occurances of ID is at least 2, use position history to estimate the expected position of the feature at time (n)
!!!!!!!!!!!!!!!FAIL TEST - Position of feature (f,n,s) is more than 1000km different than the expected position
!!!!!!!!!!!!!!!FAIL TEST - delta-x between feature (f,n,s) and feature (f-1,n,t) is greater than 300 km in the opposite direction of the expected position
!!!!!!!!!!!!!!!FAIL TEST - delta-y between feature (f,n,s) and feature (f-1,n,t) is greater than 600 km in the opposite direction of the expected position
!!!!!!!!!!PENALTY - Large change in So between (f-1,n,t) and (f,n,s) ---- PTY = abs(So(f,n,s)-So(f-1,n,t))/10.0
!!!!!!!!!!PENALTY - Large change in BGo between (f-1,n,t) and (f,n,s) ---- PTY = abs(BGo(f,n,s)-BGo(f-1,n,t))/10.0
!!!!!!!!!!PENALTY - Large change in Ro between (f-1,n,t) and (f,n,s) ---- PTY = abs(Ro(f,n,s)-Ro(f-1,n,t))/900.0
!!!!!!!!!!FAIL TEST - Total penalty term is greater than 1.5
!!!!!!!!!!If multiple features (n-1,t) pass all tests - Keep the one with the smallest value of PENALTY-%OVERLAP
!!!!!!If the same (n-1,t) feature is a best match to multiple (n,s) features, keep the match with feature with the longest history
!!!!!!!!!!If the same history, keep the match with the smallest value of PENALTY-%OVERLAP
!!!!!!!!!!Reassign the "unmatched" feature (n,s) an ID of -1*ID(n-1,t) so that it may be tested on subsequent pass to determine if it is an acceptable match for a different point 
!!!!!!If a feature (n,s) has no match from (n-1,t), repeat the above examination for features at (n-2,t) that were not propagated to time (n-1)
call cpu_time(ST)




!!!!!Set some failure thresholds
CSOMIN = 5.0		! Minimum So for a cutoff to be identified
TSOMIN = 5.0		! Minimum So for a trough to be identified
OMIN   = 0.0		! Minimum overlap percentage for matching features
PMAX   = user_pmax	! 1.5		! Maximum allowable penalty points before a feature is disqualified
DMAX   = user_dmax	! 1500.0		! Maximum distance between matches allowed
EMAX   = user_emax	! 1500.0		! Maximum "error" distance between predicted and candidate points
OXMAX  = user_xoppmax	! 1000.0       ! Maximum distance a matching point can be in the opposite x-direction of the predicted point (revised from 300km to 1000km since this seems to have caused problems with the september 2020 colorado cutoff track
OYMAX  = user_yoppmax	! 1000.0       ! ".."y-direction"" More lenient to allow for features moving around larger scale troughs to move to the north (revised from 600, same as above)
SNORM  = user_normso	! 15.0         ! Normalization for the So penalty ! Motivated by Sep 2020 case, changed to 15 from 10)
BNORM  = user_normbgo	! 40.0         ! Normalization for the BGo penalty     ! motivated by Sep 2020 case, changed to 35 from 15)
RNORM  = user_normro	! 1000.0       ! Normalization for the Ro penalty      ! motivated by the Sep 2020 case, changed to 1000 from 900

!!!!!!!thresholds for matching scheme

!!! Pre june 2022
!Z850NORM = 75.0		!100.0
!Z500NORM = 150.0 	!200.0
!Z200NORM = 225.0 	!300.0

!T850NORM = 7.5		!10.0
!T500NORM = 7.5		!10.0
!T200NORM = 7.5		!10.0

!U850NORM = 10.0
!U500NORM = 10.0
!U200NORM = 10.0

!V850NORM = 5.0		! All were 10.0
!V500NORM = 5.0
!V200NORM = 5.0

!Q850NORM = 1.5
!Q500NORM = 0.75
!Q200NORM = 0.015

!BNORMXY_m = 6.0
!SNORM_m = 6.0
!RNORM_m = 750.0
!PMAX_m = 1.0
!DMAX_m = 750.0
!METMAX = 1.5

Z500NORM = user_Mnorm_z500
T500NORM = user_Mnorm_t500
U500NORM = user_Mnorm_u500
V500NORM = user_Mnorm_v500
Q500NORM = user_Mnorm_mr500
BNORMXY_m = user_Mnorm_bgo
SNORM_m = user_Mnorm_so
RNORM_m = user_Mnorm_ro
DMAX_m = user_Mnorm_dist
PMAX_m = user_Mpmax



do f=2,nf
  if(f.gt.nt) exit															! Used only during the final month of the dataset. Forecasts initialized within the last 10 days of the period of interest verify after the end of the analysis data. To prevent unexpected nonsensical tracks/matches of these features, terminate all tracking once the forecast hour exceeds the maximum analysis hour
  do ipass=1,npass															! START DO NUM PASSES  
    do s=1,nsfmax														      	! START DO NUM FEATURES s at FHOUR f
      if(itime(f,s).eq."-9999") exit  											      		! IF AT A "MISSING" itime, exit the s loop (no more features in the file)     
      if(f.eq.2) then														      	! TEST IF at F006. If f=2, we are comparing f006 to the analysis time (analysis features identified in track_kasugaEA21 code)	   

        TARGET = 9999.9 													      	! SET INITIAL VALUE/THRESHODL FOR PENALTY-OVERLAP TERM
	do t=1,nsmax														      	! START DO NUM FEATURES t AT TIME n, FHOUR 0 (the verification/initial time data)
	  if(.not.(any(vID(min(f,nt):min(f+1,nt),:).eq.ID(f,t)))) cycle							     	 	! DO NOT continue tracking an ID that no longer exists in the verification dataset
          if(vsaveDUR(f-1,t).eq.vmaxID(f-1,t)) cycle											! DO NOT continue tracking beyond the end of an analysis time features
	  if(vitime(f-1,t).eq."-9999") exit											      	! IF AT A "MISSING" itime in the verification data, exit the t loop (no more features in the file)
          if(vID(f-1,t).lt.0) cycle												      	! IF THE VALID ID IS FLAGGED AS MISSING, SKIP IT
          if(-1*vID(f-1,t).eq.ID(f,s)) cycle											      	! IF THE CURRENT ID is NEGATIVE & EQUAL TO THE NEGATIVE OF THE VALID ID, SKIP IT (forbidden ID)
        		
          DISTANCE = haversine(Slat(f,s),Slon(f,s),vSlat(f-1,t),vSlon(f-1,t))							      	! HOW FAR DID THE FORECAST FEATURE TRAVEL FROM THE FEATURE AT F000?
          DX = haversine(vSlat(f-1,t),Slon(f,s),vSlat(f-1,t),vSlon(f-1,t))  							      	! "..." IN THE X-DIR?
          DY = haversine(Slat(f,s),vSlon(f-1,t),vSlat(f-1,t),vSlon(f-1,t))  							      	! "..." IN THE Y-DIR?
          if((Slat(f,s)-vSlat(f-1,t)).lt.0) DY = -1.0*DY									      	! GET THE SIGNS OF THE X AND Y MOTION
          if(((Slon(f,s)-vSlon(f-1,t)).lt.0).and.((Slon(f,s)-vSlon(f-1,t)).gt.-100.0)) DX = -1.0*DX
          if(((Slon(f,s)-vSlon(f-1,t)).gt.0).and.((Slon(f,s)-vSlon(f-1,t)).gt.100.0))  DX = -1.0*DX     
        
          if(DISTANCE.gt.DMAX) cycle												      	! IF THE DISTANCE IS LARGER THAN THE MAXIMUM ALLOWABLE DISTANCE, SKIP THIS FEATURE
        
          OVERLAP = c1c2_overlap(DISTANCE,Ro(f,s),vRo(f-1,t))									      	! COMPUTE THE PERCENT OVERLAP OF THE FEATURE AT F000 and F006
          if(OVERLAP.le.OMIN) cycle												      	! IF THE OVERLAP IS LESS THAN or EQUAL TO THE MINUMUM OVERLAP (0.0), SKIP THIS FEATURE
       	
        	
          !if(count(vID(max(1,n-5):n,:).eq.vID(n,t)).gt.1) then  								        ! TEST IF THE FEATURE IN THE ANALYSIS DATA HAS EXISTED FOR MORE THAN 1 TIMESTEP 
          !															      	! .....Do we have an analysis history to use to help get the first forecast feature? (Is the # vID(n,t) from 1->n gt 1?)
          !																	  
          !  if(any(vID(n-1,:).eq.vID(n,t))) then										      	! DETERMINE IF THE MOST RECENT ANALYSIS HISTORY SHOULD USE TIME n-1 or n-2
          !    prior = n-1
          !  else
          !    prior = n-2
          !  end if														      ! END IF DETERMINING PRIOR
          !  
          !  do ss=1,nsmax													      ! START DO NUM FEATURES ss AT TIME prior TO GET THE HISTORY INDEX
          !    if(vID(prior,ss).ne.vID(n,t)) cycle										      ! IF THIS IS NOT THE FEATURE WE ARE INTERESTED IN, SKIP IT		    
          !    pSlat = vSlat(prior,ss)
          !    pSlon = vSlon(prior,ss)
          !  end do														      ! END DO NUM FEATURES ss
          ! 
          !  PDX       = haversine(pSlat,pSlon,pSlat,vSlon(n,t)) 								      ! ASSUME PERSISTENCE AND GET THE PREDICTED DX AND DY FROM THE FEATURE HISTORY
          !  PDY       = haversine(pSlat,pSlon,vSlat(n,t),pSlon) 								      ! GET THE SIGNS OF THE PREDICTED X AND Y MOTION
          !  if((vSlat(n,t)-pSlat).lt.0) PDY = -1.0*PDY  	
          !  if(((vSlon(n,t)-pSlon).lt.0).and.((vSlon(n,t)-pSlon).gt.-100.0)) PDX = -1.0*PDX
          !  if(((vSlon(n,t)-pSlon).gt.0).and.((vSlon(n,t)-pSlon).gt.100.0))  PDX = -1.0*PDX													 
         ! 
         !   EDIST     = ((PDX-DX)**2+(PDY-DY)**2)**(0.5)									      ! GET THE POSTITION ERROR BETWEEN THE PREDICTED POINT AND THE TEST POINT
         !   OXRAT     = DX/PDX  												      ! DETERMINE IF THE TEST POINT MOVED IN THE OPOSITE DIRECTION OF THE PREDICTED POINT
         !   OYRAT     = DY/PDY
          
          
            
          !  if(OXRAT.lt.0.0.and.abs(PDX-DX).gt.OXMAX) cycle									      ! Skip this (t) because it moved too far in the opposite of the expected x-direction
          !  if(OYRAT.lt.0.0.and.abs(PDY-DY).gt.OYMAX) cycle									      ! Skip this (t) because it moved too far in the opposite of the expected y-direction
          !  if(EDIST.gt.EMAX) cycle												      ! Skip this (t) because it is too far from its expected position
        			  
         ! end if														      ! END IF DETERMINING AN ANALYSIS HISTORY
        
          validind = findloc(vID(f,:),vID(f-1,t),1)										      ! GET ACTUAL ERROR INFORMATION FOR THE TEST POINT s at the valid time corresponding to the fhour assoicted with test point s
          if(validind.gt.0) then
            ERROR_DIST = haversine(Slat(f,s),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
            ERROR_DX = haversine(vSlat(f,validind),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
            ERROR_DY = haversine(Slat(f,s),vSlon(f,validind),vSlat(f,validind),vSlon(f,validind))
            
	    VAL_LAT = vSlat(f,validind)
	    VAL_LON = vSlon(f,validind)
	    VAL_SO = vSo(f,validind)
            VAL_RO = vRo(f,validind)
	    VAL_ZMIN = vZmin(f,validind)	
	        
            if((Slat(f,s)-vSlat(f,validind)).lt.0) ERROR_DY = -1.0*ERROR_DY 						      		! GET THE SIGN OF THE FORECAST ERROR    
            if(((Slon(f,s)-vSlon(f,validind)).lt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.-100.0)) ERROR_DX = -1.0*ERROR_DX
            if(((Slon(f,s)-vSlon(f,validind)).gt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.100.0))  ERROR_DX = -1.0*ERROR_DX
            
            holdSoFlag = visFMaxSo(f,validind)										      		! Save whether the verification feature is at its MaxSo or MinZ time
            holdZFlag = visFMinZ(f,validind)
          else
            ERROR_DIST = -99999.9
	    ERROR_DX = -99999.9
            ERROR_DY = -99999.9
            
	    VAL_LAT = -9999.9
	    VAL_LON = -9999.9
	    VAL_SO = -9999.9
            VAL_RO = -9999.9
	    VAL_ZMIN = -9999.9 
	    
            holdSoFlag = -1
            holdZFlag = -1
          end if														      ! END IF COMPARISON TO THE ANALYSIS
        
          SPTY = abs(So(f,s)-vSo(f-1,t))/SNORM  										      ! DETERMINE PENALTIES BASED ON FEATURE CHARACTERISTICS
          BPTY = abs(BGo(f,s)-vBGo(f-1,t))/BNORM
          RPTY = abs(Ro(f,s)-vRo(f-1,t))/RNORM
          PSUM = SPTY+BPTY+RPTY
          !if((PSUM-OVERLAP).gt.PMAX) cycle												      ! IF THE PENALTY IS TOO LARGE, SKIP THIS FEATURE  	
          if(PSUM.gt.PMAX) cycle
          
	  !if((PSUM-OVERLAP).lt.TARGET) then											      ! IF PENALTY-OVERLAP IS LESS THAN THE TARGET(default or another feature), THIS (n,t) feature is the current best match to feature(f,n,s)
          if(PSUM.lt.TARGET) then
	    ID(f,s)		      	= vID(f-1,t)
            nid(f,s)  	     		= 2
            TrackMatch(f,s)		= 0
            !hstry_So(f,s,1)	    	= vSo(f-1,t)  										    ! START A HISTORY for FEATURE (f,n,s)
            !hstry_So(f,s,2)	    	= So(f,s)
            !hstry_BG(f,s,1)	    	= vBGo(f-1,t)
            !hstry_BG(f,s,2)	    	= BGo(f,s)
            !hstry_Ro(f,s,1)	    	= vRo(f-1,t)
            !hstry_Ro(f,s,2)	    	= Ro(f,s)
          
            hstry_la(f,s,1)	    	= vSlat(f-1,t)
            hstry_la(f,s,2)	    	= Slat(f,s)
            hstry_lo(f,s,1)	    	= vSlon(f-1,t)
            hstry_lo(f,s,2)	    	= Slon(f,s)
            hstry_dx(f,s,2)	    	= DX
            hstry_dy(f,s,2)	    	= DY
          
            !runav_So(f,n,s)		= sum(hstry_So(f,s,:))/2.0 
            !runav_BG(f,n,s)		= sum(hstry_BG(f,s,:))/2.0 
            !runav_Ro(f,n,s)	      	= sum(hstry_Ro(f,s,:))/2.0 

	    !runav_la(f,n,s)	      	= sum(hstry_la(f,s,:))/2.0
	    !runav_lo(f,n,s)	      	= sum(hstry_lo(f,s,:))/2.0
	    !runav_dx(f,n,s)	      	= hstry_dx(f,s,2)
	    !runav_dy(f,n,s)	      	= hstry_dy(f,s,2)
          
            saveDX(f,s)	      		= DX
            saveDY(f,s)	      		= DY
            saveDIST(f,s)	      	= DISTANCE
            saveDT(f,s)	      		= 6.0
            saveVLat(f,s)		= VAL_LAT
	    saveVLon(f,s)		= VAL_LON
	    saveVSo(f,s)		= VAL_SO
	    saveVRo(f,s)		= VAL_RO
	    saveVZmin(f,s)		= VAL_ZMIN
	    
            hstry_dt(f,s,1)	      	= vsaveDUR(f-1,t)	  
            hstry_dt(f,s,2)	      	= saveDT(f,s)
	    
	    saveDUR(f,s)	      	= sum(hstry_dt(f,s,1:nid(f,s)))
	    TARGET		      	= PSUM	!-OVERLAP
            FTARGET(f,s)	      	= TARGET
            saveFTARG(f,s)	      	= TARGET

            FERRY(f,s)	      		= ERROR_DY
            FERRX(f,s)	      		= ERROR_DX
            FERR(f,s) 	      		= ERROR_DIST
            isFMaxSo(f,s)	      	= holdSoFlag
            isFMinZ(f,s)	      	= holdZFlag
            maxID(f,s)	     		= vmaxID(f-1,t)
            
          end if														      ! END TARGET if()
        end do  														      ! END DO NUM FEATURES t do()

      else															      ! if f.gt.1, we are comparing fxxx to fxxx-previous (tracking between forecast hours)

        prior = f-1
        TARGET = 9999.9
        
        do t=1,nsfmax
          if(.not.(any(vID(min(f,nt):min(f+1,nt),:).eq.ID(prior,t)))) cycle						      ! DO NOT continue tracking an ID that no longer exists in the verification dataset
	  if(saveDUR(prior,t).eq.maxID(prior,t)) cycle										! DO NOT continue tracking beyond the end of an analysis time features
          if(itime(f,t).eq."-9999") exit
          if(ID(prior,t).lt.0) cycle
          if(-1*ID(prior,t).eq.ID(f,s)) cycle
          
	  if(count(hstry_dt(prior,t,(max(1,nid(prior,t)-3)):nid(prior,t)).eq.12).gt.1) cycle					! No longer tracking a feature if at least 2 of the last 4 times were skipped
	  
          DISTANCE = haversine(Slat(f,s),Slon(f,s),Slat(prior,t),Slon(prior,t))
          DX = haversine(Slat(prior,t),Slon(f,s),Slat(prior,t),Slon(prior,t))
          DY = haversine(Slat(f,s),Slon(prior,t),Slat(prior,t),Slon(prior,t))	
        											      ! Get the signs of the test point motion
          if((Slat(f,s)-hstry_la(prior,t,nid(prior,t))).lt.0) DY = -1.0*DY	
          if(((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).lt.0).and.((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.-100.0)) DX = -1.0*DX
          if(((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.0).and.((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.100.0))  DX = -1.0*DX   
          if(DISTANCE.gt.DMAX) cycle
        
          OVERLAP = c1c2_overlap(DISTANCE,Ro(f,s),Ro(prior,t))
          if(OVERLAP.le.OMIN) cycle
          
          if(nid(prior,t).gt.1) then
            PDX       = abs(hstry_dx(prior,t,nid(prior,t)))
            PDY       = abs(hstry_dy(prior,t,nid(prior,t))) 	  
        											      ! Get the signs of the predicted motion
            if((hstry_la(prior,t,nid(prior,t))-hstry_la(prior,t,nid(prior,t)-1)).lt.0) PDY = -1.0*PDY	
            if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).lt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.-100.0)) PDX = -1.0*PDX
            if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.100.0))  PDX = -1.0*PDX
        														  
            EDIST     = ((PDX-DX)**2+(PDY-DY)**2)**(0.5)
            OXRAT     = DX/PDX
            OYRAT     = DY/PDY
              
            if(OXRAT.lt.0.0.and.abs(PDX-DX).gt.OXMAX) cycle						      ! Skip this (t) because it moved too far in the opposite of the expected x-direction
            if(OYRAT.lt.0.0.and.abs(PDY-DY).gt.OYMAX) cycle						      ! Skip this (t) because it moved too far in the opposite of the expected y-direction
            if(EDIST.gt.EMAX) cycle								      ! Skip this (t) because it is too far from its expected position
          end if
         
          validind = findloc(vID(f,:),ID(prior,t),1)
          
          if(validind.gt.0) then
            ERROR_DIST = haversine(Slat(f,s),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
            ERROR_DX = haversine(vSlat(f,validind),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
            ERROR_DY = haversine(Slat(f,s),vSlon(f,validind),vSlat(f,validind),vSlon(f,validind))
            
	    VAL_LAT = vSlat(f,validind)
	    VAL_LON = vSlon(f,validind)	  
	    VAL_SO = vSo(f,validind)
            VAL_RO = vRo(f,validind)
	    VAL_ZMIN = vZmin(f,validind)
	    
            if((Slat(f,s)-vSlat(f,validind)).lt.0) ERROR_DY = -1.0*ERROR_DY       
            if(((Slon(f,s)-vSlon(f,validind)).lt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.-100.0)) ERROR_DX = -1.0*ERROR_DX
            if(((Slon(f,s)-vSlon(f,validind)).gt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.100.0))  ERROR_DX = -1.0*ERROR_DX
            
            holdSoFlag = visFMaxSo(f,validind)										      ! Save whether the verification feature is at its MaxSo or MinZ time
            holdZFlag = visFMinZ(f,validind)
          else
            ERROR_DIST = -99999.9
	    ERROR_DX = -99999.9
            ERROR_DY = -99999.9
            
	    VAL_LAT = -9999.9
	    VAL_LON = -9999.9
	    VAL_SO = -9999.9
            VAL_RO = -9999.9
	    VAL_ZMIN = -9999.9
	    
            holdSoFlag = -1										      
            holdZFlag = -1
          end if 
          SPTY = abs(So(f,s)-So(prior,t))/SNORM
          BPTY = abs(BGo(f,s)-BGo(prior,t))/BNORM
          RPTY = abs(Ro(f,s)-Ro(prior,t))/RNORM
          PSUM = SPTY+BPTY+RPTY
          
          !if((PSUM-OVERLAP).gt.PMAX) cycle
          if(PSUM.gt.PMAX) cycle
          				    							       ! Compare the combined penalty & overlap term to the previous "TARGET" combined term to beat
	  !if((PSUM-OVERLAP).lt.TARGET) then
	  if(PSUM.lt.TARGET) then
            ID(f,s) = ID(prior,t)								      ! Assign a new history and ID to feature (n,s)
            nid(f,s)        = nid(prior,t) + 1
            TrackMatch(f,s)		= 0
	    										      ! Assign characteristic history
	    !hstry_So(f,s,:)		= hstry_So(prior,t,:)
	    !hstry_BG(f,s,:)		= hstry_BG(prior,t,:)
	    !hstry_Ro(f,s,:)		= hstry_Ro(prior,t,:)
            !hstry_So(f,s,nid(f,s))     = So(f,s)
	    !hstry_BG(f,s,nid(f,s))     = BGo(f,s)
	    !hstry_Ro(f,s,nid(f,s))     = Ro(f,s)
        											      ! Assign position history
            hstry_la(f,s,:)		= hstry_la(prior,t,:)
	    hstry_lo(f,s,:)		= hstry_lo(prior,t,:)
	    hstry_dx(f,s,:)		= hstry_dx(prior,t,:)
            hstry_dy(f,s,:)		= hstry_dy(prior,t,:)	  
            hstry_la(f,s,nid(f,s))      = Slat(f,s)
	    hstry_lo(f,s,nid(f,s))      = Slon(f,s)
	    hstry_dx(f,s,nid(f,s))      = DX
            hstry_dy(f,s,nid(f,s))      = DY
        	
            !runav_So(f,s)	      	= sum(hstry_So(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 
            !runav_BG(f,s)	      	= sum(hstry_BG(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 
            !runav_Ro(f,s)	      	= sum(hstry_Ro(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 

	    !runav_la(f,s)	      	= sum(hstry_la(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/float(min(nid(f,s),4))
	    !runav_lo(f,s)	      	= sum(hstry_lo(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/float(min(nid(f,s),4))
	    !runav_dx(f,s)	      	= sum(hstry_dx(f,s,(max(2,nid(f,s)-2)):nid(f,s)))/float(min(nid(f,s),3))
	    !runav_dy(f,s)	      	= sum(hstry_dy(f,s,(max(2,nid(f,s)-2)):nid(f,s)))/float(min(nid(f,s),3))
          
            saveDX(f,s)	      		= DX
            saveDY(f,s)	      		= DY
            saveDIST(f,s)	      	= DISTANCE
            saveDT(f,s)	      		= 6.0
	    saveVLat(f,s)		= VAL_LAT
	    saveVLon(f,s)		= VAL_LON
	    saveVSo(f,s)		= VAL_SO
	    saveVRo(f,s)		= VAL_RO
	    saveVZmin(f,s)		= VAL_ZMIN
	    
            hstry_dt(f,s,:)	     	= hstry_dt(prior,t,:)
            hstry_dt(f,s,nid(f,s))    	= saveDT(f,s)
            saveDUR(f,s)	      	= sum(hstry_dt(f,s,1:nid(f,s)))
            TARGET		      	= PSUM !-OVERLAP
            FTARGET(f,s)	      	= TARGET
            saveFTARG(f,s)	      	= TARGET
          
            FERRY(f,s)	      		= ERROR_DY
            FERRX(f,s)	      		= ERROR_DX
            FERR(f,s) 	      		= ERROR_DIST
            isFMaxSo(f,s)	      	= holdSoFlag
            isFMinZ(f,s)	      	= holdZFlag
            
            maxID(f,s)	      		= maxID(prior,t)
          
          end if !TARGET if()
        end do ! t do()

              
        if(f.gt.3.and.ID(f,s).lt.0) then
          prior = f-2
          TARGET = 9999.9
          
          do t=1,nsfmax 
            if(.not.(any(vID(min(f,nt):min(f+1,nt),:).eq.ID(prior,t)))) cycle						      ! DO NOT continue tracking an ID that no longer exists in the verification dataset
	    if(saveDUR(prior,t).eq.maxID(prior,t)) cycle									! DO NOT continue tracking beyond the end of an analysis time features
	    if(itime(f,t).eq."-9999") exit
            if(ID(prior,t).lt.0) cycle
            if(-1*ID(prior,t).eq.ID(f,s)) cycle
            
	    if(count(hstry_dt(prior,t,(max(1,nid(prior,t)-3)):nid(prior,t)).eq.12).gt.1) cycle					! No longer tracking a feature if at least 2 of the last 4 times were skipped

	    
            DISTANCE = haversine(Slat(f,s),Slon(f,s),Slat(prior,t),Slon(prior,t))
            DX = haversine(Slat(prior,t),Slon(f,s),Slat(prior,t),Slon(prior,t))
            DY = haversine(Slat(f,s),Slon(prior,t),Slat(prior,t),Slon(prior,t)) 		  
        												      ! Get the signs of the test point motion
            if((Slat(f,s)-hstry_la(prior,t,nid(prior,t))).lt.0) DY = -1.0*DY	
            if(((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).lt.0).and.((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.-100.0)) DX = -1.0*DX
            if(((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.0).and.((Slon(f,s)-hstry_lo(prior,t,nid(prior,t))).gt.100.0))  DX = -1.0*DX
            
            if(DISTANCE.gt.DMAX) cycle
        
            OVERLAP = c1c2_overlap(DISTANCE,Ro(f,s),Ro(prior,t))
            if(OVERLAP.le.OMIN) cycle
        
            if(nid(prior,t).gt.1) then
              PDX     = abs(hstry_dx(prior,t,nid(prior,t)))
              PDY     = abs(hstry_dy(prior,t,nid(prior,t))) 	  
        											      ! Get the signs of the predicted motion
              if((hstry_la(prior,t,nid(prior,t))-hstry_la(prior,t,nid(prior,t)-1)).lt.0) PDY = -1.0*PDY 	
              if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).lt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.-100.0)) PDX = -1.0*PDX
              if(((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.0).and.((hstry_lo(prior,t,nid(prior,t))-hstry_lo(prior,t,nid(prior,t)-1)).gt.100.0))  PDX = -1.0*PDX
        										      
          
              EDIST   = ((PDX-DX)**2+(PDY-DY)**2)**(0.5)
              OXRAT   = DX/PDX
              OYRAT   = DY/PDY
        	
              if(OXRAT.lt.0.0.and.abs(PDX-DX).gt.OXMAX) cycle						      ! Skip this (t) because it moved too far in the opposite of the expected x-direction
              if(OYRAT.lt.0.0.and.abs(PDY-DY).gt.OYMAX) cycle						      ! Skip this (t) because it moved too far in the opposite of the expected y-direction
              if(EDIST.gt.EMAX) cycle								      ! Skip this (t) because it is too far from its expected position
            end if
          
            validind = findloc(vID(f,:),ID(prior,t),1)							      ! Get forecast error information for this ID	 
            if(validind.gt.0) then
              ERROR_DIST = haversine(Slat(f,s),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
              ERROR_DX = haversine(vSlat(f,validind),Slon(f,s),vSlat(f,validind),vSlon(f,validind))
              ERROR_DY = haversine(Slat(f,s),vSlon(f,validind),vSlat(f,validind),vSlon(f,validind))
              
	      VAL_LAT = vSlat(f,validind)
	      VAL_LON = vSlon(f,validind)
	      VAL_SO = vSo(f,validind)
              VAL_RO = vRo(f,validind)
	      VAL_ZMIN = vZmin(f,validind)
	    
              if((Slat(f,s)-vSlat(f,validind)).lt.0) ERROR_DY = -1.0*ERROR_DY     
              if(((Slon(f,s)-vSlon(f,validind)).lt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.-100.0)) ERROR_DX = -1.0*ERROR_DX
              if(((Slon(f,s)-vSlon(f,validind)).gt.0).and.((Slon(f,s)-vSlon(f,validind)).gt.100.0))  ERROR_DX = -1.0*ERROR_DX
              
              holdSoFlag = visFMaxSo(f,validind)										      ! Save whether the verification feature is at its MaxSo or MinZ time
              holdZFlag = visFMinZ(f,validind)
            else
              ERROR_DIST = -99999.9
	      ERROR_DX = -99999.9
              ERROR_DY = -99999.9
            
	      VAL_LAT = -9999.9
	      VAL_LON = -9999.9 
	      VAL_SO = -9999.9
              VAL_RO = -9999.9
	      VAL_ZMIN = -9999.9
	      
              holdSoFlag = -1										      
              holdZFlag = -1
            end if
        	
            SPTY = abs(So(f,s)-So(prior,t))/SNORM
            BPTY = abs(BGo(f,s)-BGo(prior,t))/BNORM
            RPTY = abs(Ro(f,s)-Ro(prior,t))/RNORM
            PSUM = SPTY+BPTY+RPTY
            !if((PSUM-OVERLAP).gt.PMAX) cycle
            if(PSUM.gt.PMAX) cycle
            				        						       ! Compare the combined penalty & overlap term to the previous "TARGET" combined term to beat
	    !if((PSUM-OVERLAP).lt.TARGET) then
	    if(PSUM.lt.TARGET) then
              ID(f,s)       = ID(prior,t)								      ! Assign a new history and ID to feature (n,s)
              nid(f,s)      = nid(prior,t) + 1
              TrackMatch(f,s)		= 0
	      										      ! Assign characteristic history
              hstry_So(f,s,:) 	      	= hstry_So(prior,t,:)
	      hstry_BG(f,s,:) 	      	= hstry_BG(prior,t,:)
	      hstry_Ro(f,s,:) 	      	= hstry_Ro(prior,t,:)
              hstry_So(f,s,nid(f,s))    = So(f,s)
	      hstry_BG(f,s,nid(f,s))    = BGo(f,s)
	      hstry_Ro(f,s,nid(f,s))    = Ro(f,s)
        											      ! Assign position history
              hstry_la(f,s,:) 	      	= hstry_la(prior,t,:)
	      hstry_lo(f,s,:) 	      	= hstry_lo(prior,t,:)
	      hstry_dx(f,s,:) 	      	= hstry_dx(prior,t,:)
              hstry_dy(f,s,:) 	      	= hstry_dy(prior,t,:)	
              hstry_la(f,s,nid(f,s))    = Slat(f,s)
	      hstry_lo(f,s,nid(f,s))    = Slon(f,s)
	      hstry_dx(f,s,nid(f,s))    = DX
              hstry_dy(f,s,nid(f,s))    = DY
        											      ! Compute running averages of histories
              !runav_So(f,s)  	      	= sum(hstry_So(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 
              !runav_BG(f,s)  	      	= sum(hstry_BG(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 
              !runav_Ro(f,s)  	      	= sum(hstry_Ro(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/min(nid(f,s),4) 

	     ! runav_la(f,s)  	      	= sum(hstry_la(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/float(min(nid(f,s),4))
	      !runav_lo(f,s)  	      	= sum(hstry_lo(f,s,(max(1,nid(f,s)-3)):nid(f,s)))/float(min(nid(f,s),4))
	      !runav_dx(f,s)  	      	= sum(hstry_dx(f,s,(max(2,nid(f,s)-2)):nid(f,s)))/float(min(nid(f,s),3))
	      !runav_dy(f,s)  	      	= sum(hstry_dy(f,s,(max(2,nid(f,s)-2)):nid(f,s)))/float(min(nid(f,s),3))
            
              saveDX(f,s)	      	= DX
              saveDY(f,s)	      	= DY
              saveDIST(f,s)	      	= DISTANCE
              saveDT(f,s)	      	= 12.0
	      saveVLat(f,s)		= VAL_LAT
	      saveVLon(f,s)		= VAL_LON            
	      saveVSo(f,s)		= VAL_SO
	      saveVRo(f,s)		= VAL_RO
	      saveVZmin(f,s)		= VAL_ZMIN
	      
              hstry_dt(f,s,:) 	      	= hstry_dt(prior,t,:)
              hstry_dt(f,s,nid(f,s)) 	= saveDT(f,s)
              saveDUR(f,s)	      	= sum(hstry_dt(f,s,1:nid(f,s)))
              TARGET		      	= PSUM	!-OVERLAP
              FTARGET(f,s)	      	= TARGET
              saveFTARG(f,s)  	      	= TARGET
            
              FERRY(f,s)	      	= ERROR_DY
              FERRX(f,s)	      	= ERROR_DX
              FERR(f,s)	      		= ERROR_DIST
              isFMaxSo(f,s)	      	= holdSoFlag
              isFMinZ(f,s)	      	= holdZFlag
              maxID(f,s)	      	= maxID(prior,t)
	
            end if !TARGET if()
          end do ! t do()
        end if ! Second attempt if()
      end if ! f.gt.2 if()	    
    
      if(count(ID(f,:).eq.ID(f,s)).gt.1.and.ID(f,s).gt.0) then				      ! Check if multiple (f,n,s) features matched the same (f-1,n,t) feature
	TESTID = ID(f,s)
        a = minloc(FTARGET(f,:),1,MASK = ID(f,:).eq.ID(f,s))
          
        do t=1,nsfmax
          if(ID(f,t).eq.ID(f,s).and.t.ne.a) then
            ID(f,t)		      	= -1*TESTID
	    TrackMatch(f,t)		= -1
            nid(f,t)  	      		= 0
            FTARGET(f,t)    		= 9999.9
            d2prev(f,t)     		= 0.0
            o2prev(f,t)     		= 0.0
            !hstry_So(f,t,:)	  	= 0.0
            !hstry_BG(f,t,:)	  	= 0.0
            !hstry_Ro(f,t,:)	  	= 0.0
            hstry_dx(f,t,:)	  	= 0.0
            hstry_dy(f,t,:)	  	= 0.0
            hstry_la(f,t,:)	  	= 0.0
            hstry_lo(f,t,:)	  	= 0.0
           ! runav_So(f,n,t)   		= 0.0
           ! runav_BG(f,n,t)   		= 0.0
           ! runav_Ro(f,n,t)   		= 0.0
           ! runav_dx(f,n,t)   		= 0.0
           ! runav_dy(f,n,t)   		= 0.0
           ! runav_la(f,n,t)   		= 0.0
           ! runav_lo(f,n,t)   		= 0.0
          
            saveDX(f,t)	      		= -9999.9
            saveDY(f,t)	      		= -9999.9
            saveDIST(f,t)	      	= -9999.9
	    saveVLat(f,t)		= -9999.9
	    saveVLon(f,t)		= -9999.9
            saveVSo(f,t)		= -9999.9
	    saveVRo(f,t)		= -9999.9
	    saveVZmin(f,t)		= -9999.9
	    
	    saveDT(f,t)	      		= 0.0
            hstry_dt(f,t,:)	      	= 0.0
            saveDUR(f,t)    		= 0.0
            saveFTARG(f,t)  		= 9999.9
            
            FERRY(f,t)       		= -99999.9
            FERRX(f,t)       		= -99999.9
            FERR(f,t) 	     		= -99999.9 
            
            isFMaxSo(f,t)   		= -1
            isFMinZ(f,t)    		= -1
            
            maxID(f,t)      		= -9999.9 
          end if
	end do
      end if  
    
      if(ipass.eq.npass.and.ID(f,s).le.-1) then						      ! End of the line. If the feature has not been tagged (no match whatsoever, only matched to a forbidden feature, or no match to a future analysis time), it will tested any further at this step. 	    
	ID(f,s) = -1
      end if
    
    end do    ! s      
  end do      ! pass
  
  do s=1,nsfmax
    if(itime(f,s).eq."-9999") exit
    if(ID(f,s).le.-1) then									      ! Do this again after the pass loop in case some feature (f,n,1:s-1) was set to a negative value in the comparison step at ipass=npass.	    
      ID(f,s) = -1
    end if
  end do      ! s
  
  
  
  !!! For features that have not successfully been tracked from f-1 (or f-2) to f, attempt to match them to corresponding valid time features so they can be tracked to the next fhour
  
  
  
  
  
  
  
  do ipass=1,npass															! START DO NUM PASSES    
      do s=1,nsfmax															! START DO NUM FEATURES s AT TIME n, FHOUR f
        if(itime(f,s).eq."-9999") exit     												! IF AT A "MISSING" itime, exit the s loop (no more features in the file)
       	if(ID(f,s).gt.0) cycle
        TARGET = 9999.9													      	! A PERFECT MATCH WOULD BE 1.0
        do t=1,nsmax        
	  if(vitime(f,t).eq."-9999") exit											      	! EXIT IF NO MORE VERIF FEATURES
	  if(vID(f,t).lt.0) cycle												      	! CYCLE IF VERIF ID IS LT 0 (impossible in v7, but should keep)
	  if(-1*vID(f,t).eq.ID(f,s)) cycle										      	! CYCLE IF VERIF ID IS A FORBIDDEN ID OF FEATURE
	  if(any(ID(1:f,:).eq.vID(f,t))) cycle											! CYCLE IF VERIF ID has already been tracked to the current time or was tracked and lost during previous forecast hours
	  
	  if(abs(Slat(f,s)-vSlat(f,t)).gt.30.0)then										! This is a safeguard to prevent potentially absurd matches in polar regions (but still allow for reasonable matches near 0 deg)
	    if(.not.(((Slat(f,s).ge.345.0).or.(Slat(f,s).le.15.0)).and.((vSlat(f,t).ge.345.0).or.(vSlat(f,t).le.15.0)))) cycle
	  end if
	  
	  ERROR_DIST = haversine(Slat(f,s),Slon(f,s),vSlat(f,t),vSlon(f,t))
	  ERROR_DX = haversine(vSlat(f,t),Slon(f,s),vSlat(f,t),vSlon(f,t))
	  ERROR_DY = haversine(Slat(f,s),vSlon(f,t),vSlat(f,t),vSlon(f,t))
	  
	  VAL_LAT = vSlat(f,t)
	  VAL_LON = vSlon(f,t)
	  VAL_SO = vSo(f,t)
          VAL_RO = vRo(f,t)
	  VAL_ZMIN = vZmin(f,t)
	  
	  
	  if((Slat(f,s)-vSlat(f,t)).lt.0) ERROR_DY = -1.0*ERROR_DY      
	  if(((Slon(f,s)-vSlon(f,t)).lt.0).and.((Slon(f,s)-vSlon(f,t)).gt.-100.0)) ERROR_DX = -1.0*ERROR_DX
	  if(((Slon(f,s)-vSlon(f,t)).gt.0).and.((Slon(f,s)-vSlon(f,t)).gt.100.0))  ERROR_DX = -1.0*ERROR_DX
			
	  OVERLAP = c1c2_overlap(ERROR_DIST,Ro(f,s),vRo(f,t))
	  
	 
	
	
	!!!! Pre-June 2022 correction    
	  !!!!!! COMPUTE PENALTY TERM BASED ON KS21 DIAGNOSTICS (So, BGox, BGoy, Ro)
	    
	!!  SPTY = abs(So(f,s)-vSo(f,t))/SNORM_m
	!!  BPTY = (abs(BGoy(f,s)-vBGoy(f,t))/BNORMXY_m) + (abs(BGox(f,s)-vBGox(f,t))/BNORMXY_m)
	!!  RPTY = abs(Ro(f,s)-vRo(f,t))/RNORM_m
	    
	!!  KPTY = (SPTY + BPTY + RPTY)/4.0
	    
	!!  !!!!!!!if(KPTY.gt.KMAX) cycle
	    
	!!  !!!!!! COMPUTE PENALTY TERM BASED ON POSITION ERROR (ERROR_DIST, OVERLAP)
	    
	!!  !!!!!!!!if(OVERLAP.gt.0.0)then												      	! IF THERE IS NO OVERLAP, PENALIZE BY THE DISTANCE BETWEEN THE "EDGE" OF THE FEATURES' RADII, NORMALIZED BY THE 
	!!  !!!!!!!  DPTY = 0.0													      	! ...VERIFICATION RADIUS. DPTY=1.0 IS AN EDGE-TO-EDGE DISTANCE EQUAL TO THE VERIFICATION RADIUS
	!!  !!!!!!!!else
	!!   DPTY = (ERROR_DIST/DMAX_m) !-(1.0*OVERLAP) !(2.0*(vRo(VER,t)+Ro(f,n,s))) !(ERROR_DIST-(Ro(f,n,s)+vRo(VER,t)))/(0.5*(vRo(VER,t)+Ro(f,n,s)))			      
	  !!!!!!!!!!!end if
	    
	    
	  !!!!!! COMPUTE PENALTY TERM BASED ON METEOROLOGICAL CHARACTERISTICS (850,500,200 Z,T,U,V,Q)
	    
	 !!   Z850PTY = abs(z850(f,s)-vz850(f,t))/Z850NORM
	 !!   T850PTY = abs(t850(f,s)-vt850(f,t))/T850NORM
	 !!   U850PTY = abs(u850(f,s)-vu850(f,t))/U850NORM
	 !!   V850PTY = abs(v850(f,s)-vv850(f,t))/V850NORM
	 !!   Q850PTY = abs(mr850(f,s)-vmr850(f,t))/Q850NORM
	    
	 !!   Z500PTY = abs(z500(f,s)-vz500(f,t))/Z500NORM
	 !!   T500PTY = abs(t500(f,s)-vt500(f,t))/T500NORM
	 !!   U500PTY = abs(u500(f,s)-vu500(f,t))/U500NORM
	 !!   V500PTY = abs(v500(f,s)-vv500(f,t))/V500NORM
	 !!   Q500PTY = abs(mr500(f,s)-vmr500(f,t))/Q500NORM
	    
	 !!   Z200PTY = abs(z200(f,s)-vz200(f,t))/Z200NORM
	 !!   T200PTY = abs(t200(f,s)-vt200(f,t))/T200NORM
	 !!   U200PTY = abs(u200(f,s)-vu200(f,t))/U200NORM
	 !!   V200PTY = abs(v200(f,s)-vv200(f,t))/V200NORM
	 !!   Q200PTY = abs(mr200(f,s)-vmr200(f,t))/Q200NORM 
	    
	 !!   ZPTY = Z850PTY + Z500PTY + Z200PTY
	 !!   TPTY = T850PTY + T500PTY + T200PTY
	 !!   UPTY = U850PTY + U500PTY + U200PTY
	 !!   VPTY = V850PTY + V500PTY + V200PTY
	 !!   QPTY = Q850PTY + Q500PTY + Q200PTY
	    
	 !!   METPTY = (ZPTY + TPTY + UPTY + VPTY + QPTY)/15.0
	    
	    !!!!!!!!if(METPTY.gt.METMAX) cycle
	    
	 !!   TOTPTY = ((1.0*METPTY) + (0.8*KPTY) + (1.2*DPTY))/3.0
	 !!   if(TOTPTY.gt.PMAX_m) cycle										      
	 
	 
	 
	 !!!!!June 2022 correction - 500 hPa only over a fixed radius, adjusted distance and diagnostic penalty terms
	   !!!!!! COMPUTE PENALTY TERM BASED ON KS21 DIAGNOSTICS (So, BGo, Ro)
	    
	    SPTY = abs(So(f,s)-vSo(f,t))/SNORM_m
	    BPTY = abs(BGo(f,s)-vBGo(f,t))/BNORMXY_m
	    RPTY = abs(Ro(f,s)-vRo(f,t))/RNORM_m	    
	    KPTY = (SPTY + BPTY + RPTY)/3.0
	    
	    !!!!!! COMPUTE NORMALIZED DISTANCE PENALTY
	    DPTY = (ERROR_DIST/DMAX_m)
	     
	    !!!!!! COMPUTE NORMALIZED MET PENALTY
            Z500PTY = abs(FIXEDz500(f,s)-vFIXEDz500(f,t))/Z500NORM
	    T500PTY = abs(FIXEDt500(f,s)-vFIXEDt500(f,t))/T500NORM
	    U500PTY = abs(FIXEDu500(f,s)-vFIXEDu500(f,t))/U500NORM
	    V500PTY = abs(FIXEDv500(f,s)-vFIXEDv500(f,t))/V500NORM
	    Q500PTY = abs(FIXEDmr500(f,s)-vFIXEDmr500(f,t))/Q500NORM
	    METPTY = (Z500PTY+T500PTY+U500PTY+V500PTY+Q500PTY)/5.0	 
	    TOTPTY = METPTY + KPTY + DPTY
	    if(TOTPTY.gt.PMAX_m) cycle
	    
	    
	    ! Find a valid point for the verification ID at the previous timestep. If this exists, use this point and the prospective matched point to compute DX,DY,DT
	    validind = findloc(vID(f-1,:),vID(f,t),1)							     	 
            if(validind.gt.0) then
              DISTANCE = haversine(Slat(f,s),Slon(f,s),vSlat(f-1,validind),vSlon(f-1,validind))
              DX = haversine(vSlat(f-1,validind),Slon(f,s),vSlat(f-1,validind),vSlon(f-1,validind))
              DY = haversine(Slat(f,s),vSlon(f-1,validind),vSlat(f-1,validind),vSlon(f-1,validind))
            
    
            else
              DISTANCE = -9999.9
	      DX = -9999.9
              DY = -9999.9	 
              
            end if
	    
	    
	    
	    
	    if(TOTPTY.lt.TARGET) then
	      ID(f,s) 	      	= vID(f,t)
	      nid(f,s)		= 1
	      TrackMatch(f,s)		= 1	  
	      TARGET		      	= TOTPTY !OfLAP-PSUM
	      FTARGET(f,s)	      	= TARGET
	      saveFTARG(f,s)        	= TARGET
	    
	      FERRY(f,s)	      	= ERROR_DY
	      FERRX(f,s)	      	= ERROR_DX
	      FERR(f,s)	      		= ERROR_DIST	    
	      maxID(f,s)	      	= vmaxID(f,t)
	      isFMaxSo(f,s)         	= visFMaxSo(f,t)
	      isFMinZ(f,s)	      	= visFMinZ(f,t)
	      saveDUR(f,s)		= vsaveDUR(f,t)
	      hstry_dt(f,s,1)		= saveDUR(f,s)		! Save the accumulated dt to this point to it can be applied to features tracked from the match
	      saveDT(f,s)		= vsaveDT(f,t)
	      saveDX(f,s)		= DX 	!vsaveDX(f,t)
	      saveDY(f,s)		= DY	!vsaveDY(f,t)
	      saveDIST(f,s)		= DISTANCE	!vsaveDIST(f,t)
	      saveVLat(f,s)		= VAL_LAT
	      saveVLon(f,s)		= VAL_LON
	      saveVSo(f,s)		= VAL_SO
	      saveVRo(f,s)		= VAL_RO
	      saveVZmin(f,s)		= VAL_ZMIN
	      
	    end if ! TARGET if()
	  end do ! t do()
               
	
	
        if(count(ID(f,:).eq.ID(f,s)).gt.1.and.ID(f,s).gt.0) then									! Check if multiple (f,s) features matched the same valid(VER,t) feature
          TESTID = ID(f,s)
	  a = minloc(FTARGET(f,:),1,MASK = ID(f,:).eq.ID(f,s))
	    
	  do t=1,nsfmax
	    if(ID(f,t).eq.ID(f,s).and.t.ne.a) then
	      ID(f,t) 		= -1*TESTID
	      nID(f,t)		= 0
	      FTARGET(f,t) 		= 9999.9
	      TrackMatch(f,t)		= -1
	      saveDT(f,t)		= 0.0
	      saveDUR(f,t)		= 0.0
	      saveFTARG(f,t)		= 9999.9
	      saveDX(f,t)		= -9999.9
	      saveDY(f,t)		= -9999.9
	      saveDIST(f,t)		= -9999.9
	      hstry_dt(f,t,:)		= 0.0
	      
	      FERRY(f,t)       		= -99999.9
	      FERRX(f,t)       		= -99999.9
	      FERR(f,t)	      		= -99999.9 
	      saveVLat(f,t)		= -9999.9
	      saveVLon(f,t)		= -9999.9
	      saveVSo(f,t)		= -9999.9
	      saveVRo(f,t)		= -9999.9
	      saveVZmin(f,t)		= -9999.9
	      isFMaxSo(f,t)		= -1
	      isFMinZ(f,t)		= -1	      
	      maxID(f,t)		= -9999.9  
	    end if
          end do
        end if  
      
        if(ipass.eq.npass.and.ID(f,s).le.-1) then											! End of the line. If the feature has not been tagged (no match whatsoever, only matched to a forbidden feature, or no match to a future analysis time), it will tested any further at this step.	      
          ID(f,s) = -1
        end if
      
      end do 	! s      
    end do 	! pass
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end do        ! f
	


call cpu_time(FN)
print '("Time to ID forecast features = ",f9.3," seconds.")',FN-ST


    
    
!!!!Write to text file!!!!
call cpu_time(ST)
do f=2,nf
  !do t=1,nt
    write(*,*)'outfile=',trim(outfile(f))
    open(10, file=trim(outfile(f)), status="unknown")

    fstr = "(A7,1x,A7,1x,A7,1x,A11,1x,A7,1x,A7,1x,A8,1x,A7,1x,A4,1x,A15,1x,A16,1x,A16,1x,A7,1x,A7,1x,A7,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A17,1x,A17,1x,A17,1x,A17,1x,A17,1x,A6,1x,A6,1x,A8,1x,A6,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)"
    write(10,fstr) "ITIME","FHOUR","ID","So(m/100km)","LAT(N)","LON(E)","SoFlag","Ro(km)","SR","BGo(m/100km)","BGo-lat(m/100km)","BGo-lon(m/100km)","ZMIN(m)","ZLAT(N)","ZLON(E)","ZFlag","Z850(m)","Z500(m)","Z200(m)","T850(K)","T500(K)","T200(K)","U850(m/s)","U500(m/s)","U200(m/s)","V850(m/s)","V500(m/s)","V200(m/s)","MR850(g/kg)","MR500(g/kg)","MR200(g/kg)","600kmZ500(m)","600kmT500(K)","600kmU500(m/s)","600kmV500(m/s)","600kmMR500(g/kg)","DY(km)","DX(km)","DIST(km)","DT(h)","DUR(h)","MAXDUR(h)","PTY-OVR", "FERRY(km)","FERRX(km)","FERR(km)","T(0)/M(1)/N(-1)","VLat(N)","VLon(E)","VSo","VRo","VZmin"

    fstr = "(A11,1x,A4,1x,I9,1x,F6.2,1x,F6.2,1x,F6.2,1x,I3,1x,F10.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,I3,1x,F8.2,1x,F8.2,1x,F8.2,1x,F5.1,1x,F5.1,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.4,1x,F7.4,1x,F7.4,1x,F8.2,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.4,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F10.2,1x,F10.2,1x,F10.2,1x,I4,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2)"
    do s=1,nsfmax
      if(itime(f,s).eq."-9999") exit  
      write(10,fstr) itime(f,s), fhour(f,s), ID(f,s), So(f,s), Slat(f,s), Slon(f,s), isFMaxSo(f,s), Ro(f,s), SR(f,s), BGo(f,s), BGoy(f,s), BGox(f,s), Zmin(f,s), Zlat(f,s), Zlon(f,s), isFMinZ(f,s), z850(f,s), z500(f,s), z200(f,s), t850(f,s), t500(f,s), t200(f,s), u850(f,s), u500(f,s), u200(f,s), v850(f,s), v500(f,s), v200(f,s), mr850(f,s), mr500(f,s), mr200(f,s), FIXEDz500(f,s), FIXEDt500(f,s), FIXEDu500(f,s), FIXEDv500(f,s), FIXEDmr500(f,s), saveDY(f,s), saveDX(f,s), saveDIST(f,s), saveDT(f,s),  saveDUR(f,s),maxID(f,s),saveFTARG(f,s),FERRY(f,s),FERRX(f,s),FERR(f,s),TrackMatch(f,s),saveVLat(f,s),saveVLon(f,s),saveVSo(f,s),saveVRo(f,s),saveVZmin(f,s)
    end do
    close(10)
  !end do
end do
call cpu_time(FN)
print '("Time to write forecast track files = ",f9.3," seconds.")',FN-ST


















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


