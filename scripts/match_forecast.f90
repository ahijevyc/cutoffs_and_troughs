!match_forecast.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Kevin Lupo
!klupo@ucar.edu
!match cutoff lows and preexisting troughs identified using the KS21 algorithm at different fhours to tracked analysis features
!
!15 Feb 2022
!	Starting from 13-14 Feb code from track_forecast.f90
!
!16 Feb 2022
!	Adjusted the background penalty term to account for x- and y- components rather than the the total BGo
!
!19 May 2022
!	Adjusted the distance normalization to 750 km (was 1500km) to match the track-match code. 
!	Adjust weighting of kasuga, met, and distance penalty terms (1.0 met, 0.8 kasuga, 1.2 distance) - penalize the distance most harshly, kasuga least, and met middle
!	Some definitions - KPTY is the average of the normalized So, Ro, BGx, and BGy penalties (4 terms)
!			 - METPTY is the average of the normalized penalties for all met terms at each level (15 terms)
!			 - DPTY is the normalized distance penalty (1 term)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
program match_forecast
use netcdf

implicit none

integer				:: nx, ny, nf, ns, nt, nsmax, nsfmax, s, f, t, tt, ss, n
integer				:: i, ii, j, k, x, y, r, a, b, np1, np2, nz1, np
integer				:: i0, i1, j0, j1, xx0, xx1
integer				:: istat,istat2,istat3
real				:: d, hgtmin, dmin,dthresh
real				:: haversine, c1c2_overlap
real				:: ST, FN
integer, allocatable		:: ID(:,:,:),TrackMatch(:,:,:)		!, nid(:,:,:)
!real, allocatable		:: d2prev(:,:,:),o2prev(:,:,:)
real, allocatable		:: So(:,:,:), Slat(:,:,:), Slon(:,:,:), Ro(:,:,:), SR(:,:,:), BGo(:,:,:), BGoy(:,:,:), BGox(:,:,:), Zmin(:,:,:), Zlat(:,:,:), Zlon(:,:,:)
real, allocatable		:: z850(:,:,:), z500(:,:,:), z200(:,:,:), t850(:,:,:), t500(:,:,:), t200(:,:,:), u850(:,:,:), u500(:,:,:), u200(:,:,:), v850(:,:,:), v500(:,:,:), v200(:,:,:), mr850(:,:,:), mr500(:,:,:), mr200(:,:,:)
integer, allocatable		:: isFMaxSo(:,:,:),isFMinZ(:,:,:)

integer, allocatable		:: vID(:,:)
real, allocatable		:: vSo(:,:), vSlat(:,:), vSlon(:,:), vRo(:,:), vSR(:,:), vBGo(:,:), vBGoy(:,:), vBGox(:,:), vZmin(:,:), vZlat(:,:), vZlon(:,:)
real, allocatable		:: vz850(:,:), vz500(:,:), vz200(:,:), vt850(:,:), vt500(:,:), vt200(:,:), vu850(:,:), vu500(:,:), vu200(:,:), vv850(:,:), vv500(:,:), vv200(:,:), vmr850(:,:), vmr500(:,:), vmr200(:,:)
integer, allocatable		:: visFMaxSo(:,:),visFMinZ(:,:)
integer				:: holdSoFlag, holdZFlag

!real, allocatable		:: runav_So(:,:,:),runav_BG(:,:,:),runav_Ro(:,:,:),runav_la(:,:,:),runav_lo(:,:,:),runav_dx(:,:,:),runav_dy(:,:,:)
!real, allocatable		:: hstry_So(:,:,:,:),hstry_BG(:,:,:,:),hstry_Ro(:,:,:,:), hstry_la(:,:,:,:),hstry_lo(:,:,:,:),hstry_dx(:,:,:,:),hstry_dy(:,:,:,:)

!real, allocatable		:: hstry_dt(:,:,:,:)
real, allocatable		:: saveDUR(:,:,:),FERRX(:,:,:),FERRY(:,:,:),FERR(:,:,:),saveFTARG(:,:,:)

integer, parameter		:: npass = 4
integer				:: ipass
real				:: CSOMIN,TSOMIN,OMIN,PMAX,DMAX,EMAX,OXMAX,OYMAX,SNORM,BNORM,RNORM,BNORMXY
real				:: OVERLAP,DISTANCE
real				:: PDX,PDY,DX,DY,EDIST,OYRAT,OXRAT,pSlat,pSlon
real				:: SPTY,BPTY,RPTY,DPTY,PSUM
real				:: TARGET
real, allocatable		:: FTARGET(:,:,:)
integer				:: prior,validind,TESTID
real, allocatable		:: saveDX(:,:,:),saveDY(:,:,:),saveDIST(:,:,:),saveDT(:,:,:),saveVLat(:,:,:),saveVLon(:,:,:),saveVSo(:,:,:),saveVRo(:,:,:),saveVZmin(:,:,:)
real, allocatable		:: vsaveDX(:,:),vsaveDY(:,:),vsaveDIST(:,:),vsaveDT(:,:)
real, allocatable		:: vsaveDUR(:,:),vFERRX(:,:),vFERRY(:,:),vFERR(:,:),vsaveFTARG(:,:)
real				:: ERROR_DIST,ERROR_DY,ERROR_DX,VAL_LAT,VAL_LON,VAL_SO,VAL_RO,VAL_ZMIN
real, allocatable		:: vmaxID(:,:),maxID(:,:,:)

real				:: Z850NORM, Z500NORM, Z200NORM, Z850PTY, Z500PTY, Z200PTY
real				:: T850NORM, T500NORM, T200NORM, T850PTY, T500PTY, T200PTY
real				:: U850NORM, U500NORM, U200NORM, U850PTY, U500PTY, U200PTY
real				:: V850NORM, V500NORM, V200NORM, V850PTY, V500PTY, V200PTY
real				:: Q850NORM, Q500NORM, Q200NORM, Q850PTY, Q500PTY, Q200PTY
real				:: ZPTY, TPTY, UPTY, VPTY, QPTY, METPTY, METMAX
real				:: KPTY, KMAX, TOTPTY


integer				:: fnum, VER
character(len=10)		:: fnum_char
character(len=80)		:: name
character(len=270)		:: ffilelist
character(len=270)		:: line, line2, line3
character(len=270)		:: vfile, listfile, infile, ffile
character(len=20)		:: varstr
!character(len=10)		:: YYYY
!character(len=10)		:: MM
!character(len=5)		:: DD
!character(len=5)		:: hh
!character(len=10)		:: hemis
!character(len=5)		:: hext
character(len=20),allocatable   :: fhour(:,:,:),vfhour(:,:)
character(len=300),allocatable	:: outfile(:,:)
character(len=400)		:: fstr
character(len=20),allocatable	:: itime(:,:,:),vitime(:,:)
character(len=300)		:: outdir

!!!!!!Get some command line arguments
call get_command_argument(1,vfile)
call get_command_argument(2,ffile)
call get_command_argument(3,varstr)
!call get_command_argument(4,hemis)
!call get_command_argument(5,hext)
call get_command_argument(4,fnum_char)
read(fnum_char , * ) fnum		! Since this is now a batch submission script, we need the correct number of steps to add to n to get to the valid time. In this configuration, nf will always be 1.


!!!!!Name the output text file
outdir = "/glade/work/klupo/postdoc/kasugaEA21/version8/" // trim(varstr) // "/" 
print *, outdir

!!!!!Get dimensions for forecast lists!!!!!
call cpu_time(ST)
nf    = 0
nsfmax = 0
open(unit = 10, file = trim(ffile), status='old', access='sequential', form='formatted', action='read')
do
  read(10,'(a256)',iostat=istat) listfile      
  if(istat /= 0)exit

  nt=0
  open(unit = 11, file = trim(listfile), status='old', access='sequential',form='formatted', action='read')
  do
    read(11,'(a256)',iostat=istat2) infile
    if(istat2 /=0)exit
    
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
    nt=nt+1
  end do
  
  close(11)  
  nf=nf+1
end do
close(10)


allocate(itime(nf,nt,nsfmax),fhour(nf,nt,nsfmax),ID(nf,nt,nsfmax),TrackMatch(nf,nt,nsfmax),So(nf,nt,nsfmax),Slat(nf,nt,nsfmax),Slon(nf,nt,nsfmax),Ro(nf,nt,nsfmax),SR(nf,nt,nsfmax),BGo(nf,nt,nsfmax),BGoy(nf,nt,nsfmax),BGox(nf,nt,nsfmax),Zmin(nf,nt,nsfmax),Zlat(nf,nt,nsfmax),Zlon(nf,nt,nsfmax))
allocate(outfile(nf,nt))
allocate(z850(nf,nt,nsfmax), z500(nf,nt,nsfmax), z200(nf,nt,nsfmax), t850(nf,nt,nsfmax), t500(nf,nt,nsfmax), t200(nf,nt,nsfmax), u850(nf,nt,nsfmax), u500(nf,nt,nsfmax), u200(nf,nt,nsfmax), v850(nf,nt,nsfmax), v500(nf,nt,nsfmax), v200(nf,nt,nsfmax), mr850(nf,nt,nsfmax), mr500(nf,nt,nsfmax), mr200(nf,nt,nsfmax))
!allocate(nid(nf,nt,nsfmax))
!allocate(runav_So(nf,nt,nsfmax),runav_BG(nf,nt,nsfmax),runav_Ro(nf,nt,nsfmax))
allocate(FTARGET(nf,nt,nsfmax))
allocate(saveDX(nf,nt,nsfmax),saveDY(nf,nt,nsfmax),saveDIST(nf,nt,nsfmax),saveDT(nf,nt,nsfmax),saveVLat(nf,nt,nsfmax),saveVLon(nf,nt,nsfmax),saveVSo(nf,nt,nsfmax),saveVRo(nf,nt,nsfmax),saveVZmin(nf,nt,nsfmax))
!allocate(hstry_So(nf,nt,nsfmax,300),hstry_BG(nf,nt,nsfmax,300),hstry_Ro(nf,nt,nsfmax,300))
!allocate(runav_la(nf,nt,nsfmax),runav_lo(nf,nt,nsfmax),runav_dx(nf,nt,nsfmax),runav_dy(nf,nt,nsfmax))
!allocate(hstry_la(nf,nt,nsfmax,300),hstry_lo(nf,nt,nsfmax,300),hstry_dx(nf,nt,nsfmax,300),hstry_dy(nf,nt,nsfmax,300))
allocate(isFMaxSo(nf,nt,nsfmax),isFMinZ(nf,nt,nsfmax))

allocate(saveDUR(nf,nt,nsfmax))
allocate(saveFTARG(nf,nt,nsfmax))
allocate(FERRX(nf,nt,nsfmax),FERRY(nf,nt,nsfmax),FERR(nf,nt,nsfmax))
allocate(maxID(nf,nt,nsfmax))

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

allocate(vitime(nt,nsmax),vfhour(nt,nsmax),vID(nt,nsmax),vSo(nt,nsmax),vSlat(nt,nsmax),vSlon(nt,nsmax),vRo(nt,nsmax),vSR(nt,nsmax),vBGo(nt,nsmax),vBGoy(nt,nsmax),vBGox(nt,nsmax),vZmin(nt,nsmax),vZlat(nt,nsmax),vZlon(nt,nsmax))
allocate(vsaveDX(nt,nsmax),vsaveDY(nt,nsmax),vsaveDIST(nt,nsmax),vsaveDT(nt,nsmax))
allocate(vz850(nt,nsmax), vz500(nt,nsmax), vz200(nt,nsmax), vt850(nt,nsmax), vt500(nt,nsmax), vt200(nt,nsmax), vu850(nt,nsmax), vu500(nt,nsmax), vu200(nt,nsmax), vv850(nt,nsmax), vv500(nt,nsmax), vv200(nt,nsmax), vmr850(nt,nsmax), vmr500(nt,nsmax), vmr200(nt,nsmax))
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
FERRY = -9999.9
FERRX = -9999.9
FERR = -9999.9
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
vFERRY = -9999.9
vFERRX = -9999.9
vFERR = -9999.9




!!!!!Read Verif vars from .track files!!!!!
call cpu_time(ST)
open(unit = 11, file = trim(vfile), status='old', access='sequential', form='formatted', action='read')
do t=1,nt
  read(11,'(a256)',iostat=istat) infile     
  if(istat /= 0)exit
  
  open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
  do s=0,nsmax
    if(s.eq.0)then
      read(12,*,iostat=istat2) line
    else
      read(12,*,iostat=istat2) vitime(t,s), vfhour(t,s), vID(t,s), vSo(t,s), vSlat(t,s), vSlon(t,s), visFMaxSo(t,s), vRo(t,s), vSR(t,s), vBGo(t,s), vBGoy(t,s), vBGox(t,s), vZmin(t,s), vZlat(t,s), vZlon(t,s), visFMinZ(t,s), vz850(t,s), vz500(t,s), vz200(t,s), vt850(t,s), vt500(t,s), vt200(t,s), vu850(t,s), vu500(t,s), vu200(t,s), vv850(t,s), vv500(t,s), vv200(t,s), vmr850(t,s), vmr500(t,s), vmr200(t,s), vsaveDY(t,s), vsaveDX(t,s), vsaveDIST(t,s), vsaveDT(t,s), vsaveDUR(t,s),vmaxID(t,s), vsaveFTARG(t,s),vFERRY(t,s),vFERRX(t,s),vFERR(t,s)
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
  read(10,'(a256)',iostat=istat) listfile      
  if(istat /= 0)exit

  open(unit = 11, file = trim(listfile), status='old', access='sequential',form='formatted', action='read')
  do t=1,nt
    read(11,'(a256)',iostat=istat2) infile
    if(istat2 /=0)exit
    
    open(unit = 12, file = trim(infile), status='old', access='sequential',form='formatted', action='read')
    do s=0,nsfmax
      if(s.eq.0)then
        read(12,*,iostat=istat3) line
      else
        read(12,*,iostat=istat3) itime(f,t,s), fhour(f,t,s), So(f,t,s), Slat(f,t,s), Slon(f,t,s), Ro(f,t,s), SR(f,t,s), BGo(f,t,s), BGoy(f,t,s), BGox(f,t,s), Zmin(f,t,s), Zlat(f,t,s), Zlon(f,t,s), z850(f,t,s), z500(f,t,s), z200(f,t,s), t850(f,t,s), t500(f,t,s), t200(f,t,s), u850(f,t,s), u500(f,t,s), u200(f,t,s), v850(f,t,s), v500(f,t,s), v200(f,t,s), mr850(f,t,s), mr500(f,t,s), mr200(f,t,s)
      end if    
    
      if(istat3 /=0)exit
    end do
    
    close(12)
    outfile(f,t) = trim(outdir) // trim(infile(57:80)) // ".match" 
  end do
  
  close(11)  
end do

close(10)
call cpu_time(FN)
print '("Time gather forecast data = ",f8.3," seconds.")',FN-ST


!!!!!BELOW TEXT BLOCK IS MOSTLY UNUSED FOR SIMPLE MATCHING!!!!!!

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
PMAX   = 1.0		! Maximum allowable penalty points before a feature is disqualified
KMAX   = 1.0
DMAX   = 750.0		! Maximum distance between matches allowed (was 1700.0, was 1500 after that)
EMAX   = 1000.0		! Maximum "error" distance between predicted and candidate points
OXMAX  = 300.0		! Maximum distance a matching point can be in the opposite x-direction of the predicted point
OYMAX  = 600.0		! "..."y-direction"" More lenient to allow for features moving around larger scale troughs to move to the north
SNORM  = 6.0 !10.0		! Normalization for the So penalty
BNORM  = 15.0		! Normalization for the BGo penalty
BNORMXY= 6.0		! Normalization for the Zonal and Meridional BGo penalties
RNORM  = 750.0 !900.0		! Normalization for the Ro penalty


Z850NORM = 75.0		!100.0
Z500NORM = 150.0 	!200.0
Z200NORM = 225.0 	!300.0

T850NORM = 7.5		!10.0
T500NORM = 7.5		!10.0
T200NORM = 7.5		!10.0

U850NORM = 10.0
U500NORM = 10.0
U200NORM = 10.0

V850NORM = 5.0		! All were 10.0
V500NORM = 5.0
V200NORM = 5.0

Q850NORM = 1.5
Q500NORM = 0.75
Q200NORM = 0.015

METMAX = 1.5

do n=1,nt																! START DO NUM TIMES
  print *, itime(1,n,1)
  VER = n+fnum																! GET THE VALID TIME INDEX
  do f=1,nf																! START DO NUM FHOURS
    do ipass=1,npass															! START DO NUM PASSES    
      do s=1,nsfmax															! START DO NUM FEATURES s AT TIME n, FHOUR f
        if(itime(f,n,s).eq."-9999") exit     												! IF AT A "MISSING" itime, exit the s loop (no more features in the file)
        
	if((VER).le.nt) then
	
          TARGET = 9999.9													      	! A PERFECT MATCH WOULD BE 1.0
          do t=1,nsmax        
	    if(vitime(VER,t).eq."-9999") exit											      	! EXIT IF NO MORE VERIF FEATURES
	    if(vID(VER,t).lt.0) cycle												      	! CYCLE IF VERIF ID IS LT 0 (impossible in v7, but should keep)
	    if(-1*vID(VER,t).eq.ID(f,n,s)) cycle										      	! CYCLE IF VERIF ID IS A FORBIDDEN ID OF FEATURE
	
	    ERROR_DIST = haversine(Slat(f,n,s),Slon(f,n,s),vSlat(VER,t),vSlon(VER,t))
	    ERROR_DX = haversine(vSlat(VER,t),Slon(f,n,s),vSlat(VER,t),vSlon(VER,t))
	    ERROR_DY = haversine(Slat(f,n,s),vSlon(VER,t),vSlat(VER,t),vSlon(VER,t))
	    	    
	    VAL_LAT = vSlat(VER,t)
	    VAL_LON = vSlon(VER,t)
	    VAL_SO = vSo(VER,t)
            VAL_RO = vRo(VER,t)
	    VAL_ZMIN = vZmin(VER,t)
	    
	    if((Slat(f,n,s)-vSlat(VER,t)).lt.0) ERROR_DY = -1.0*ERROR_DY      
	    if(((Slon(f,n,s)-vSlon(VER,t)).lt.0).and.((Slon(f,n,s)-vSlon(VER,t)).gt.-100.0)) ERROR_DX = -1.0*ERROR_DX
	    if(((Slon(f,n,s)-vSlon(VER,t)).gt.0).and.((Slon(f,n,s)-vSlon(VER,t)).gt.100.0))  ERROR_DX = -1.0*ERROR_DX
			
	    OVERLAP = c1c2_overlap(ERROR_DIST,Ro(f,n,s),vRo(VER,t))
	    
	    
	    !!!!!! COMPUTE PENALTY TERM BASED ON KS21 DIAGNOSTICS (So, BGox, BGoy, Ro)
	    
	    SPTY = abs(So(f,n,s)-vSo(VER,t))/SNORM
	    BPTY = (abs(BGoy(f,n,s)-vBGoy(VER,t))/BNORMXY) + (abs(BGox(f,n,s)-vBGox(VER,t))/BNORMXY)
	    RPTY = abs(Ro(f,n,s)-vRo(VER,t))/RNORM
	    
	    KPTY = (SPTY + BPTY + RPTY)/4.0
	    
	    !if(KPTY.gt.KMAX) cycle
	    
	    !!!!!! COMPUTE PENALTY TERM BASED ON POSITION ERROR (ERROR_DIST, OVERLAP)
	    
	    !if(OVERLAP.gt.0.0)then												      	! IF THERE IS NO OVERLAP, PENALIZE BY THE DISTANCE BETWEEN THE "EDGE" OF THE FEATURES' RADII, NORMALIZED BY THE 
	    !  DPTY = 0.0													      	! ...VERIFICATION RADIUS. DPTY=1.0 IS AN EDGE-TO-EDGE DISTANCE EQUAL TO THE VERIFICATION RADIUS
	    !else
	      DPTY = (ERROR_DIST/DMAX) !-(1.0*OVERLAP) !(2.0*(vRo(VER,t)+Ro(f,n,s))) !(ERROR_DIST-(Ro(f,n,s)+vRo(VER,t)))/(0.5*(vRo(VER,t)+Ro(f,n,s)))			      
	    !end if
	    
	    
	    !!!!!! COMPUTE PENALTY TERM BASED ON METEOROLOGICAL CHARACTERISTICS (850,500,200 Z,T,U,V,Q)
	    
	    Z850PTY = abs(z850(f,n,s)-vz850(VER,t))/Z850NORM
	    T850PTY = abs(t850(f,n,s)-vt850(VER,t))/T850NORM
	    U850PTY = abs(u850(f,n,s)-vu850(VER,t))/U850NORM
	    V850PTY = abs(v850(f,n,s)-vv850(VER,t))/V850NORM
	    Q850PTY = abs(mr850(f,n,s)-vmr850(VER,t))/Q850NORM
	    
	    Z500PTY = abs(z500(f,n,s)-vz500(VER,t))/Z500NORM
	    T500PTY = abs(t500(f,n,s)-vt500(VER,t))/T500NORM
	    U500PTY = abs(u500(f,n,s)-vu500(VER,t))/U500NORM
	    V500PTY = abs(v500(f,n,s)-vv500(VER,t))/V500NORM
	    Q500PTY = abs(mr500(f,n,s)-vmr500(VER,t))/Q500NORM
	    
	    Z200PTY = abs(z200(f,n,s)-vz200(VER,t))/Z200NORM
	    T200PTY = abs(t200(f,n,s)-vt200(VER,t))/T200NORM
	    U200PTY = abs(u200(f,n,s)-vu200(VER,t))/U200NORM
	    V200PTY = abs(v200(f,n,s)-vv200(VER,t))/V200NORM
	    Q200PTY = abs(mr200(f,n,s)-vmr200(VER,t))/Q200NORM 
	    
	    ZPTY = Z850PTY + Z500PTY + Z200PTY
	    TPTY = T850PTY + T500PTY + T200PTY
	    UPTY = U850PTY + U500PTY + U200PTY
	    VPTY = V850PTY + V500PTY + V200PTY
	    QPTY = Q850PTY + Q500PTY + Q200PTY
	    
	    METPTY = (ZPTY + TPTY + UPTY + VPTY + QPTY)/15.0
	    
	    !if(METPTY.gt.METMAX) cycle
	    
	    TOTPTY = ((1.0*METPTY) + (0.8*KPTY) + (1.2*DPTY))/3.0
	    if(TOTPTY.gt.PMAX) cycle										      
	
	    if(TOTPTY.lt.TARGET) then
	      ID(f,n,s) 	      	= vID(VER,t)
	      TrackMatch(f,n,s)		= 1	  
	      TARGET		      	= TOTPTY !OVERLAP-PSUM
	      FTARGET(f,n,s)	      	= TARGET
	      saveFTARG(f,n,s)        	= TARGET
	    
	      FERRY(f,n,s)	      	= ERROR_DY
	      FERRX(f,n,s)	      	= ERROR_DX
	      FERR(f,n,s)	      	= ERROR_DIST	    
	      maxID(f,n,s)	      	= vmaxID(VER,t)
	      isFMaxSo(f,n,s)         	= visFMaxSo(VER,t)
	      isFMinZ(f,n,s)	      	= visFMinZ(VER,t)
	      saveDUR(f,n,s)		= vsaveDUR(VER,t)
	      saveDT(f,n,s)		= vsaveDT(VER,t)
	      saveDX(f,n,s)		= vsaveDX(VER,t)
	      saveDY(f,n,s)		= vsaveDY(VER,t)
	      saveDIST(f,n,s)		= vsaveDIST(VER,t)
	      saveVLat(f,n,s)		= VAL_LAT
	      saveVLon(f,n,s)		= VAL_LON
	      saveVSo(f,n,s)		= VAL_SO
	      saveVRo(f,n,s)		= VAL_RO
	      saveVZmin(f,n,s)		= VAL_ZMIN
	    end if ! TARGET if()
	  end do ! t do()
	end if
               
	
	
        if(count(ID(f,n,:).eq.ID(f,n,s)).gt.1.and.ID(f,n,s).gt.0) then									! Check if multiple (f,n,s) features matched the same valid(VER,t) feature
          TESTID = ID(f,n,s)
	  a = minloc(FTARGET(f,n,:),1,MASK = ID(f,n,:).eq.ID(f,n,s))
	    
	  do t=1,nsfmax
	    if(ID(f,n,t).eq.ID(f,n,s).and.t.ne.a) then
	      ID(f,n,t) 		= -1*TESTID
	      FTARGET(f,n,t) 		= 9999.9
	      TrackMatch(f,n,t)		= -1
	      saveDT(f,n,t)		= 0.0
	      saveDUR(f,n,t)		= 0.0
	      saveFTARG(f,n,t)		= 9999.9
	      saveDX(f,n,t)		= -9999.9
	      saveDY(f,n,t)		= -9999.9
	      saveDIST(f,n,t)		= -9999.9
	      saveVLat(f,n,t)	       	= -9999.9
	      saveVLon(f,n,t)	       	= -9999.9
	      saveVSo(f,n,t)	       	= -9999.9
	      saveVRo(f,n,t)	       	= -9999.9
	      saveVZmin(f,n,t)         	= -9999.9
	      
	      FERRY(f,n,t)       	= -9999.9
	      FERRX(f,n,t)       	= -9999.9
	      FERR(f,n,t)	      	= -9999.9 
	      
	      isFMaxSo(f,n,t)		= -1
	      isFMinZ(f,n,t)		= -1	      
	      maxID(f,n,t)		= -9999.9  
	    end if
          end do
        end if  
      
        if(ipass.eq.npass.and.ID(f,n,s).le.-1) then											! End of the line. If the feature has not been tagged (no match whatsoever, only matched to a forbidden feature, or no match to a future analysis time), it will tested any further at this step.	      
          ID(f,n,s) = -1
        end if
      
      end do 	! s      
    end do 	! pass
    
    do s=1,nsfmax
      if(itime(f,n,s).eq."-9999") exit
      if(ID(f,n,s).le.-1) then														! Do this again after the pass loop in case some feature (f,n,1:s-1) was set to a negative value in the comparison step at ipass=npass.	      
        ID(f,n,s) = -1
      end if
    end do	! s
    
  end do 	! f
end do 		! n
	


call cpu_time(FN)
print '("Time to match forecast features = ",f9.3," seconds.")',FN-ST


    
    
!!!!Write to text file!!!!
call cpu_time(ST)
do f=1,nf
  do t=1,nt
    print *,outfile(f,t)
    !print *,ID(f,t,:)
    !if(t.lt.100)cycle
    !stop
    open(10, file=trim(outfile(f,t)), status="unknown")
    !print *,"hello"
    fstr = "(A7,1x,A7,1x,A7,1x,A11,1x,A7,1x,A7,1x,A8,1x,A7,1x,A4,1x,A15,1x,A16,1x,A16,1x,A7,1x,A7,1x,A7,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A6,1x,A6,1x,A8,1x,A6,1x,A8,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)"
    !print *,"hello1"
    write(10,fstr) "ITIME","FHOUR","ID","So(m/100km)","LAT(N)","LON(E)","SoFlag","Ro(km)","SR","BGo(m/100km)","BGo-lat(m/100km)","BGo-lon(m/100km)","ZMIN(m)","ZLAT(N)","ZLON(E)","ZFlag","Z850(m)","Z500(m)","Z200(m)","T850(K)","T500(K)","T200(K)","U850(m/s)","U500(m/s)","U200(m/s)","V850(m/s)","V500(m/s)","V200(m/s)","MR850(g/kg)","MR500(g/kg)","MR200(g/kg)","VDY(km)","VDX(km)","VDIST(km)","DT(h)","DUR(h)","MAXDUR(h)","PTY-OVR","FERRY(km)","FERRX(km)","FERR(km)","T(0)/M(1)/N(-1)","VLat(N)","VLon(E)","VSo","VRo","VZmin"
    !print *,"hello2"
    fstr = "(A11,1x,A4,1x,I9,1x,F6.2,1x,F6.2,1x,F6.2,1x,I3,1x,F10.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F6.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,I3,1x,F8.2,1x,F8.2,1x,F8.2,1x,F5.1,1x,F5.1,1x,F5.1,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.2,1x,F7.4,1x,F7.4,1x,F7.4,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,I4,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2,1x,F8.2)"
    !print *,"hello3"
    do s=1,nsfmax
      !print *,"help"
      if(itime(f,t,s).eq."-9999") exit  
      !print *,"hello4"
      write(10,fstr) itime(f,t,s), fhour(f,t,s), ID(f,t,s), So(f,t,s), Slat(f,t,s), Slon(f,t,s), isFMaxSo(f,t,s), Ro(f,t,s), SR(f,t,s), BGo(f,t,s), BGoy(f,t,s), BGox(f,t,s), Zmin(f,t,s), Zlat(f,t,s), Zlon(f,t,s), isFMinZ(f,t,s), z850(f,t,s), z500(f,t,s), z200(f,t,s), t850(f,t,s), t500(f,t,s), t200(f,t,s), u850(f,t,s), u500(f,t,s), u200(f,t,s), v850(f,t,s), v500(f,t,s), v200(f,t,s), mr850(f,t,s), mr500(f,t,s), mr200(f,t,s), saveDY(f,t,s), saveDX(f,t,s), saveDIST(f,t,s), saveDT(f,t,s),  saveDUR(f,t,s),maxID(f,t,s),saveFTARG(f,t,s),FERRY(f,t,s),FERRX(f,t,s),FERR(f,t,s),TrackMatch(f,t,s),saveVLat(f,t,s),saveVLon(f,t,s),saveVSo(f,t,s),saveVRo(f,t,s),saveVZmin(f,t,s)
      !print *,"hello5"
    end do
    close(10)
  end do
end do
call cpu_time(FN)
print '("Time to write forecast match files = ",f8.3," seconds.")',FN-ST


end program match_forecast


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
	    overlap = 0.0		!At most overlap at a single point
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
