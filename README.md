# Identifying Cutoff Lows &amp; Troughs in GFS Forecasts

1. Identification and Tracking of Operational GFS Analysis (Truth) and Forecast Features
- Location of gridded operational GFS data:
  - `/glade/collections/rda/data/ds084.1`
- Location of text data identified cutoff lows and troughs: 
  - `/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb`
  - *version10 was a small test to compare “cutoffs” and “closed lows” using a criteria from Muñoz et al. (2020). Can ignore this directory.
- Location of “long lists” of text data organized by fhour:
  - `/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/longlists`
- Relevant scripts:
  - Location: `scripts/`
  - `compile.csh`
    - source this on casper to compile executables `identification_algorithm_globe` and `identification_algorithm_globe_cases`
  - `driver_IdentifyFeatures.csh`
    - Submit each month and fhour as a batch job to casper to identify features. The end result of running this script is a long list of .dat files (one for each initial time (4x daily) and forecast hour (41 x itimes; 0, 6, 12,…,240)
  - `run_IdentifyFeatures.csh`
    - Submitted by driver script, receives year, month, fhour information from the driver. 
    - For operational GFS data, this scripts 850, 500, and 200 hPa Z, T, U, V, and RH from grib2 files, converts to smaller netcdf files.
    - For each itime and fhour, the identification_algorithm_globe is submitted with command line arguments setting the output file, search radii, slope ratio threshold, and other parameters
  - `identification_algorithm_global_noDisambigSteps.f90`
    - This is the fortran code that is used to identify “cutoff lows” and “preexisting troughs”
    - Refer to comments in the code for documentation. Also see Lupo et al. (2023) and Kasuga et al. (2021)
  - `identification_algorithm_globe `
    - Compiled identification_algorithm_globe_ noDisambigSteps.f90. This code is run using a combination of driver_IdentifyFeatures.csh run_IdentifyFeatures.csh. 
  - `driver_TrackAnalysis.csh`
    - Driver script to run track_analysis (note that there is no “run_” script). 
    - Tracks analysis features in serial (from the first valid time to the last), run on an interactive casper node.
    - User can set normalization values for the penalty terms (some testing configurations are provided. The current configuration is “pmax1.5_2stdnorms_munozDmax1200_oppmax700” and is based on normalization values used by Lupo et al. (2023)
    - The list of analysis/verification time files is provided to the fortran code, which outputs .track files
  - `track_analysis.f90`
    - This is the fortran code that is used to track analysis-time features. See .f90 file for comments.  
  - `track_analysis`
    - Compiled `track_analysis.f90`. This code is run using driver_TrackAnalysis.csh
  - `driver_TrackForecast.csh`
    - Driver script to run `run_TrackForecast.csh` and `track_forecast.f90`
    - Submit batch jobs to track forecast features for from f000-f240 for each initial time during each month. Submit each month as a separate casper job
    - User can set normalization values for the penalty terms (some testing configurations are provided. The current configuration is “pmax1.5_2stdnorms_munozDmax1200_oppmax700” and is based on normalization values used by Lupo et al. (2023)
    - The list of analysis/verification time files is provided to the fortran code, which outputs .track files
  - `run_TrackForecast.csh`   
    - Submits the compiled fortran code track_forecast  with several command line arguments. 
    - All forecast files matching a use-specified file naming convention are added to a list, and the list is provided to track_forecast. 
    - The list of valid-time track files is subset to find the valid time files corresponding to the forecast hours being submitted to the Fortran code. This subset of the verification data is also submitted to track_forecast
    - The list of fx/verification files is provided to the fortran code, which outputs .track files
    - See fortran code for documentation
  - `utility_MergeLists.csh`
    - For cleaner organization and to more easily aggregate metrics for a specific forecast hour, this utility merges all track files valid at a specified fhhh (note that all header/column label information is also removed)
    - User sets which experiment is being merged, and which set of forecast hours the utility should iterate over.

- General order of operations:
  1. Run identify features code for analysis/verif data and forecast data (order does not matter here). This will generate .dat files of feature characteristics (and optional debugging netcdf files)
  2. Run “track_analysis” code. This will assign IDs to all verification features and add tracking information to new .track files (for the verification data, “error” columns are still included, but set to missing values for consistency with the forecast .track files). This code must be run prior to track_forecast, as track_forecast uses the ID information and feature characteristics from the verification dataset to match forecast and verification features (for computing forecast errors)
  3. Run “track_forecast” code. This will match forecast and verification features at the initial time, and then track forecast features through subsequent fhours. Features that are not tracked from the initial time are attempted to be matched to verification features at the corresponding valid time, and then are subsequently tracked.
  4. Merge lists. For ease of aggregating data to compute forecast errors as a function of forecast lead time.
- For GFS analyses and forecasts from Jan 2015–August 2022, the above steps are completed, with the output stored in `~klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb` and `~klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/longlists`

2. Identification and Tracking of Experimental UFS-MRW Forecast Features
- Location of verification/truth data:
   - `/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/*f000*`
-	Location of UFS model output data:
    - `/glade/campaign/mmm/parc/mwong/ufs-mrw/2019102206.F240.C768`
-	Relevant scripts:
   - Location: /glade/scratch/klupo/interp
   - interp_ufs_output.csh
     - Adapted version of Craig’s tile stitching code.
     - Note changes to “top_dir”, “storage”, and the command line arg “infhr” to allow these forecast hours to be stitched in parallel
     - The storage directory is later used by the feature ID/tracking code
  - Location: scripts/
  - `netcdf_routines_mod.mod`
    - Necessary script to read 2d lat/lon from gridded UFS data (user doesn’t need to modify)
  -	`netcdf_routines_mod.o`  
    -	Necessary script to read 2d lat/lon from gridded UFS data (user doesn’t need to modify)
  - `netcdf_routines_mod.f90`
    - Necessary script to read 2d lat/lon from gridded UFS data (user doesn’t need to modify)
  - `compile.csh`     
    - Run this script to compile identification_algorithm_global_noDisambigSteps_cases.f90 if modifications are made (compiles with netcdf routines needed to read lat/lon
  -	`driver_IdentifyFeatures_cases.csh`
    -	Submit each itime and fhour as a batch job to casper to identify features. The end result of running this script is a list of .dat files for each initial time and forecast hour (41 x itimes; 0, 6, 12,…,240)
  -	`run_IdentifyFeatures_cases.csh`
    -	Submitted by driver script, receives year, month, fhour information from the driver. 
    -	For each itime and fhour, the identification_algorithm_globe_cases is submitted with command line arguments setting the output file, search radii, slope ratio threshold, and other parameters
  -	`identification_algorithm_global_noDisambigSteps_cases.f90`
    -	This is the fortran code that is used to identify “cutoff lows” and “preexisting troughs”
    -	Refer to comments in the code for documentation. Also see Lupo et al. (2023) and Kasuga et al. (2021)
    -	Output of UFS feature scheme is swapped latitude order from GFS verif. Shouldn’t be a problem though, I don’t think the fortran code cares which order the features are read in.
  -	`identification_algorithm_globe_cases`
    -	Compiled `identification_algorithm_globe_noDisambigSteps_cases.f90`. This code is run using a combination of `driver_IdentifyFeatures_cases.csh` `run_IdentifyFeatures_cases.csh`. 
  -	`driver_TrackForecast_cases.csh`
    -	Driver script to run `run_TrackForecast_cases.csh` and `track_forecast_cases.f90` 
    -	Submit batch jobs to track forecast features for from f000-f240 for each initial time. Submit each itime as a separate casper job (or, could probably get away with running this on a login node..it’s fast and just using text data)
    -	User can set normalization values for the penalty terms (some testing configurations are provided. The current configuration is “pmax1.5_2stdnorms_munozDmax1200_oppmax700” and is based on normalization values used by Lupo et al. (2023)
    - The list of analysis/verification time files is provided to the fortran code, which outputs .track files
  -	`run_TrackForecast_cases.csh`
    -	Symlink verif track data (GFS analyses) to case directory. Also symlink the f000 forecast at the initial time as a substitute UFS f000 forecast (should be basically identical, but UFS doesn’t output data at f000, and fortran code expects the first dat file to be filled with analysis time features)
    -	Submits the compiled fortran code track_forecast_cases  with several command line arguments. 
    -	All forecast files matching a use-specified file naming convention are added to a list, and the list is provided to track_forecast. 
    -	The list of valid-time track files is subset to find the valid time files corresponding to the forecast hours being submitted to the Fortran code. This subset of the verification data is also submitted to track_forecast
    -	The list of fx/verification files is provided to the fortran code, which outputs .track files
    -	See fortran code for documentation
  -	`track_forecast_cases.f90`
    -	Fortran code to track/match forecast features. 
    -	Reads both UFS fx and GFS vx files
    -	See fortran code for documentation
  - `track_forecast_cases`
    -	Compiled version of the above code

- Description of output
  * YYYYMMDDhh.fhhh.nc: Stitched output files
  * diag_TroughsCutoffs_*.dat: Output of feature ID scheme
  * diag_TroughsCutoffs_*.track: Output of feature tracking scheme
  * diag_TroughsCutoffs_*.nc: Debugging output of feature ID scheme
  * gfs.0p25.*: Verification files from original cutoff low climo
 
 
- Output file format (.dat columns)
1. `ITIME` 
    *	Initial time of forecast
    *	YYYYMMDDhh
    *	Can treat as integer
2. `FHOUR`
    *	Forecast hour from ITIME
    *	Treat as string and subset when necessary
3. `So(m/100km)`
    *	Optimal slope of the identified feature with units of m/100km
    *	Treat as floating point value
4. `LAT(N)`
    *	Latitude of the feature (position of optimal slope), north is positive
    *	Treat as floating point value
5. `LON(E)`
    *	Longitude (degrees east) of the feature (position of the optimal slope)
    *	Treat as floating point value
6. `Ro(km)`
    *	Optimal radius of the feature (km)
    *	Treat as floating point value
7. `SR`
    *	Slope ratio; ratio of BGo to So
    *	This value is limited to a user specified value. Features with SR>thresh are excluded from output
    *	Treat as floating point value
8. `BGo(m/100km)`
    *	Background slope at the location of the identified feature
    *	Units of m/100km
    *	Treat as floating point value
9. `BGo-lat(m/100km)`
    *	Latitudinal component of background slope at the location of the identified feature
    *	Units of m/100km
    *	Treat as floating point value
10. `BGo-lon(m/100km)`
    *	Longitudinal component of background slope at the location of the identified feature
    *	Units of m/100km
    *	Treat as floating point value
11. `ZMIN(m)`
    *	Local geopotential height minimum within the radius of a cutoff low, units of m
    *	This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    *	Treat as floating point value
12. `ZLAT(N)`
    *	Latitude of the ZMIN associated with the feature, north is positive
    *	This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    *	Treat as floating point value
13. `ZLON(E)`
    *	Longitude (degrees east) of the ZMIN associated with the feature 
    *	This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    *	Treat as floating point value
14. `Z850(m)`
    *	Area-averaged 850hPa geopotential height over a feature’s Ro (units of m)
     *	Treat as floating point value
15. `Z500(m)`
    *	Area-averaged 500hPa geopotential height over a feature’s Ro (units of m)
    *	Treat as floating point value
16. `Z200(m)`
    *	Area-averaged 200hPa geopotential height over a feature’s Ro (units of m)
    *	Treat as floating point value
17. `T850(K)`
    *	Area-averaged 850hPa temperature over a feature’s Ro (units of K)
    *	Treat as floating point value
18. `T500(K)`
    *	Area-averaged 500hPa temperature over a feature’s Ro (units of K)
    *	Treat as floating point value
19. `T200(K)`
    *	Area-averaged 200hPa temperature over a feature’s Ro (units of K)
    *	Treat as floating point value
20. `U850(m/s)`
    *	Area-averaged 850hPa zonal wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
21. `U500(m/s)`
    *	Area-averaged 500hPa zonal wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
22. `U200(m/s)`
    *	Area-averaged 200hPa zonal wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
23. `V850(m/s)`
    *	Area-averaged 850hPa meridional wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
24. `V500(m/s)`
    *	Area-averaged 500hPa meridional wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
25. `V200(m/s)`
    *	Area-averaged 200hPa meridional wind over a feature’s Ro (units of m/s)
    *	Treat as floating point value
26. `MR850(g/kg)`
    *	Area-averaged 850hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    *	Treat as floating point value
27. `MR500(g/kg)`
    *	Area-averaged 500hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    *	Treat as floating point value
28. `MR200(g/kg)`
    *	Area-averaged 200hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    *	Treat as floating point value
29. `600kmZ500(m)`
    *	Area-averaged 500hPa geopotential height over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of m)
    *	Used for analysis-forecast matching
    *	Treat as floating point value
30. `600kmT500(K)`
    *	Area-averaged 500hPa temperature over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature)o (units of K)
    *	Used for analysis-forecast matching
    *	Treat as floating point value
31. `600kmU500(m/s)`
    * Area-averaged 500hPa zonal wind over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of m/s)
    *	Used for analysis-forecast matching
    *	Treat as floating point value
32. `600kmV500(m/s)`
    *	Area-averaged 500hPa meridional wind over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature)  (units of m/s)
    *	Used for analysis-forecast matching
    *	Treat as floating point value
33. `600kmMR500(g/kg)`
    *	Area-averaged 500hPa water vapor mixing ratio over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of g/kg)
    *	Used for analysis-forecast matching
    *	Treat as floating point value

Output file format (.track columns)
1. `ITIME` 
   * Initial time of forecast
   * YYYYMMDDhh
   * Can treat as integer
2. `ID`
   * ID of the forecast or verif feature. Forecast features receive their ID from a matched or tracked verification feature.
   * Values of ID= –1 indicate the forecast feature is not a match to an analysis time feature or a member of any forecast feature track
   * Treat as integer
3. `FHOUR`
   * Forecast hour from ITIME
   * Treat as string and subset when necessary
4. `So(m/100km)`
   * Optimal slope of the identified feature with units of m/100km
   * Treat as floating point value
5. `LAT(N)`
   * Latitude of the feature (position of optimal slope), north is positive
   * Treat as floating point value
6. `LON(E)`
   * Longitude (degrees east) of the feature (position of the optimal slope)
   * Treat as floating point value
7. `SoFlag`
   * Flag to indicate if this is the time of So-max for the verification feature
   * 0: Not the time
   * 1: Is the time of So max
   * –1: Fx feature has no corresponding verif feature
8. `Ro(km)`
   * Optimal radius of the feature (km)
   * Treat as floating point value
9. `SR`
   * Slope ratio; ratio of BGo to So
   * This value is limited to a user specified value. Features with SR&gt;thresh are excluded from output
   * Treat as floating point value
10. `BGo(m/100km)`
    * Background slope at the location of the identified feature
    * Units of m/100km
    * Treat as floating point value
11. `BGo-lat(m/100km)`
    * Latitudinal component of background slope at the location of the identified feature
    * Units of m/100km
    * Treat as floating point value
12. `BGo-lon(m/100km)`
    * Longitudinal component of background slope at the location of the identified feature
    * Units of m/100km
    * Treat as floating point value
13. `ZMIN(m)`
    * Local geopotential height minimum within the radius of a cutoff low, units of m
    * This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    * Treat as floating point value
14. `ZLAT(N)`
    * Latitude of the ZMIN associated with the feature, north is positive
    * This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    * Treat as floating point value
15. `ZLON(E)`
    * Longitude (degrees east) of the ZMIN associated with the feature 
    * This value will be set to -9999.9 if there is no local geopotential height min within a feature’s Ro (i.e., the feature is a trough)
    * Treat as floating point value
16. `ZFlag`
    * Flag to indicate if this is the time of Z-min for the verification feature
    * 0: Not the time
    * 1: Is the time of Z min
    * –1: Fx feature has no corresponding verif feature
17. `Z850(m)`
    * Area-averaged 850hPa geopotential height over a feature’s Ro (units of m)
    * Treat as floating point value
18. `Z500(m)`
    * Area-averaged 500hPa geopotential height over a feature’s Ro (units of m)
    * Treat as floating point value
19. `Z200(m)`
    * Area-averaged 200hPa geopotential height over a feature’s Ro (units of m)
    * Treat as floating point value
20. `T850(K)`
    * Area-averaged 850hPa temperature over a feature’s Ro (units of K)
    * Treat as floating point value
21. `T500(K)`
    * Area-averaged 500hPa temperature over a feature’s Ro (units of K)
    * Treat as floating point value
22. `T200(K)`
    * Area-averaged 200hPa temperature over a feature’s Ro (units of K)
    * Treat as floating point value
23. `U850(m/s)`
    * Area-averaged 850hPa zonal wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
24. `U500(m/s)`
    * Area-averaged 500hPa zonal wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
25. `U200(m/s)`
    * Area-averaged 200hPa zonal wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
26. `V850(m/s)`
    * Area-averaged 850hPa meridional wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
27. `V500(m/s)`
    * Area-averaged 500hPa meridional wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
28. `V200(m/s)`
    * Area-averaged 200hPa meridional wind over a feature’s Ro (units of m/s)
    * Treat as floating point value
29. `MR850(g/kg)`
    * Area-averaged 850hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    * Treat as floating point value
30. `MR500(g/kg)`
    * Area-averaged 500hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    * Treat as floating point value
31. `MR200(g/kg)`
    * Area-averaged 200hPa water vapor mixing ratio over a feature’s Ro (units of g/kg)
    * Treat as floating point value
32. `600kmZ500(m)`
    * Area-averaged 500hPa geopotential height over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of m)
    * Used for analysis-forecast matching
    * Treat as floating point value
33. `600kmT500(K)`
    * Area-averaged 500hPa temperature over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature)o (units of K)
    * Used for analysis-forecast matching
    * Treat as floating point value
34. `600kmU500(m/s)`
    * Area-averaged 500hPa zonal wind over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of m/s)
    * Used for analysis-forecast matching
    * Treat as floating point value
35. `600kmV500(m/s)`
    * Area-averaged 500hPa meridional wind over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature)  (units of m/s)
    * Used for analysis-forecast matching
    * Treat as floating point value
36. `600kmMR500(g/kg)`
    * Area-averaged 500hPa water vapor mixing ratio over a fixed 600km radius surrounding the location of a feature’s So (the location of the feature) (units of g/kg)
    * Used for analysis-forecast matching
    * Treat as floating point value
37. `DY(km)`
    * Meridional displacement from the previous location of the feature (units of km)
    * Treat as floating point value
38. `DX(km)`
    * Zonal displacement from the previous location of the feature (units of km)
    * Treat as floating point value
39. `DIST(km)`
    * Total displacement from the previous location of the feature (units of km)
    * Treat as floating point value
40. `DT(h)`
    * Time (hours) since the last time the feature was identified
    * Can be up to 12 hours, since a single missed step is allowed
    * Treat as floating point value
41. `DUR(h)`
    * Duration of the feature up to this valid time in forecasts/verif (hours)
    * Value is `0` if the forecast feature is not a match for any analysis feature
    * Treat as floating point value
42. `MAXDUR(h)`
    * Saved value referring to how long (hours) this feature track lasts for in the entire dataset
    * Can be used to determine if a feature is `non transient` e.g., part of a track that lasts for at least 36 hours
    * Value is -9999.9 if the forecast feature is not a match for any analysis feature
    * Treat as floating point value
43. `PTY-OVR`
    * Ignore the `-OVR` part of this variable name
    * For features that are tracked, this is the tracking penalty term
    * For features that are matched, this is the matching penalty term
    * 9999.9 for forecast features that do not correspond to verifying track
    * Treat as floating point value
44. `FERRY(km)`
    * Meridional displacement error (km) of the forecast feature
    * -99999.9 for forecast features that do not correspond to verifying track
    * Treat as floating point value
45. `FERRX(km)`
    * Zonal displacement error (km) of the forecast feature
    * -99999.9 for forecast features that do not correspond to verifying track
    * Treat as floating point value
46. `FERR(km)`
    * Total displacement error (km) of the forecast feature
    * -99999.9 for forecast features that do not correspond to verifying track
    * Treat as floating point value
47. `T(0)/M(1)/N`
    * Flag to determine if a forecast feature is a continuation of a forecast track (0), was matched to a verification feature (1), or does not correspond to a verification feature (-1)
48. `VLat(N)`
    * Latitude of the corresponding verification feature (position of optimal slope), north is positive
    * -9999.9 if no corresponding verification feature
    * Treat as floating point value
49. `VLon(E)`
    * Longitude (degrees east) of the corresponding verification feature (position of the optimal slope)
    * -9999.9 if no corresponding verification feature
    * Treat as floating point value
50. `VSo`
    * Optimal slope of the corresponding verification feature feature with units of m/100km
    * -9999.9 if no corresponding verification feature
    * Treat as floating point value
51. `VRo`
    * Optimal radius of the corresponding verification feature (km)
    * -9999.9 if no corresponding verification feature
    * Treat as floating point value
52. `VZmin`
    * Local height min of the corresponding verification feature (m)
    * -9999.9 if no corresponding verification feature (of if a Trough)
    * Treat as floating point value
