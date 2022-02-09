;+
; NAME:
;    c20c_dtos_v2_adjust_sic_stonepall2018
;
; PURPOSE:
;    This procedure returns the sea ice concentration versus sea surface 
;    temperature function used in Stone and Pall (2018) and used in the C20C+ 
;    D&A project Nat-Hist/CMIP5-est1 estimate of a world without human 
;    interference with the climate system.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_adjust_sic_stonepall2018, sic_file=sic_file, tos_file=tos_file, fit_method=fit_method, period=period 
;
; INPUTS:
;    FIT_METHOD, OCEAN_MASK_DATA, OCEAN_MASK_THRESH, PERIOD, SIC_FILE, TOS_FILE
;      TOS_FREEZE
;
; KEYWORD PARAMETERS:
;    FIT_METHOD:  A required input string vector providing instructions for the 
;        method of estimating the function.  One of the following must be 
;        included to specify how data should be pooled along the longitude 
;        dimension and should be returned:
;        * 'lon global':  Along each latitude division (see below), data is 
;          pooled from all ocean grid cells and a single global function is 
;          returned.  This setting was used for the C20C+ D&A project's 
;          Nat-Hist/CMIP5-est1 sea ice coverage estimate.
;        * 'lon hemispheres':  Along each latitude division (see below), data 
;          is pooled separately for the Eastern and Western Hemispheres, and a 
;          function is returned separately for each hemisphere.
;        * 'lon grid':  Along each latitude division (see below), data is 
;          pooled separately for each longitude in the input data and a 
;          separate function is returned for each longitude division.
;        One of the following must be included to specify how data should be 
;        pooled along the latitude dimension and should be returned:
;        * 'lat global':  Along each longitude division (see above), data is 
;          pooled from all ocean grid cells and a single global function is 
;          returned.
;        * 'lat hemispheres':  Along each longitude division (see above), data 
;          is pooled separately for the Northern and Southern Hemispheres, and 
;          a function is returned separately for each hemisphere.  This setting 
;          was used for the C20C+ D&A project's Nat-Hist/CMIP5-est1 sea ice 
;          coverage estimate.
;        * 'lat grid':  Along each longitude division (see above), data is 
;          pooled separately for each longitude in the input data and a 
;          separate function is returned for each latitude division.
;        One of the following must be included to specify how data should be 
;        pooled through the calendar year and should be returned:
;        * 'time annual':  Data from throughout the year is pooled and a single 
;          function is returned (for each longitude and latitude division).  
;          This setting was used for the C20C+ D&A project's 
;          Nat-Hist/CMIP5-est1 sea ice coverage estimate.
;        * 'time seasonal':  Data is pooled across years for each date within 
;          the year and a separate function is retuned for each date within the 
;          year.  For example, for monthly data only data from Januaries is 
;          used in calculating a function for January.
;        One of the following must be included to specify the type of function 
;        to calculate.
;        * 'linear anchored':  A linear fit is calculated, with one point set 
;          to the freezing point with full ice coverage, and the other 
;          determined by the average of the medians in various bins of sea ice 
;          coverage (with the average coverage across the bins being 
;          half-covered).  This setting was used for the C20C+ D&A project's 
;          Nat-Hist/CMIP5-est1 sea ice coverage estimate.
;        * 'linear tls':  A total least squares linear fit is calculated.
;        * 'bins':  The median sea surface temperature values for each of a 
;          number of bins of sea ice coverage are returned.  These are the same 
;          bins and medians used for the 'linear anchored' method.
;    FIT_LAT:  Returns a floating-point vector of length N_FIT_LAT containing 
;        the latitude dimension for FIT_SIC_DATA and FIT_TOS_DATA.
;    FIT_LON:  Returns a floating-point vector of length N_FIT_LON containing 
;        the longitude dimension for FIT_SIC_DATA and FIT_TOS_DATA.
;    FIT_SIC_DATA:  Returns a floating-point array containing the sea ice 
;        concentration data describing the function.  The array is of size 
;        N_FIT_POINT*N_FIT_LON*N_FIT_LAT*N_FIT_TIME, where the second, third, 
;        and fourth dimensions are described by FIT_LON, FIT_LAT, and FIT_TIME 
;        respectively.  The first dimension contains the number of sea surface 
;        temperature/sea ice concentration points described by this array and 
;        FIT_TOS_DATA.  In units of % (0 to 100).
;    FIT_TIME:  Returns a floating-point vector of length N_FIT_TIME containing 
;         the time dimension for FIT_SIC_DATA and FIT_TOS_DATA.  The format is 
;         'yyyymmdd'.
;    FIT_TOS_DATA:  Returns a floating-point array containing the sea surface 
;        temperature data describing the function.  The array is of size 
;        N_FIT_POINT*N_FIT_LON*N_FIT_LAT*N_FIT_TIME, where the second, third, 
;        and fourth dimensions are described by FIT_LON, FIT_LAT, and FIT_TIME 
;        respectively.  The first dimension contains the number of sea surface 
;        temperature/sea ice concentration points described by this array and 
;        FIT_SIC_DATA.  In units of Kelvin.
;    FIT_TYPE:  Returns a description for c20c_dtos_v2_adjust_sic.pro on how 
;        to interpret the data in FIT_TOS_DATA and FIT_SIC_DATA returned by 
;        this procedure.  Possible values are:
;        * 'straight line':  The two values for the first dimension in 
;          FIT_TOS_DATA and FIT_SIC_DATA describe the end-points of a straight 
;          full-ice to no-ice line.  This setting was used for the C20C+ D&A 
;          project's Nat-Hist/CMIP5-est1 sea ice coverage estimate.
;        * 'bins':  The first dimension in FIT_SIC_DATA is for the mean of the 
;          bins for sea ice concentration.  The first dimension in FIT_TOS_DATA 
;          is for the sea surface temperature values corresponding to those 
;          bins.
;    OCEAN_MASK_DATA:  An optional floating-point array defining the fraction 
;        the area of each grid cells that is ocean.  Of size N_LON*N_LAT, where 
;        these are the same sizes of the spatial dimensions of the data read 
;        from SIC_FILE and TOS_FILE.  Values should be 1 over ocean, 0 over 
;        land, and between 0 and 1 for partially-oceanic cells.
;    OCEAN_MASK_THRESH:  An optional floating-point scalar defining the 
;        threshold for the fraction of each cell that is ocean when determining 
;        a binary land-or-ocean distinction.  The default is 0.1.
;    PERIOD:  A required two-element string vector defining the beginning and 
;        end dates, respectively, of the data to use for calculating the 
;        functional relation.  Format is ['yyyymmdd','yyyymmdd'], with 
;        ['yyyymm','yyyymm'] also working for monthly data..
;    SIC_FILE:  A required string scalar or array listing the NetCDF files 
;        from which to read the sea ice concentration data which will be used 
;        to calculate the functional relation with sea surface temperature.  
;        This can include partial or full directory paths, and also can include 
;        wildcards.  The dimensions of the data must be the same as for the 
;        data specified in TOS_FILE.
;    TOS_FILE:  A required string scalar or array listing the NetCDF files 
;        from which to read the sea surface temperature data which will be used 
;        to calculate the functional relation with sea ice concentration.  
;        This can include partial or full directory paths, and also can include 
;        wildcards.  The dimensions of the data must be the same as for the 
;        data specified in SIC_FILE.
;    TOS_FREEZE:  An optional floating-point scalar defining the freezing point 
;        in Kelvin.  The default is 271.35K.
;    V1: If set then the routine reproduces the output of the v1-0 code.  The 
;        differences are that the v1-0 had bin_sample_min=1 and used a time-
;        varying land-sea mask, which for the NOAA-OI.v2 calculations was 
;        actually time-varying and based on the 2001-2010 period (note that 
;        this set only reproduces the method of using the time-varying mask, 
;        not the period of fit).
;
; OUTPUTS:
;    FIT_LAT, FIT_LON, FIT_SIC_DATA, FIT_TIME, FIT_TOS_DATA, FIT_TYPE
;
; USES:
;    c20c_dtos_v2_unnan.pro
;    convert_time_format.pro
;    netcdf_read_geo_multitime.pro
;    plus.pro
;    regtls.pro
;
; PROCEDURE:
;    This defines the sea ice adjustment function developed in:
;      * Stone, D. A., and P. Pall.  2018.  A benchmark estimate of the effect 
;        of anthropogenic emissions on the ocean surface.  Geosci. Model. Dev., 
;        submitted.
;
; EXAMPLES:
;    See c20c_dtos_v2_adjust_sic.pro.
;    This call should reproduce Nat-Hist/CMIP5-est1 calculation for the C20C+ 
;    D&A project, provided the appropriate sic_file and tos_file entries are 
;    provided and no ocean masking is required:
;      c20c_dtos_v2_adjust_sic, sic_file=sic_file, tos_file=tos_file, fit_method=['lon global','lat hemispheres','time annual','linear anchored'], period=['20010101','20101231'], fit_tos_data=fit_tos_data, fit_sic_data=fit_sic_data, fit_type=fit_type
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-06-18, as 
;        c20c_adjust_sic_pall.pro.
;    Modified:  DAS, 2017-10-13 (Extracted from c20c_adjust_sic_pall.pro into 
;        c20c_dtos_v2_adjust_sic_stonepall2018.pro;  added to IDL routine 
;        library)
;    Modified:  DAS, 2017-10-28 (Fixed problem with time vector comparison 
;        between tos and sic files)
;    Modified:  DAS, 2018-05-08 (Added V1 keyword option)
;    Modified:  DAS, 2018-08-17 (Added time-varying mask implementation to V1 
;        option;  Added V1 option implementation on bin selection)
;-

;***********************************************************************

PRO C20C_DTOS_V2_ADJUST_SIC_STONEPALL2018, $
    SIC_FILE=sic_file, $
    TOS_FILE=tos_file, $
    FIT_METHOD=fit_method, $
    OCEAN_MASK_DATA=ocean_mask_data, OCEAN_MASK_THRESH=ocean_mask_thresh, $
    PERIOD=period, $
    TOS_FREEZE=tos_freeze_0, $
    FIT_TOS_DATA=fit_tos_data, FIT_SIC_DATA=fit_sic_data, $
    FIT_LON=lon_out, FIT_LAT=lat_out, FIT_TIME=time_out, $
    FIT_TYPE=fit_type, $
    V1=v1_opt

;***********************************************************************
; Constants

; Option to reproduce v1 method
v1_opt = keyword_set( v1_opt )

; Ensure required inputs are entered
if not( keyword_set( tos_file ) ) then stop
if not( keyword_set( sic_file ) ) then stop

; The not-a-number flag
nan = !values.f_nan

; The freezing temperature
if n_elements( tos_freeze_0 ) eq 1 then tos_freeze = tos_freeze_0

; Define full and empty ice coverage
sic_full = 100.
sic_none = 0.

; The ice coverage histogram bin size
bin_sic_size = 1.0
; Define the bins
n_bin_sic = round( ( sic_full - sic_none ) / bin_sic_size ) + 1
bin_sic_min = findgen( n_bin_sic ) * bin_sic_size
bin_sic_max = bin_sic_min + bin_sic_size
bin_sic_mean = bin_sic_min + bin_sic_size / 2.
; The minimum number of samples in a bin in order to accept the calculation
if v1_opt eq 1 then begin
  bin_sample_min = 1
endif else begin
  bin_sample_min = 5
endelse

; The ocean fraction threshold for considering an area to be ocean
if not( keyword_set( ocean_mask_thresh ) ) then ocean_mask_thresh = 0.1
; Ensure that threshold is enforced on any input mask
if keyword_set( ocean_mask_data ) and ( v1_opt eq 0 ) then begin
  premask_data = float( plus( ocean_mask_data - ocean_mask_thresh ) )
  id = where( premask_data eq 0, n_id )
  if n_id gt 0 then premask_data[id] = nan
endif

;***********************************************************************
; Load fit data

; Load fit sea surface temperature data from file
tos_data = netcdf_read_geo_multitime( tos_file, 'tos', period_time=period, $
    lon=lon_data, lat=lat_data, time=time_data, premask=premask_data, $
    units_time=time_units, calendar=time_calendar, quiet=1 )
if not( keyword_set( tos_data ) ) then stop
n_lon = n_elements( lon_data )
n_lat = n_elements( lat_data )
n_time = n_elements( time_data )
tos_data = reform( tos_data, n_lon, n_lat, n_time )
; Convert time dimension to date format
time_data = convert_time_format( time_data, time_units, 'yyyymmdd', $
    calendar=time_calendar )

; Remove highly negative values in HadISST1
id = where( tos_data lt -500, n_id )
if n_id gt 1 then tos_data[id] = tos_freeze
; Ensure the data looks okay
if max( finite( tos_data ) ) eq 0 then stop
if min( tos_data, nan=1 ) lt -273.15 - 2. then stop
if max( tos_data, nan=1 ) gt 343.15 then stop
temp = min( tos_data, nan=1 )
if n_elements( tos_freeze ) eq 1 then begin
  if abs( tos_freeze - temp ) gt 0.01 then stop
endif else begin
  ; Take the lowest tos_data value as the default freezing point
  tos_freeze = temp
endelse

; Load fit SIC data
sic_data = netcdf_read_geo_multitime( sic_file, 'sic', period_time=period, $
    lon=temp_lon, lat=temp_lat, time=temp_time, premask=premask_data, quiet=1 )
if not( keyword_set( sic_data ) ) then stop
if n_elements( temp_lon ) ne n_lon then stop
if max( abs( temp_lon - lon_data ) ) gt 0 then stop
if n_elements( temp_lat ) ne n_lat then stop
if max( abs( temp_lat - lat_data ) ) gt 0 then stop
if n_elements( temp_time ) ne n_time then stop
temp_time = convert_time_format( temp_time, time_units, 'yyyymmdd', $
    calendar=time_calendar )
if max( abs( long( temp_time ) - long( time_data ) ) ) gt 0 then stop
sic_data = reform( sic_data, n_lon, n_lat, n_time )
; Ensure the data looks okay
if max( finite( sic_data ) ) eq 0 then stop
if min( sic_data ) lt -0.01 then stop
if min( sic_data ) gt 1.01 then stop

; Ensure that ocean mask is identical for both data sets
id = where( finite( sic_data ) + finite( tos_data ) eq 1, n_id )
if n_id gt 0 then begin
  sic_data[id] = !values.f_nan
  tos_data[id] = !values.f_nan
  if keyword_set( premask_data ) then begin
    id = where( finite( sic_data[*,*,0] ) + finite( premask_data ) eq 1, n_id )
    if n_id gt 0 then premask_data[id] = !values.f_nan
  endif
endif

; Confirm that approximately 30% of the values in either tos_data or 
; fit_data_sic are NaNs (presumably corresponding to land values)
id_tos = where( finite( tos_data ) eq 0, n_id_tos )
id_sic = where( finite( sic_data ) eq 0, n_id_sic )
if ( n_id_tos eq 0 ) and ( n_id_sic eq 0 ) then begin
  stop
endif else if ( n_id_tos eq 0 ) and ( n_id_sic ne 0 ) then begin
  tos_data[id_sic] = nan
  n_id_tos = n_id_sic
endif else if ( n_id_tos ne 0 ) and ( n_id_sic eq 0 ) then begin
  sic_data[id_tos] = nan
  n_id_sic = n_id_tos
endif else if n_id_tos ne n_id_sic then begin
  stop
endif else begin
  if max( abs( id_tos - id_sic ) ) ne 0 then stop
endelse
temp = float( n_id_tos ) / n_lon / n_lat / n_time
if ( temp lt 0.25 ) or ( temp gt 0.35 ) then stop
id_tos = -1
id_sic = -1

;***********************************************************************
; Determine output dimensions

; Determine the number of different longitudes for which to perform separate 
; fits
if max( fit_method eq 'lon global' ) eq 1 then begin
  n_lon_out = 1
  lon_out = 0.
endif else if max( fit_method eq 'lon hemispheres' ) eq 1 then begin
  n_lon_out = 2
  lon_out = [ -90., 90. ]
endif else if max( fit_method eq 'lon grid' ) eq 1 then begin
  n_lon_out = n_lon
  lon_out = lon_data
endif else begin
  stop
endelse
; Determine the number of different latitudes for which to perform separate fits
if max( fit_method eq 'lat global' ) eq 1 then begin
  n_lat_out = 1
  lat_out = 0.
endif else if max( fit_method eq 'lat hemispheres' ) eq 1 then begin
  n_lat_out = 2
  lat_out = [ -45., 45. ]
endif else if max( fit_method eq 'lat grid' ) eq 1 then begin
  n_lat_out = n_lat
  lat_out = lat_data
endif
; Determine the number of different dates within a year for which to perform 
; separate fits
if max( fit_method eq 'time annual' ) eq 1 then begin
  n_time_out = 1
  time_out = ''
endif else if max( fit_method, 'time seasonal' ) eq 1 then begin
  ; Identify the time of recurrence of the first mmdd string
  id = where( strmid( time_data, 4, 4 ) eq strmid( time_data[0], 4, 4 ), n_id )
  if n_id lt 2 then stop
  n_time_out = id[1]
  time_out = strmid( time_data[0:n_time_out-1], 4, 4 )
endif else begin
  stop
endelse
; Determine the number of fit values to return for each fit
if max( fit_method eq 'bins' ) eq 1 then begin
  n_point_out = n_bin_sic
endif else if ( max( fit_method eq 'linear anchored' ) eq 1 ) $
    or ( max( fit_method eq 'linear tls' ) eq 1 ) then begin
  n_point_out = 2
endif else begin
  stop
endelse

; Determine sea mask for the output data
if n_lon_out + n_lat_out gt 2 then begin
  if keyword_set( premask_data ) then begin
    ocean_mask_data_out = float( finite( premask_data ) )
  endif else begin
    ocean_mask_data_out = float( finite( tos_data[*,*,0] ) )
  endelse
  id = where( ocean_mask_data_out eq 0, n_id )
  if n_id eq 0 then stop
  if n_lon_out eq 1 then begin
    ocean_mask_data_out = reform( mean( ocean_mask_data_out, dimension=1 ), 1, $
        n_lat )
  endif else if n_lon_out eq 2 then begin
    ocean_mask_data_out[0,*] = reform( $
        mean( ocean_mask_data_out[0:n_lon/2-1,*], dimension=1 ), 1, n_lat )
    ocean_mask_data_out[1,*] = reform( $
        mean( ocean_mask_data_out[n_lon/2:n_lon-1,*], dimension=1 ), 1, n_lat )
    ocean_mask_data_out = ocean_mask_data_out[0:1,*]
  endif else if n_lon_out ne n_lon then begin
    stop
  endif
  if n_lat_out eq 1 then begin
    ocean_mask_data_out = reform( mean( ocean_mask_data_out, dimension=2 ), $
    n_lon_out, 1 )
  endif else if n_lat_out eq 2 then begin
    ocean_mask_data_out[*,0] = reform( $
        mean( ocean_mask_data_out[*,0:n_lat/2-1], dimension=2 ), n_lon_out, 1 )
    ocean_mask_data_out[*,1] = reform( $
        mean( ocean_mask_data_out[*,n_lat/2:n_lat-1], dimension=2 ), $
        n_lon_out, 1 )
    ocean_mask_data_out = ocean_mask_data_out[*,0:1]
  endif else if n_lat_out ne n_lat then begin
    stop
  endif
  ocean_mask_data_out = float( plus( ocean_mask_data_out - ocean_mask_thresh ) )
  id = where( ocean_mask_data_out eq 0, n_id )
  if n_id gt 0 then ocean_mask_data_out[id] = nan
endif

; Initialise arrays containing fit information
fit_sic_data = nan * fltarr( n_point_out, n_lon_out, n_lat_out, n_time_out )
fit_tos_data = nan * fltarr( n_point_out, n_lon_out, n_lat_out, n_time_out )

;***********************************************************************
; Determine functions in the partially-ice-covered regime

; Iterate through dates within the year
for i_time = 0, n_time_out - 1 do begin
  ; Initialise array containing mean bin sea surface temperature for this time 
  ; step
  hist_tos_median = nan * fltarr( n_bin_sic, n_lon_out, n_lat_out )
  ; Identify time elements to take
  if max( fit_method eq 'time annual' ) eq 1 then begin
    n_id_time = n_time
    id_time = indgen( n_id_time )
  endif else if max( fit_method eq 'time seasonal' ) eq 1 then begin
    id_time = where( strmid( time_data, 4, 4 ) eq time_out[i_time], n_id_time )
  endif else begin
    stop
  endelse
  ; Iterate through latitudes to output
  for i_lat = 0, n_lat_out - 1 do begin
    ; Identify latitudes to extract
    if n_lat_out eq 1 then begin
      id_lat = indgen( n_lat )
    endif else if n_lat_out eq 2 then begin
      id_lat = indgen( ( n_lat + i_lat ) / 2 ) + i_lat * n_lat / 2
    endif else if n_lat_out eq n_lat then begin
      id_lat = i_lat
    endif else begin
      stop
    endelse
    ; Iterate through longitudes to output
    for i_lon = 0, n_lon_out - 1 do begin
      ; Identify longitudes to extract
      if n_lon_out eq 1 then begin
        id_lon = indgen( n_lon )
      endif else if n_lon_out eq 2 then begin
        id_lon = indgen( ( n_lon + i_lon ) / 2 ) + i_lon * n_lon / 2
      endif else if n_lon_out eq n_lon then begin
        id_lon = i_lon
      endif else begin
        stop
      endelse
      ; Extract data to be sampled
      temp_tos_data = tos_data[*,*,id_time]
      temp_tos_data = temp_tos_data[*,id_lat,*]
      temp_tos_data = temp_tos_data[id_lon,*,*]
      temp_sic_data = sic_data[*,*,id_time]
      temp_sic_data = temp_sic_data[*,id_lat,*]
      temp_sic_data = temp_sic_data[id_lon,*,*]
      ; Calculate median sea surface temperature for each sea ice concentration 
      ; bin
      if max( finite( temp_tos_data ) ) eq 1 then begin
        id = where( finite( temp_sic_data ) eq 1, n_id )
        if n_id ge bin_sample_min then begin
          temp_tos_data = temp_tos_data[id]
          temp_sic_data = temp_sic_data[id]
          ; Note the "+0.001" term is required for numerical reasons
          hist_bin_index $
              = floor( ( temp_sic_data - sic_none ) / bin_sic_size + 0.001 )
          id = where( hist_bin_index eq n_bin_sic, n_id )
          if n_id gt 0 then hist_bin_index[id] = n_bin_sic - 1
          for i_bin = 1, n_bin_sic - 2 do begin
            id = where( hist_bin_index eq i_bin, n_id )
            if n_id ge bin_sample_min then begin
              hist_tos_median[i_bin,i_lon,i_lat] = median( temp_tos_data[id] )
            endif
          endfor
        endif
      endif
    endfor
  endfor
  ; Fill ocean bins with insufficient data with values from nearest neighbours
  if n_lon_out + n_lat_out gt 3 then begin
    for i_bin = 0, n_bin_sic - 1 do begin
      temp_hist = reform( hist_tos_median[i_bin,*,*], n_lon_out, n_lat_out )
      id = where( finite( temp_hist ) eq 0, n_id )
      if max( ocean_mask_data_out[id], nan=1 ) eq 1 then begin
        c20c_dtos_v2_unnan, temp_hist, wrap_x=1
        id = where( finite( ocean_mask_data_out ) eq 0, n_id )
        if n_id gt 0 then begin
          hist_tos_median[i_bin,*,*] = reform( $
              temp_hist * ocean_mask_data_out, 1, n_lon_out, n_lat_out )
        endif
      endif
    endfor
  endif
  ; If recording of the bins is requested
  if max( fit_method eq 'bins' ) eq 1 then begin
    ; Record values
    fit_tos_data[*,*,*,i_time] = hist_tos_median
  ; If a linear fit is requested
  endif else if ( max( fit_method eq 'linear anchored' ) eq 1 ) $
      or ( max( fit_method eq 'linear tls' ) eq 1 ) then begin
    ; Iterate through latitude
    for i_lat = 0, n_lat_out - 1 do begin
      ; Iterate through longitude
      for i_lon = 0, n_lon_out - 1 do begin
        ; If an unconstrained linear fit is requested
        if max( fit_method eq 'linear tls' ) eq 1 then begin
          ; Identify bins that have data, excluding the ice-free and ice-full 
          ; bins
          id = where( finite( hist_tos_median[1:n_bin_sic-2,i_lon,i_lat]  ) $
              eq 1, n_id )
          if n_id eq 0 then stop
          id = id + 1
          temp_tos = hist_tos_median[id,i_lon,i_lat]
          temp_sic = bin_sic_mean[id]
          ; Calculate the slope and intercept of the fit
          temp_fit_tos_slope = regtls( temp_tos - mean( temp_tos ), $
              temp_sic - mean( temp_sic ) )
          temp_fit_tos_intercept = mean( temp_sic ) $
              - temp_tos_fit_slope * mean( temp_tos )
        ; If a linear fit anchored to the freezing/full-coverage point is 
        ; requested
        endif else if max( fit_method eq 'linear anchored' ) eq 1 then begin
          ; Identify bins that have data, excluding the ice-free and ice-full 
          ; bins (except with the V1 option then include ice-full)
          if v1_opt eq 1 then begin
            id = where( finite( hist_tos_median[1:n_bin_sic-1,i_lon,i_lat]  ) $
                eq 1, n_id )
            id = id + 1
          endif else begin
            id = where( finite( hist_tos_median[1:n_bin_sic-2,i_lon,i_lat]  ) $
                eq 1, n_id )
            id = id + 1
          endelse
          if n_id eq 0 then stop
          temp_tos = hist_tos_median[id,i_lon,i_lat]
          temp_sic = bin_sic_mean[id]
          ; Calculate the slope and intercept of the fit
          temp_fit_tos_slope = ( mean( temp_sic ) - sic_full ) $
              / ( mean( temp_tos, nan=1 ) - tos_freeze )
          temp_fit_tos_intercept = ( ( mean( temp_tos, nan=1 ) * sic_full ) $
              - ( tos_freeze * mean( temp_sic ) ) ) $
              / ( mean( temp_tos, nan=1 ) - tos_freeze )
        endif else begin
          stop
        endelse
        ; Record values
        fit_sic_data[*,i_lon,i_lat,i_time] = [ sic_full, sic_none ]
        fit_tos_data[*,i_lon,i_lat,i_time] $
            = ( fit_sic_data[*,i_lon,i_lat,i_time] - temp_fit_tos_intercept ) $
            / temp_fit_tos_slope
      endfor
    endfor
  endif else begin
    stop
  endelse
endfor

;***********************************************************************
; Generate other output variables

; Define fit type
if ( max( fit_method eq 'linear anchored' ) eq 1 ) $
    or ( max( fit_method eq 'linear tls' ) eq 1 ) then begin
  fit_type = 'straight line'
endif else if max( fit_method eq 'bins' ) eq 1 then begin
  fit_type = 'bins'
  ; Define the bins
  fit_sic_data = bin_sic_mean
endif else begin
  stop
endelse

;***********************************************************************
; The end

;stop
return
END
