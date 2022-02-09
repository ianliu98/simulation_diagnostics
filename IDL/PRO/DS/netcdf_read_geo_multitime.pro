;+
; NAME:
;    netcdf_read_geo_multitime.pro
;
; PURPOSE:
;   This function reads and returns the desired geographical data variable 
;    from a collection of NetCDF files which differ only in the time period 
;    included.  It also returns dimension variables and can do some data 
;    manipulation, thus acting as a driver of netcdf_read.pro for geographical 
;    data.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    var_data = netcdf_read_geo_multitime( file_name, var_label )
;
; INPUTS:
;    FILE_NAME:  A required string containing the name(s) of the NetCDF 
;        files to read.  This can be scalar or vector, its entries can 
;        include the '*' wildcard, and/or its can be directories (with a 
;        wildcard at the end '.../*').
;    VAR_LABEL:  An optional string scalar containing the label of the variable 
;        to read.  If no value is input then only dimension variables are 
;        returned.
;    CALENDAR, FILE_NETCDF_READ_GEO_VARINFO, GLOBAL_LABEL, LABEL_IN_HEIGHT, 
;      LABEL_IN_LAT, LABEL_IN_LON, LABEL_IN_MISSING, LABEL_IN_REALIZATION, 
;      LABEL_IN_TIME, LABEL_IN_VAR, MASK_LON, MASK_LAT, ORIGIN_TIME, 
;      PERIOD_DATA, PERIOD_REF, REGION, SHIFT_TIME, SUBSTITUTE
;
; KEYWORD PARAMETERS:
;    ANOMALY:  If set then anomalies from the PERIOD_REF period are returned.  
;        The default is not to calculate anomalies.
;    CALENDAR:  An optional scalar string naming the calendar type used in the 
;        TIME dimension.  The default is determined from the TIME variable's 
;        attributes.
;    FILE_NETCDF_READ_GEO_VARINFO:  An optional scalar string containing the 
;        name of the file, including directory, containing information about 
;        each variable, including standard units and mapping to alternate 
;        labels.  The default is netcdf_read_geo_varinfo.txt, which can be 
;        found automatically by the code if it is within the $IDL_PATH 
;        directory tree.
;    FIX_TIME:  If set, then the routine will attempt to produce a CF-standard 
;        time vector.  This is useful for instance in situations where monthly 
;        mean data has a time vector that is only approximately in the middle 
;        of the month, instead of exactly in the middle.  The default is to 
;        return the time values as they are in the input files.
;    GLOBAL_LABEL:  An optional vector string containing the labels of the 
;        global attributes which should have their values returned in 
;        GLOBAL_VALUE.
;    GLOBAL_TYPE:  If GLOBAL_LABEL is input, then this returns the variable 
;        type of the values returned in GLOBAL_VALUE, for instance "7" for 
;        a string.  Of same size as GLOBAL_LABEL, and of type integer.
;    GLOBAL_VALUE:  If GLOBAL_LABEL is input, then this returns the values of 
;        the global attributes requested in GLOBAL_LABEL.  Of same size as 
;        GLOBAL_LABEL.
;    HEIGHT:  Returns a vector containing the vertical coordinate variable.  Of 
;        length N_HEIGHT.  If the output field variable does not use a vertical 
;        dimension, then N_HEIGHT=1 and HEIGHT returns NaN.
;    INTEGRATE:  An optional string expression providing an instruction for 
;        integrating across dimensions.  See process_lonlatmonth.pro for 
;        details.  For example, 'integrate=1,2' integrates across the longitude 
;        and latitude dimensions, thus returning the global average (for each 
;        height level, time step, and realization).
;    LABEL_IN_HEIGHT:  An optional scalar string containing possible labels for 
;        the vertical dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'height,plev'.  The default is to use the 
;        list defined in netcdf_read_geo_varinfo.txt.
;    LABEL_IN_LAT:  An optional scalar string containing possible labels for 
;        the latitude dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'lat,latitude'.  The default is to use 
;        the list defined in netcdf_read_geo_varinfo.txt.
;    LABEL_IN_LON:  An optional scalar string containing possible labels for 
;        the longitude dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'lon,longitude'.  The default is to use 
;        the list defined in netcdf_read_geo_varinfo.txt.
;    LABEL_IN_MISSING:  An optional scalar string containing the label of the 
;        missing value attribute for the field variable.
;    LABEL_IN_REALIZATION:  An optional scalar string containing possible 
;        labels for the realization dimension and variable that may be used in 
;        the file.  If more that one label is provided they should be given in 
;        a comma-delimited format.  The default is to assume no realization 
;        dimension or variable.
;    LABEL_IN_TIME:  An optional scalar string containing possible labels for 
;        the time dimension and variable that may be used in the file.  If more 
;        that one label is provided they should be given in a comma-delimited 
;        format, e.g. 'time,t'.  The default is to use the list defined in 
;        netcdf_read_geo_varinfo.txt.
;    LABEL_IN_VAR:  An optional scalar string containing possible labels for 
;        the field variable that may be used in the file.  If more that one 
;        label is provided they should be given in a comma-delimited format, 
;        e.g. 'tas,TREFHT,T2m'.  The default is to use the list defined in 
;        netcdf_read_geo_varinfo.txt, or if such a list is absent then the 
;        standard variable label as provided in VAR_LABEL.
;    NO_UNITS_CONVERSION:  If set, then the function does not try to ensure 
;        that all variables, including dimensions, are in standard units.  The 
;        default is to check and try to convert if necessary.
;    LAT:  Returns a vector containing the latitude coordinate variable.  Of 
;        length N_LAT.  If the output field variable does not use a latitude 
;        dimension, then N_LAT=1 and LAT returns NaN.
;    LON:  Returns a vector containing the longitude coordinate variable.  Of 
;        length N_LON.  If the output field variable does not use a longitude 
;        dimension, then N_LON=1 and LON returns NaN.
;    MASK_INTERPOLATE_FRAC:  An optional floating point scalar defining the 
;        minimum fraction of the area in a new MASK_LON-MASK_LAT grid box 
;        required to have data in the original LON-LAT grid.  If the condition 
;        is not satisfied then the new value is defined as missing.  The 
;        default is 0.5.
;    MASK_LAT:  An optional floating point vector of latitude values for a new 
;        grid onto which to interpolate the data.
;    MASK_LON:  An optional floating point vector of longitude values for a new 
;        grid onto which to interpolate the data.
;    ORIGIN_TIME:  An optional scalar string containing the time origin value 
;        to use when converting the values of the time dimension to the 
;        "days since ORIGIN_TIME" format.  The default is to obtain the origin 
;        value from the file's metadata.  This should be in one of the 
;        following formats:  "yyyy-mm-dd", "yyyy-mm-dd-hh", 
;        "yyyy-mm-dd-hh-mm", or "yyyy-mm-dd-hh-mm-ss".
;    NO_SHRINK:  If set, then extraction of regional data using the REGION 
;        keyword input is returned on the original grid, but with areas outside 
;        the region(s) specified in REGION being assigned NaN values.  The 
;        default is to restrict the spatial extent of the VAR_DATA output 
;        to include the smallest possible area fully including the regions(s) 
;        specified in REGION.  This has no effect if REGION is not input.
;    NO_TIME_CHECK:  If set, no check is performed to ensure that the time 
;        dimension is not missing steps.  The default is to perform the check, 
;        which assumes that the data is regularly sampled in time.
;    NO_UNITS_CONVERSION:  If set, then the function does not try to ensure 
;        that all variables, including dimensions, are in standard units.  The 
;        default is to check and try to convert if necessary.
;    PREMASK:  On optional array of the same type as the VAR_DATA output array 
;        containing a data mask on the same longitude-latitude-height 
;        dimensions as the input NetCDF data.  It can either be of size 
;        N_LON*N_LAT*N_HEIGHT or N_LON*N_LAT*N_HEIGHT*N_TIME.  This is useful 
;        for instance for applying a land-sea mask.
;    QUIET:  If set then non-critical messages are not reported.  The default 
;        is to report all messages.
;    REALIZATION:  Returns a vector containing the realization coordinate 
;        variable.  Of length N_REALIZATION.  If the output field variable does 
;        not use a realization dimension, then N_REALIZATION=1 and REALIZATION 
;        returns NaN.
;    REF_PERIOD:  An optional string vector containing the end-points of the 
;        reference time period to use for calculating anomalies.  The default 
;        is TIME_PERIOD.  See TIME_PERIOD for the format to use.  Currently 
;        REF_PERIOD must fit fully within TIME_PERIOD.
;    REF_FRAC:  An optional floating point scalar defining the minimum fraction 
;        of the REF_PERIOD period required to have data when calculating 
;        anomalies.  If there are fewer valid data values available for a given 
;        longitude-latitude-height-realization, then all time steps are set to 
;        the missing value.  The default is 0.5.
;    REGION:  An optional input describing regional extraction of data to 
;        be performed.  See extract_region.pro for acceptable input formats.
;    QUIET:  If set then non-critical messages are not reported.  The default 
;        is to report all messages.
;    SHIFT_TIME:  A scalar number defining a forward (negative for backward) 
;        shift in time values.  This can be useful for instance if the time 
;        stamp at the end of the interval was used for mean values, instead of 
;        the standard middle of the interval.
;    SUBSTITUTE:  An optional 2*N_SUBSTITUTE array listing values to be 
;        replaced, and values to replace with, in the output VAR_DATA data 
;        array.  The value in element [0,i] is replaced with the value in 
;        [1,i].  This must be of the same type as VAR_DATA.
;    TIME:  Returns a vector containing the time coordinate variable.  Of 
;        length N_TIME.  If the output field variable does not use a time 
;        dimension, then N_TIME=1 and TIME returns NaN.
;    PERIOD_TIME:  An optional two-element string vector containing the 
;        end-points of the time period for which to return data.  The format 
;        should be in YYYYMMDDHHMMSS, although right-hand components (e.g. 
;        'SS') can be omitted.  YYYY is required as a minimum.  If right-hand 
;        components are omitted, it is assumed that the earliest and latest 
;        time, respectively, within the specified value applies.  For instance, 
;        if ['2016','2016'] is entered, then this is interpreted as 00:00:00 on 
;        1 January 2016 through 23:59:59 on 31 December 2016.  If available 
;        data in the files does not fully cover the period, then TIME and 
;        VAR_DATA are expanded such that they do fully cover the period, with 
;        the padded VAR_DATA elements being assigned NaN values.  If no value 
;        is input, then TIME_PERIOD is defined by the earliest and latest time 
;        values in the available data.  Note this is ignored if the time 
;        dimension is of length 1, because we then have no way of determining 
;        the time step.
;    UNITS_HEIGHT:  Returns a string containing the units of the height 
;        variable.
;    UNITS_LAT:  Returns a string containing the units of the lat variable.
;    UNITS_LON:  Returns a string containing the units of the lon variable.
;    UNITS_REALIZATION:  Returns a string containing the units of the 
;        realization variable.
;    UNITS_TIME:  An optional string containing the units attribute for the 
;        time dimension.  Also returns a string containing the attribute.
;    UNITS_VAR:  Returns a string containing the units of the field variable 
;        returned in VAR_DATA.
;
; OUTPUTS:
;    VAR_DATA:  Returns an array containing the data from the requested 
;        variable in the requested files.  Of size 
;        N_LON*N_LAT*N_HEIGHT*N_TIME*N_REALIZATION.
;    GLOBAL_TYPE, GLOBAL_VALUE, HEIGHT, LAT, LON,, REALIZATON, TIME, UNITS_VAR, 
;      UNITS_HEIGHT, UNITS_LAT, UNITS_LON, UNITS_REALIZATION, UNITS_TIME
;
; USES:
;    ncdump
;    convert_time_format.pro
;    mask_lonlattime.pro
;    month_day.pro
;    netcdf_read_geo.pro
;    process_lonlatmonth.pro
;    str.pro
;    var_type.pro
;
; PROCEDURE:
;    This function finds available NetCDF files according to the request, 
;    determines which file covers which time interval, and iterates through 
;    netcdf_read_geo.pro to load data from each file into the final output 
;    variable.
;
; EXAMPLE:
;    var_data = netcdf_read_geo_multitime( 'rain_data_19*.nc', 'pr', lon=lon, 
;        lat=lat, time=time )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-04-25, as 
;        climate_setread.pro.
;    Modified:  DAS, 2017-10-10 (Branched from climate_setread.pro;  $
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2017-10-27 (Added GLOBAL_LABEL and GLOBAL_VALUE keyword 
;        parameters;  fixed bug when FILE_NAME is input as a comma-delimited 
;        scalar;  fixed bug when working with Mac systems;  ensured that 
;        UNITS_TIME is returned;  fixed issues restricting time with 
;        TIME_PERIOD is input;  added FIX_TIME option)
;    Modified:  DAS, 2017-11-24 (Allowed absence of VAR_LABEL input)
;    Modified:  DAS, 2017-11-30 (Allowed case of undefined standard dimensions 
;        in file)
;    Modified:  DAS, 2017-12-24 (Added GLOBAL_TYPE keyword)
;    Modified:  DAS, 2018-02-15 (Corrected invalid keyword in call to 
;        process_lonlatmonth.pro)
;    Modified:  DAS, 2018-03-06 (Corrected bug in time restriction enforcement 
;        when there are multiple input NetCDF files)
;    Modified:  DAS, 2018-03-15 (Fixed confusion over time vector length when 
;        PERIOD_TIME is not defined)
;    Modified:  DAS, 2018-05-07 (Fixed bug in which PREMASK was ignored)
;    Modified:  DAS, 2018-09-13 (Added some error messages when code stops;  
;        Corrected error in retrieval of global attribute values)
;    Modified:  DAS, 2018-09-18 (Added "noleap" to list of accepted calendars)
;    Modified:  DAS, 2018-10-11 (Implemented padding of the time dimension 
;        when PERIOD_TIME is specified;  Fixed crashes when time dimension is 
;        of length one)
;    Modified:  DAS, 2018-10-16 (Fixed bug in padding at start of time period)
;    Modified:  DAS, 2018-11-07 (Fixed failure of MASK_LON and MASK_LAT 
;        interpolation to be performed)
;    Modified:  DAS, 2018-11-12 (Modified the spacers in the default time 
;        origin definition to 'yyyy-mm-dd dd:mm:ss')
;    Modified:  DAS, 2018-11-15 (Removed some ambiguity in time calendar search)
;    Modified:  DAS, 2018-11-28 (Switched order of time sorting and checking;  
;        Removed forcing of first file's calendar on other files)
;    Modified:  DAS, 2018-11-29 (Completed padding at start of period for 
;        monthly data)
;-

;***********************************************************************

FUNCTION NETCDF_READ_GEO_MULTITIME, $
    FILE_NAME, $
    VAR_LABEL, $
    ANOMALY=anomaly, $
    CALENDAR=time_calendar, ORIGIN=time_origin_in, $
    FILE_NETCDF_READ_GEO_VARINFO=file_netcdf_read_geo_varinfo, $
    GLOBAL_LABEL=global_label, GLOBAL_TYPE=global_type, $
      GLOBAL_VALUE=global_value, $
    INTEGRATE=integrate, $
    LABEL_IN_HEIGHT=height_label_in, LABEL_IN_LAT=lat_label_in, $
        LABEL_IN_LON=lon_label_in, LABEL_IN_REALIZATION=realization_label_in, $
        LABEL_IN_TIME=time_label_in, LABEL_IN_VAR=var_label_in, $
    LABEL_IN_MISSING=missing_label_in, $
    MASK_LON=mask_lon, MASK_LAT=mask_lat, $
      MASK_INTERPOLATE_FRAC=MASK_interpolate_frac, $
    PERIOD_TIME=time_period, $
    PREMASK=premask_data, $
    REF_PERIOD=ref_period, REF_FRAC=ref_frac, $
    REGION=region, $
    SHIFT_TIME=time_shift, $
    SUBSTITUTE=substitute, $
    LON=lon_data, LAT=lat_data, HEIGHT=height_data, TIME=time_data, $
      REALIZATION=realization_data, $
    UNITS_VAR=var_units, UNITS_LON=lon_units, UNITS_LAT=lat_units, $
      UNITS_HEIGHT=height_units, UNITS_TIME=time_units, $
      UNITS_REALIZATION=units_realization, $
    FIX_TIME=time_fix_opt, $
    NO_SHRINK=no_shrink_opt, $
    NO_TIME_CHECK=no_time_check_opt, $
    NO_UNITS_CONVERSION=no_units_conversion_opt, $
    QUIET=quiet_opt

;***********************************************************************
; Constants and Variables

; Option for doing things quietly (no reporting to terminal)
quiet_opt = keyword_set( quiet_opt )

; Absolute constants
nan = !values.f_nan
; An environment variable correction when running spawn on Mac systems
spawn, 'echo $OSTYPE', os_type
if strpos( strlowcase( os_type ), 'darwin' ) ge 0 then begin
  spawn_str = 'DYLD_LIBRARY_PATH="" ; '
endif else begin
  spawn_str = ''
endelse

; Ensure file and variable label specifications have been provided
if not( keyword_set( file_name ) ) then stop
if not( keyword_set( var_label ) ) then var_label = ''

; Ensure output variables are cleared
if n_elements( var_data ) ne 0 then temp = temporary( var_data )
if n_elements( lon_data ) ne 0 then temp = temporary( lon_data )
if n_elements( lat_data ) ne 0 then temp = temporary( lat_data )
if n_elements( height_data ) ne 0 then temp = temporary( height_data )
if n_elements( time_data ) ne 0 then temp = temporary( time_data )
if n_elements( realization_data ) ne 0 then temp = temporary( realization_data )

; Find the files satisfying request
temp_file_name = strsplit( file_name[0], ',', extract=1 )
for i_file = 1, n_elements( file_name ) - 1 do begin
  temp_file_name = [ temp_file_name, $
      strsplit( file_name[i_file], ',', extract=1 ) ]
endfor
file_list = file_search( temp_file_name, count=n_file )
if n_file eq 0 then begin
  stop, 'Error netcdf_read_geo_multitime.pro:  No files found.'
endif

; Copy any input calendar
if keyword_set( time_calendar ) then time_calendar_in = time_calendar
; Determine calendar (assume first file represents all of the files)
if not( keyword_set( time_calendar ) ) then begin
  spawn, spawn_str + 'ncdump -h ' + file_list[0] + " | grep 'time:calendar '", $
      time_calendar
  if n_elements( time_calendar ) ne 1 then stop
  if strlen( time_calendar ) eq 0 then stop
  time_calendar = strsplit( time_calendar, '=', extract=1, count=n_temp )
  if n_temp ne 2 then stop
  time_calendar = strsplit( time_calendar[1], '"', extract=1, count=n_temp )
  if n_temp ne 3 then stop
  time_calendar = strtrim( time_calendar[1], 2 )
endif
; Convert calendar to standard format
time_calendar = strlowcase( time_calendar )
if max( time_calendar eq [ 'gregorian', 'proleptic_gregorian', 'standard' ] ) $
    eq 1 then begin
  time_calendar = 'gregorian'
endif else if max( time_calendar eq [ '365_day', 'noleap' ] ) eq 1 then begin
  time_calendar = '365_day'
endif else if time_calendar eq '360_day' then begin
  time_calendar = '360_day'
endif else begin
  stop, 'Error netcdf_read_geo_multitime.pro:  Calendar not recognised (' $
      + time_calendar + ').'
endelse

; Check that valid time_period entries have been input
if keyword_set( time_period ) then begin
  if n_elements( time_period ) ne 2 then stop
  if var_type( time_period ) ne 7 then time_period = string( time_period )
  if max( strlen( time_period[0] ) eq [ 4, 6, 8, 10, 12, 14 ] ) ne 1 then stop
  if max( strlen( time_period[1] ) eq [ 4, 6, 8, 10, 12, 14 ] ) ne 1 then stop
  ; Expand abbreviated entries
  if strlen( time_period[0] ) eq 4 then time_period[0] = time_period[0] + '01'
  if strlen( time_period[0] ) eq 6 then time_period[0] = time_period[0] + '01'
  if strlen( time_period[0] ) eq 8 then time_period[0] = time_period[0] + '00'
  if strlen( time_period[0] ) eq 10 then time_period[0] = time_period[0] + '00'
  if strlen( time_period[0] ) eq 12 then time_period[0] = time_period[0] + '00'
  if strlen( time_period[1] ) eq 4 then time_period[1] = time_period[1] + '12'
  if strlen( time_period[1] ) eq 6 then begin
    if max( time_calendar eq [ '365_day', 'gregorian' ] ) eq 1 then begin
      temp_month = fix( strmid( time_period[1], 4, 2 ) )
      temp_day = month_day( temp_month - 1 )
      temp_day = temp_day[1] - temp_day[0] + 1
      if ( temp_month eq 2 ) and ( time_calendar eq 'gregorian' ) then begin
        temp_year = fix( strmid( time_period[1], 0, 4 ) )
        if temp_year mod 4 eq 0 then begin
          if ( temp_year mod 100 ne 0 ) or ( temp_year mod 400 eq 0 ) then begin
            temp_day = temp_day + 1
          endif
        endif
      endif
      time_period[1] = time_period[1] + str( temp_day, length=2, filler='0' )
    endif else if time_calendar eq '360_day' then begin
      time_period[1] = time_period[1] + '30'
    endif else begin
      stop
    endelse
  endif
  if strlen( time_period[1] ) eq 8 then time_period[1] = time_period[1] + '23'
  if strlen( time_period[1] ) eq 10 then time_period[1] = time_period[1] + '59'
  if strlen( time_period[1] ) eq 12 then time_period[1] = time_period[1] + '59'
endif
; Check that valid ref_period entries have been input
if keyword_set( ref_period ) then begin
  if n_elements( ref_period ) ne 2 then stop
  if var_type( ref_period ) ne 7 then ref_period = string( ref_period )
  if max( strlen( ref_period[0] ) eq [ 4, 6, 8, 10, 12, 14 ] ) ne 1 then stop
  if max( strlen( ref_period[1] ) eq [ 4, 6, 8, 10, 12, 14 ] ) ne 1 then stop
  ; Expand abbreviated entries
  if strlen( ref_period[0] ) eq 4 then ref_period[0] = ref_period[0] + '01'
  if strlen( ref_period[0] ) eq 6 then ref_period[0] = ref_period[0] + '01'
  if strlen( ref_period[0] ) eq 8 then ref_period[0] = ref_period[0] + '00'
  if strlen( ref_period[0] ) eq 10 then ref_period[0] = ref_period[0] + '00'
  if strlen( ref_period[0] ) eq 12 then ref_period[0] = ref_period[0] + '00'
  if strlen( ref_period[1] ) eq 4 then ref_period[1] = ref_period[1] + '12'
  if strlen( ref_period[1] ) eq 6 then begin
    if max( time_calendar eq [ '365_day', 'gregorian' ] ) eq 1 then begin
      temp_month = fix( strmid( ref_period[1], 4, 2 ) )
      temp_day = month_day( temp_month - 1 )
      temp_day = temp_day[1] - temp_day[0] + 1
      if ( temp_month eq 2 ) and ( time_calendar eq 'gregorian' ) then begin
        temp_year = fix( strmid( ref_period[1], 0, 4 ) )
        if temp_year mod 4 eq 0 then begin
          if ( temp_year mod 100 ne 0 ) or ( temp_year mod 400 eq 0 ) then begin
            temp_day = temp_day + 1
          endif
        endif
      endif
      ref_period[1] = ref_period[1] + str( temp_day, length=2, filler='0' )
    endif else if time_calendar eq '360_day' then begin
      ref_period[1] = ref_period[1] + '30'
    endif else begin
      stop
    endelse
  endif
  if strlen( ref_period[1] ) eq 8 then ref_period[1] = ref_period[1] + '23'
  if strlen( ref_period[1] ) eq 10 then ref_period[1] = ref_period[1] + '59'
  if strlen( ref_period[1] ) eq 12 then ref_period[1] = ref_period[1] + '59'
; The default reference period
endif else if keyword_set( time_period ) then begin
  ref_period = time_period
endif

; Set default requirement for good data
if n_elements( ref_frac ) eq 0 then ref_frac = 0.5
if n_elements( mask_interpolate_frac ) eq 0 then mask_interpolate_frac = 0.5

; Set option for restricting data to time_period or determining time_period
if keyword_set( time_period ) then begin
  time_period_activity = 'restrict'
endif else begin
  time_period_activity = 'determine'
endelse
file_list_period = strarr( 2, n_file )

;; Confirm any premask looks okay
;if keyword_set( premask_data ) then begin
;  if n_elements( premask_data[*,0,0,0] ) ne n_lon then stop
;  if n_elements( premask_data[0,*,0,0] ) ne n_lat then stop
;  if keyword_set( include_height_opt ) then begin
;    if n_elements( premask_data[0,0,*,0] ) ne n_height then stop
;  endif
;endif

;***********************************************************************
; Determine files, time dimension, and metadata

; Initialise file counter
ctr_file = 0

; Iterate through available files
for i_file = 0, n_file - 1 do begin
  ; Read the time vector for this file
  temp_time_data = nan
  if keyword_set( time_calendar_in ) then begin
    temp_time_calendar = time_calendar_in
  endif else begin
    temp_time_calendar = ''
  endelse
  if keyword_set( time_units ) then begin
    temp_time_units = time_units
  endif else begin
    temp_time_units = ''
  endelse
  temp = netcdf_read_geo( file_list[i_file], '', time=temp_time_data, $
      calendar=temp_time_calendar, units_time=temp_time_units, $
      shift_time=time_shift, label_in_time=time_label_in, quiet=quiet_opt )
  temp_time_range = [ min( temp_time_data, max=temp ), temp ]
  file_list_period[*,i_file] = convert_time_format( temp_time_range, $
      temp_time_units, 'yyyymmddhhmmss', calendar=temp_time_calendar )
  temp_time_data_ymd = convert_time_format( temp_time_data, $
      temp_time_units, 'yyyymmddhhmmss', calendar=temp_time_calendar )
  ;if not( keyword_set( time_calendar ) ) then begin
  ;  time_calendar = temp_time_calendar
  ;endif else begin
  ;  if temp_time_calendar ne time_calendar then stop
  ;endelse
  ; Omit the file from the list if it is outside of the requested period
  if time_period_activity eq 'restrict' then begin
    id_time = where( ( temp_time_data_ymd ge time_period[0] ) $
        and ( temp_time_data_ymd le time_period[1] ), n_temp_time )
    if n_temp_time eq 0 then begin
      file_list[i_file] = ''
    endif else begin
      temp_time_data = temp_time_data[id_time]
      temp_time_data_ymd = temp_time_data_ymd[id_time]
    endelse
  ; Refine the determination of the period based on the addition of this file
  endif else if time_period_activity eq 'determine' then begin
    if keyword_set( time_period ) then begin
      if file_list_period[0,i_file] lt time_period[0] then begin
        time_period[0] = file_list_period[0,i_file]
      endif
      if file_list_period[1,i_file] gt time_period[1] then begin
        time_period[1] = file_list_period[1,i_file]
      endif
    endif else begin
      time_period = file_list_period[*,i_file]
    endelse
    n_temp_time = n_elements( temp_time_data )
  endif else begin
    stop
  endelse
  ; Check if this file overlaps with other files
  if i_file gt 0 then begin
    id = where( $
        ( file_list_period[0,0:i_file-1] le file_list_period[1,i_file] ) and $
        ( file_list_period[1,0:i_file-1] ge file_list_period[0,i_file] ), n_id )
    if n_id ne 0 then stop
  endif
  ; Add these time values to our time dimension vector and note which file is 
  ; responsible
  if file_list[i_file] ne '' then begin
    if n_elements( time_data ) eq 0 then begin
      time_data = temp_time_data_ymd
      time_fileindex = intarr( n_temp_time ) + ctr_file + 1
    endif else begin
      time_data = [ time_data, temp_time_data_ymd ]
      time_fileindex = [ time_fileindex, intarr( n_temp_time ) + ctr_file + 1 ]
    endelse
    ctr_file = ctr_file + 1
  endif
endfor

; Refine the file list
if time_period_activity eq 'restrict' then begin
  id = where( file_list ne '', n_file )
  if n_file eq 0 then begin
    stop, 'Error netcdf_read_geo_multitime:  ' $
        + 'No files have data within requested period.'
  endif
  file_list = file_list[id]
  file_list_period = file_list_period[*,id]
endif
if ctr_file ne n_file then stop

; Sort time vector
id = sort( time_data )
time_data = time_data[id]
time_fileindex = time_fileindex[id]
id = -1
n_time = n_elements( time_data )

; Convert time dimension to 'days since' format
if keyword_set( time_origin_in ) then begin
  time_origin = time_origin_in
endif else begin
  time_origin = strmid( time_period[0], 0, 4 ) + '-01-01 00:00:00'
endelse
if not( keyword_set( time_units ) ) then begin
  time_units = 'days since ' + time_origin
endif
time_data_ymd = time_data
time_data = convert_time_format( time_data, 'yyyymmddhhmmss', time_units, $
    calendar=time_calendar )
; Convert time_period to 'days since' format
time_period_since = convert_time_format( time_period, 'yyyymmddhhmmss', $
    time_units, calendar=time_calendar )

; Check for gaps in time
; (This assumes that the biggest discrepancy from exactly regular spacing 
; should be (1-28/31), i.e. the largest with monthly data and a submonthly 
; time_units.)
if not( keyword_set( no_time_check_opt ) ) and ( n_time gt 1 ) then begin
  temp_diff = time_data[1:n_time-1] - time_data[0:n_time-2]
  temp_diff_max = max( temp_diff, min=temp_diff_min )
  if temp_diff_max / temp_diff_min gt 1.15 then stop
  if temp_diff_min eq 0 then stop
  temp_diff = 0
endif
; Restrict time vector's file link to values within a specified period
if time_period_activity eq 'restrict' then begin
  ; Identify values within specified period
  id = where( ( time_data lt time_period_since[0] ) $
      or ( time_data gt time_period_since[1] ), n_id )
  if n_id gt 1 then time_fileindex[id] = -1 * time_fileindex[id]
  if n_time gt 1 then begin
    temp = time_data[1] - time_data[0]
  endif else begin
    ; We cannot determine the type step with one time value, so this acts as a 
    ; switch to turn off padding
    temp = time_data[0] - time_period_since[0] + 1
  endelse
  ; Add padding to beginning if necessary
  if ( time_data[0] - time_period_since[0] ) gt temp then begin
    ; Work in "days since time_period_since[0]" format
    temp_time_data = convert_time_format( time_data, time_units, $
        'days since ' + time_origin, calendar=time_calendar )
    temp_time_start = convert_time_format( time_period_since[0], time_units, $
        'days since ' + time_origin, calendar=time_calendar )
    temp_time_data = temp_time_data - temp_time_start[0]
    d_temp_time_data = temp_time_data[1:n_time-1] - temp_time_data[0:n_time-2]
    ; If this is daily data
    if ( max( d_temp_time_data ) lt 1.1 ) $
        and ( min( d_temp_time_data ) gt 0.9 ) then begin
      ; Generate padding days
      temp_time_data = temp_time_data[0] $
          - reverse( 1 + findgen( floor( temp_time_data[0] ) ) ) $
          + temp_time_start[0]
      temp_time_data = convert_time_format( temp_time_data, $
          'days since ' + time_origin, time_units, calendar=time_calendar )
stop
;???Need to confirm this works
    ; If this is monthly data
    endif else if ( max( d_temp_time_data ) lt 32 ) $
        and ( min( d_temp_time_data ) gt 27 ) then begin
      ; Generate padding days
      n_temp_time = 12 * ( long( strmid( time_data_ymd[0], 0, 4 ) ) $
          - long( strmid( time_period[0], 0, 4 ) ) )
      n_temp = long( strmid( time_data_ymd[0], 4, 2 ) ) $
          - long( strmid( time_period[0], 4, 2 ) )
      n_temp_time = n_temp_time + n_temp
      temp_time_data = lonarr( n_temp_time + 1 )
      temp_time_data = long( strmid( time_period[0], 4, 2 ) ) $
          + lindgen( n_temp_time + 1 )
      temp_time_data = ( long( strmid( time_period[0], 0, 4 ) ) $
          + ( temp_time_data - 1 ) / 12l ) * 100l $
          + ( ( temp_time_data - 1 ) mod 12 ) + 1
      temp_time_data = str( temp_time_data, length=6, filler='0' ) + '01'
      temp_time_data = convert_time_format( temp_time_data, 'yyyymmdd', $
          time_units, calendar=time_calendar )
      ; If the timing is at the beginning of the month
      temp = fix( strmid( time_data_ymd, 6, 2 ) )
      if max( temp ) eq 1 then begin
        temp_time_data = temp_time_data[0:n_temp_time-1]
      ; If the timing is at the middle of the month
      endif else if ( min( temp ) ge 15 ) and ( max( temp ) le 16 ) then begin
        temp_time_data = ( temp_time_data[0:n_temp_time-1] $
          + temp_time_data[1:n_temp_time] ) / 2.
      ; Otherwise we are lost
      endif else begin
        stop
      endelse
    ; Otherwise this is not yet implemented
    endif else begin
      stop
    endelse
    ; Added padding to time vector
    time_data = [ temp_time_data, time_data ]
    time_fileindex = [ -1 + intarr( n_elements( temp_time_data ) ), $
        time_fileindex ]
    n_time = n_elements( time_data )
  endif
  ; Add padding to end if necessary
  if n_time gt 1 then begin
    temp = time_data[n_time-1] - time_data[n_time-2]
  endif else begin
    ; We cannot determine the type step with one time value, so this acts as a 
    ; switch to turn off padding
    temp = time_period_since[1] - time_data[n_time-1] + 1
  endelse
  if ( time_period_since[1] - time_data[n_time-1] ) gt temp then begin
    ; Work in "days since time_data[n_time-1]" format
    temp_time_units = convert_time_format( time_data[n_time-1], time_units, $
        'yyyymmddhhmmss', calendar=time_calendar )
    temp_time_units = 'days since ' + strmid( temp_time_units, 0, 4 ) + '-' $
        + strmid( temp_time_units, 4, 2 ) + '-' $
        + strmid( temp_time_units, 6, 2 ) + '-' $
        + strmid( temp_time_units, 8, 2 ) + '-' $
        + strmid( temp_time_units, 10, 2 ) + '-' $
        + strmid( temp_time_units, 12, 2 )
    temp_time_data = convert_time_format( time_data, time_units, $
        temp_time_units, calendar=time_calendar )
    temp_time_start = convert_time_format( time_period_since[0], time_units, $
        temp_time_units, calendar=time_calendar )
    temp_time_data = temp_time_data - temp_time_start[0]
    d_temp_time_data = temp_time_data[1:n_time-1] - temp_time_data[0:n_time-2]
    ; If this is daily data
    if ( max( d_temp_time_data ) lt 1.1 ) $
        and ( min( d_temp_time_data ) gt 0.9 ) then begin
      ; Generate padding days
      temp_time_data = temp_time_data[0] + 1 $
          + findgen( floor( temp_time_data[0] ) ) + temp_time_start[0]
      temp_time_data = convert_time_format( temp_time_data, $
          temp_time_units, time_units, calendar=time_calendar )
stop
;???Need to confirm this works
    ; If this is monthly data
    endif else if ( max( d_temp_time_data ) lt 32 ) $
        and ( min( d_temp_time_data ) gt 27 ) then begin
      ; Generate padding days
      n_temp_time = 12 * ( long( strmid( time_period[1], 0, 4 ) ) $
          - long( strmid( time_data_ymd[n_time-1], 0, 4 ) ) )
      n_temp = long( strmid( time_period[0], 4, 2 ) ) $
          - long( strmid( time_data_ymd[0], 4, 2 ) )
      n_temp_time = n_temp_time + n_temp
      temp_time_data = lonarr( n_temp_time + 1 )
      temp_time_data = long( strmid( time_data_ymd[n_time-1], 4, 2 ) ) $
          + lindgen( n_temp_time + 1 ) + 1
      temp_time_data = ( long( strmid( time_data_ymd[n_time-1], 0, 4 ) ) $
          + ( temp_time_data - 1 ) / 12l ) * 100l $
          + ( ( temp_time_data - 1 ) mod 12 ) + 1
      temp_time_data = str( temp_time_data, length=6, filler='0' ) + '01'
      temp_time_data = convert_time_format( temp_time_data, 'yyyymmdd', $
          time_units, calendar=time_calendar )
      ; If the timing is at the beginning of the month
      temp = strmid( time_data_ymd, 6, 2 )
      if max( temp ) eq 1 then begin
        temp_tie_data = temp_time_data[0:n_temp_time-1]
      ; If the timing is at the middle of the month
      endif else if ( min( temp ) ge 15 ) and ( max( temp ) le 16 ) then begin
        temp_time_data = ( temp_time_data[0:n_temp_time-1] $
            + temp_time_data[1:n_temp_time] ) / 2.
      ; Otherwise we are lost
      endif else begin
        stop
      endelse
    ; Otherwise this is not yet implemented
    endif else begin
      stop
    endelse
    ; Added padding to time vector
    time_data = [ time_data, temp_time_data ]
    time_fileindex = [ time_fileindex, $
        -1 + intarr( n_elements( temp_time_data ) ) ]
    n_time = n_elements( time_data )
  endif
endif

;***********************************************************************
; Load Data

; Initialise global attribute value vector
n_global = n_elements( global_label )
if n_global gt 0 then begin
  global_value = strarr( n_global )
  global_type = 7 + intarr( n_global )
endif

; Iterate through files
if var_label eq '' then begin
  n_temp_file = 1
endif else begin
  n_temp_file = n_file
endelse
for i_file = 0, n_temp_file - 1 do begin

  ; Copy any metadata for this segment
  if keyword_set( time_calendar_in ) then begin
    temp_time_calendar = time_calendar_in
  endif else begin
    temp_time_calendar = ''
  endelse

  ; Identify global attributes to read
  if n_global gt 0 then begin
    id_global = where( global_value eq '', n_id_global )
    if n_id_global eq 0 then begin
      temp_global_label = ''
    endif else begin
      temp_global_label = global_label[id_global]
    endelse
    temp_global_value = ''
    temp_global_type = 0
  endif else begin
    n_id_global = 0
  endelse
  ; Clear inputs for the next file read
  temp_var_units = ''
  temp_lon_units = ''
  temp_lat_units = ''
  temp_height_units = ''
  temp_time_units = ''
  temp_realization_units = ''
  if keyword_set( time_origin_in ) then begin
    temp_time_origin = time_origin_in
  endif else begin
    temp_time_origin = ''
  endelse
  temp_time_data = nan
  ; Load data from file
  temp_var_data = netcdf_read_geo( file_list[i_file], var_label, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo, $
      calendar=temp_time_calendar, origin_time=temp_time_origin, $
      label_in_var=var_label_in, label_in_lon=lon_label_in, $
      label_in_lat=lat_label_in, label_in_height=height_label_in, $
      label_in_time=time_label_in, label_in_realization=realization_label_in, $
      label_in_missing=missing_label_in, shift_time=time_shift, $
      region=region, no_shrink=no_shrink_opt, $
      no_units_conversion=no_units_conversion_opt, quiet=quiet_opt, $
      lon=temp_lon_data, lat=temp_lat_data, height=temp_height_data, $
      time=temp_time_data, realization=temp_realization_data, $
      units_var=temp_var_units, units_lon=temp_lon_units, $
      units_lat=temp_lat_units, units_height=temp_height_units, $
      units_time=temp_time_units, units_realization=temp_realization_units, $
      global_label=temp_global_label, global_type=temp_global_type, $
      global_value=temp_global_value )
  ; If this is the first file then copy non-time dimensions
  if i_file eq 0 then begin
    lon_data = temp_lon_data
    n_lon = n_elements( lon_data )
    lon_data_raw = lon_data
    n_lon_raw = n_lon
    lat_data = temp_lat_data
    n_lat = n_elements( lat_data )
    lat_data_raw = lat_data
    n_lat_raw = n_lat
    height_data = temp_height_data
    n_height = n_elements( height_data )
    height_data_raw = height_data
    n_height_raw = n_height
    realization_data = temp_realization_data
    n_realization = n_elements( realization_data )
    n_realization_raw = n_realization
    var_units = temp_var_units
    lon_units = temp_lon_units
    lat_units = temp_lat_units
    height_units = temp_height_units
    realization_units = temp_realization_units
  ; Otherwise check that non-time dimensions and units are consistent with 
  ; previous files
  endif else begin
    if n_elements( temp_lon_data ) ne n_lon_raw then stop
    if total( abs( finite( temp_lon_data ) - finite( lon_data_raw ) ) ) ne 0 $
        then stop
    id = where( finite( temp_lon_data ) eq 1, n_id )
    if n_id gt 0 then begin
      if total( abs( temp_lon_data[id] - lon_data_raw[id] ) ) ne 0 then stop
    endif
    if n_elements( temp_lat_data ) ne n_lat_raw then stop
    if total( abs( finite( temp_lat_data ) - finite( lat_data_raw ) ) ) ne 0 $
        then stop
    id = where( finite( temp_lat_data ) eq 1, n_id )
    if n_id gt 0 then begin
      ;if total( abs( temp_lat_data[id] - lat_data_raw[id] ) ) ne 0 then stop
      if total( abs( temp_lat_data[id] - lat_data_raw[id] ) ) gt 0.0002 $
          then stop
    endif
    if n_elements( temp_height_data ) ne n_height_raw then stop
    if total( abs( finite( temp_height_data ) - finite( height_data_raw ) ) ) $
        ne 0 then stop
    id = where( finite( temp_height_data ) eq 1, n_id )
    if n_id gt 0 then begin
      if total( abs( temp_height_data[id] - height_data_raw[id] ) ) ne 0 $
          then stop
    endif
    if n_elements( temp_realization_data ) ne n_realization_raw then stop
    if keyword_set( realization_data_raw ) then begin
      for i_realization = 0, n_realization - 1 do begin
        if temp_realization_data[i_realization] $
            ne realization_data_raw[i_realization] then stop
      endfor
    endif
    if temp_var_units ne var_units then stop
    if temp_lon_units ne lon_units then stop
    if temp_lat_units ne lat_units then stop
    if temp_height_units ne height_units then stop
    if keyword_set( realization_data_raw ) then begin
      if temp_realization_units ne realization_units then stop
    endif
  endelse

  ; Restrict to requested period
  temp_time_data = convert_time_format( temp_time_data, temp_time_units, $
      'yyyymmddhhmmss', calendar=temp_time_calendar )
  temp_time_data = convert_time_format( temp_time_data, 'yyyymmddhhmmss', $
      time_units, calendar=time_calendar )
  if time_period_activity eq 'restrict' then begin
    ;id_time = where( ( temp_time_data ge time_period_since[0] ) $
    ;    and ( temp_time_data le time_period_since[1] ), n_temp_time )
    id_time = where( abs( time_fileindex ) - 1 eq i_file, n_temp_time )
    if n_temp_time eq 0 then stop
    id = where( time_fileindex[id_time] - 1 eq i_file, n_temp_time )
    if n_temp_time eq 0 then stop
    id_time = id_time[id]
    temp = isin( time_data[id_time], temp_time_data )
    id = where( temp eq 1, n_id )
    if n_id ne n_temp_time then stop
    temp_time_data = temp_time_data[id]
    temp_var_data = temp_var_data[*,*,*,id,*]
  endif else begin
    n_temp_time = n_elements( temp_time_data )
  endelse

  ; Option to substitute one value for another
  if keyword_set( substitute ) then begin
    if n_elements( substitute[*,0] ) ne 2 then stop
    for i_substitute = 0, n_elements( substitute[0,*] ) - 1 do begin
      id = where( temp_var_data eq substitute[0,i_substitute], n_id )
      if n_id gt 0 then temp_var_data[id] = substitute[1,i_substitute]
    endfor
  endif

  ; Apply a pre-mask if requested (e.g. a land mask) which has the same 
  ; longitude-latitude-height grid as the source data
  if keyword_set( premask_data ) then begin
    ; If the premask is time invariant
    if n_elements( premask_data[0,0,0,*,0] ) eq 1 then begin
      ; Apply the premask
      for i_realization = 0, n_realization_raw - 1 do begin
        for i_time = 0, n_temp_time - 1 do begin
          temp_var_data[*,*,*,i_time,i_realization] $
              = temp_var_data[*,*,*,i_time,i_realization] * premask_data
        endfor
      endfor
    ; If the premask is time varying
    endif else begin
      ; Apply the premask
      for i_realization = 0, n_realization_raw - 1 do begin
        temp_var_data[*,*,*,*,i_realization] $
            = temp_var_data[*,*,*,*,i_realization] * premask_data
      endfor
    endelse
  endif

  ; Interpolate onto a different grid
  if keyword_set( mask_lon ) or keyword_set( mask_lat ) then begin
    ; Interpolate to grid (area-weighted scheme)
    temp_var_data = reform( temp_var_data, n_lon_raw, n_lat_raw, $
        n_height_raw * n_temp_time * n_realization_raw )
    mask_lonlattime, temp_var_data, lat=temp_lat_data, lon=temp_lon_data, $
        mask_lon=mask_lon, mask_lat=mask_lat, $
        frac_interpolate_thresh=mask_interpolate_frac
    if i_file eq 0 then begin
      lon_data = temp_lon_data
      n_lon = n_elements( lon_data )
      lat_data = temp_lat_data
      n_lat = n_elements( lat_data )
    endif
    temp_var_data = reform( temp_var_data, n_lon, n_lat, n_height_raw, $
        n_temp_time, n_realization_raw )
  endif

  ; Integrate along specified dimensions
  if keyword_set( integrate ) then begin
    ; Process the data
    process_lonlatmonth, temp_var_data, lat=temp_lat_data, $
        integrate=integrate, include_height=1
    if i_file eq 0 then begin
      if ( strpos( integrate, '=1' ) gt 0 ) $
          or ( strpos( integrate, ',1' ) gt 0 ) then begin
        n_lon = 1
        lon_data = 0
      endif
      if ( strpos( integrate, '=2' ) gt 0 ) $
          or ( strpos( integrate, ',2' ) gt 0 ) then begin
        n_lat = 1
        lat_data = 0
      endif
      if ( strpos( integrate, '=3' ) gt 0 ) $
          or ( strpos( integrate, ',3' ) gt 0 ) then begin
        n_height = 1
        height_data = 0
      endif
      if ( strpos( integrate, '=5' ) gt 0 ) $
          or ( strpos( integrate, ',5' ) gt 0 ) then begin
        n_realization = 1
        realization_data = ''
      endif
    endif
  endif

  ; Initialise variable's data array
  if ( i_file eq 0 ) and ( var_label ne '' ) then begin
    temp_type = var_type( temp_var_data )
    if temp_type eq 1 then begin
      var_data = bytarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 2 then begin
      var_data = intarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 3 then begin
      var_data = lonarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 4 then begin
      var_data = nan * fltarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 5 then begin
      var_data = nan * dblarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 6 then begin
      var_data = nan $
          * complexarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 7 then begin
      var_data = strarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else if temp_type eq 9 then begin
      var_data = nan $
          * dcomplexarr( n_lon, n_lat, n_height, n_time, n_realization )
    endif else begin
      stop
    endelse
  endif

  ; Copy data for the variable
  if var_label ne '' then begin
    ; If time values for this file are consecutive then copy the loaded data 
    ; into the data array in one go (greater efficiency)
    id_time = where( time_fileindex - 1 eq i_file, n_id_time )
    if n_id_time eq 0 then stop
    check_time = 0
    if n_id_time eq 1 then begin
      check_time = 1
    endif else begin
      temp = time_fileindex[id_time[1:n_id_time-1]] $
          - time_fileindex[id_time[0:n_id_time-2]]
      if max( abs( temp ) ) eq 0 then check_time = 1
      temp = temp_time_data[1:n_temp_time-1] - temp_time_data[0:n_temp_time-2]
      if min( temp ) lt 0 then check_time = 0
    endelse
    if check_time eq 1 then begin
      var_data[*,*,*,id_time[0]:id_time[n_id_time-1],*] = reform( $
          temp_var_data, n_lon, n_lat, n_height, n_temp_time, n_realization )
    ; Otherwise copy the loaded data iteratively
    endif else begin
      for i_time = 0, n_temp_time - 1 do begin
        id_time = where( time_data eq temp_time_data[i_time], n_id_time )
        if n_id_time ne 1 then stop
        var_data[*,*,*,id_time[0],*] = reform( temp_var_data[*,*,*,i_time,*], $
            n_lon, n_lat, n_height, 1, n_realization )
      endfor
    endelse
  endif

  ; Copy global attributes
  if n_id_global gt 0 then begin
    global_value[id_global] = temp_global_value
    global_type[id_global] = temp_global_type
  endif

  ; Print status
  if quiet_opt eq 0 then print, 'Loaded ' + file_list[i_file]

endfor

;; Restrict to period
;n_time = n_elements( time_data )
;id_time = where( time_fileindex gt 0, n_id_time )
;if n_time eq 0 then stop
;if n_id_time ne n_time then begin
;  time_data = time_data[id_time]
;  if keyword_set( var_data ) then var_data = var_data[*,*,*,id_time,*]
;  n_time = n_id_time
;endif

; Check if time is ascending
if n_time gt 1 then begin
  temp = min( time_data[1:n_time-1] - time_data[0:n_time-2] )
  ; If it is not ascending then fix it
  if temp lt 0 then begin
    ; Put the data into ascending order
    id = sort( time_data )
    time_data = time_data[id]
    if keyword_set( var_data ) then var_data = var_data[*,*,*,id,*]
  endif
endif
; Attempt to adjust the time vector to CF standard
if keyword_set( time_fix_opt ) then begin
  ; If this is monthly data (we need at least two time steps to determine this)
  if ( strpos( time_units, 'days since ' ) eq 0 ) and ( n_time gt 1 ) then begin
    temp_diff = time_data[1:n_time-1] - time_data[0:n_time-2]
    if ( min( temp_diff ) ge 25 ) and ( max( temp_diff ) le 35 ) then begin
      temp_diff = 0
      ; If the time vector is in the middle of the month
      temp = convert_time_format( time_data[0], time_units, 'yyyymmdd', $
          calendar=time_calendar )
      temp = fix( strmid( temp, 6, 2 ) )
      if ( temp ge 10 ) and ( temp le 20 ) then begin
        temp_time_ym = convert_time_format( time_data, time_units, $
            'yyyymm', calendar=time_calendar )
        temp_time_start = convert_time_format( temp_time_ym + '01000000', $
            'yyyymmddhhmmss', time_units, calendar=time_calendar )
        temp_time_end = str( long( temp_time_ym ) + 1, length=6, filler='0' )
        id = where( strmid( temp_time_end, 4, 2 ) eq '13', n_id )
        if n_id gt 0 then begin
          temp_time_end[id] $
              = str( long( strmid( temp_time_end[id], 0, 4 ) ) + 1, length=4, $
              filler='0' ) $
              + '01'
        endif
        temp_time_end = convert_time_format( temp_time_end + '01000000', $
            'yyyymmddhhmmss', time_units, calendar=time_calendar )
        time_data = ( temp_time_start + temp_time_end ) / 2.
      endif
    endif
  endif
endif

; Take the anomaly
if keyword_set( anomaly ) then begin
  ; Process the data
  process_lonlatmonth, var_data, time=time_data, anomaly=anomaly, $
      period_base=ref_period, frac_base_thresh=ref_frac, $
      include_height=1
endif

;***********************************************************************
; The End

if not( keyword_set( var_data ) ) then var_data = ''
return, var_data
END
