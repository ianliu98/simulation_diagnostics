;+
; NAME:
;    netcdf_read_geo.pro
;
; PURPOSE:
;    This function reads and returns the desired geographical data variable 
;    from a NetCDF file.  It also returns dimension variables and can do some 
;    data manipulation, thus acting as a driver of netcdf_read.pro for 
;    geographical data.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    var_data = netcdf_read_geo( file_name, var_label )
;
; INPUTS:
;    FILE_NAME:  A required scalar string containing the name of the NetCDF 
;        file.
;    VAR_LABEL:  A required scalar string containing the label of the variable 
;        to read.
;    ATTRIBUTE_HEIGHT_LABEL, ATTRIBUTE_LAT_LABEL, ATTRIBUTE_LON_LABEL, 
;      ATTRIBUTE_REALIZATION_LABEL, ATTRIBUTE_TIME_LABEL, ATTRIBUTE_VAR_LABEL, 
;      CALENDAR, FILE_NETCDF_READ_GEO_VARINFO, GLOBAL_LABEL, LABEL_IN_HEIGHT, 
;      LABEL_IN_LAT, LABEL_IN_LON, LABEL_IN_MISSING, LABEL_IN_REALIZATION, 
;      LABEL_IN_TIME, LABEL_IN_VAR, ORIGIN_TIME, REGION, SHIFT_TIME, UNITS_TIME
;
; KEYWORD PARAMETERS:
;    ATTRIBUTE_HEIGHT_LABEL:  An optional vector string containing the labels 
;        of the attributes for the height dimension variable which should have 
;        their values returned in ATTRIBUTE_HEIGHT_VALUE.  If 1 is input 
;        then all attribute values are returned in ATTRIBUTE_HEIGHT_VALUE, and 
;        all attribute labels are returned in this keyword.
;    ATTRIBUTE_HEIGHT_TYPE:  Returns the variable type of the height variable 
;        attributes returned in ATTRIBUTE_HEIGHT_VALUE.  Of length 
;        N_ATTRIBUTE_HEIGHT.
;    ATTRIBUTE_HEIGHT_VALUE:  Returns the values of all of the height variable 
;        attributes requested in ATTRIBUTE_HEIGHT_LABEL.  Of length 
;        N_ATTRIBUTE_HEIGHT.
;    ATTRIBUTE_LAT_LABEL:  An optional vector string containing the labels 
;        of the attributes for the latitude dimension variable which should 
;        have their values returned in ATTRIBUTE_LAT_VALUE.  If 1 is input 
;        then all attribute values are returned in ATTRIBUTE_LAT_VALUE, and 
;        all attribute labels are returned in this keyword.
;    ATTRIBUTE_LAT_TYPE:  Returns the variable type of the latitude variable 
;        attributes returned in ATTRIBUTE_LAT_VALUE.  Of length N_ATTRIBUTE_LAT.
;    ATTRIBUTE_LAT_VALUE:  Returns the values of all of the latitude variable 
;        attributes requested in ATTRIBUTE_LAT_LABEL.  Of length 
;        N_ATTRIBUTE_LAT.
;    ATTRIBUTE_LON_LABEL:  An optional vector string containing the labels 
;        of the attributes for the longitude dimension variable which should 
;        have their values returned in ATTRIBUTE_LON_VALUE.  If 1 is input 
;        then all attribute values are returned in ATTRIBUTE_LON_VALUE, and 
;        all attribute labels are returned in this keyword.
;    ATTRIBUTE_LON_TYPE:  Returns the variable type of the longitude variable 
;        attributes returned in ATTRIBUTE_LON_VALUE.  Of length N_ATTRIBUTE_LON.
;    ATTRIBUTE_LON_VALUE:  Returns the values of all of the longitude variable 
;        attributes requested in ATTRIBUTE_LON_LABEL.  Of length 
;        N_ATTRIBUTE_LON.
;    ATTRIBUTE_REALIZATION_LABEL:  An optional vector string containing the 
;        labels of the attributes for the realization dimension variable which 
;        should have their values returned in ATTRIBUTE_REALIZATION_VALUE.  If 
;        1 is input then all attribute values are returned in 
;        ATTRIBUTE_REALIZATION_VALUE, and all attribute labels are returned in 
;        this keyword.
;    ATTRIBUTE_REALIZATION_TYPE:  Returns the variable type of the realization 
;        variable attributes returned in ATTRIBUTE_REALIZATION_VALUE.  Of 
;        length N_ATTRIBUTE_REALIZATION.
;    ATTRIBUTE_REALIZATION_VALUE:  Returns the values of all of the realization 
;        variable attributes requested in ATTRIBUTE_REALIZATION_LABEL.  Of 
;        length N_ATTRIBUTE_REALIZATION.
;    ATTRIBUTE_TIME_LABEL:  An optional vector string containing the labels 
;        of the attributes for the time dimension variable which should have 
;        their values returned in ATTRIBUTE_TIME_VALUE.  If 1 is input 
;        then all attribute values are returned in ATTRIBUTE_TIME_VALUE, and 
;        all attribute labels are returned in this keyword.
;    ATTRIBUTE_TIME_TYPE:  Returns the variable type of the time variable 
;        attributes returned in ATTRIBUTE_TIME_VALUE.  Of length 
;        N_ATTRIBUTE_TIME.
;    ATTRIBUTE_TIME_VALUE:  Returns the values of all of the time variable 
;        attributes requested in ATTRIBUTE_TIME_LABEL.  Of length 
;        N_ATTRIBUTE_TIME.
;    ATTRIBUTE_VAR_LABEL:  An optional vector string containing the labels 
;        of the attributes for the data variable which should have their values 
;        returned in ATTRIBUTE_VAR_VALUE.  If 1 is input then all attribute 
;        values are returned in ATTRIBUTE_VAR_VALUE, and all attribute labels 
;        are returned in this keyword.
;    ATTRIBUTE_VAR_TYPE:  Returns the variable type of the data variable 
;        attributes returned in ATTRIBUTE_VAR_VALUE.  Of length N_ATTRIBUTE_VAR.
;    ATTRIBUTE_VAR_VALUE:  Returns the values of all of the data variable 
;        attributes requested in ATTRIBUTE_VAR_LABEL.  Of length 
;        N_ATTRIBUTE_VAR.
;    CALENDAR:  An optional scalar string naming the calendar type used in the 
;        TIME dimension.  The default is determined from the TIME variable's 
;        attributes.
;    FILE_NETCDF_READ_GEO_VARINFO:  An optional scalar string containing the 
;        name of the file, including directory, containing information about 
;        each variable, including standard units and mapping to alternate 
;        labels.  The default is netcdf_read_geo_varinfo.xml, which can be 
;        found automatically by the code if it is within the $IDL_PATH 
;        directory tree.
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
;    KEEP_TIME_FORMAT:  If set then the time dimension values are returned 
;        exactly as they are in the file.  The default is to convert, if 
;        necessary, to the "days since ORIGIN_TIME" format, where ORIGIN_TIME 
;        is specified as a keyword input or is obtained from the file's 
;        metadata.
;    LABEL_IN_HEIGHT:  An optional scalar string containing possible labels for 
;        the vertical dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'height,plev'.  The default is to use the 
;        list defined in netcdf_read_geo_varinfo.xml.
;    LABEL_IN_LAT:  An optional scalar string containing possible labels for 
;        the latitude dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'lat,latitude'.  The default is to use 
;        the list defined in netcdf_read_geo_varinfo.xml.
;    LABEL_IN_LON:  An optional scalar string containing possible labels for 
;        the longitude dimension and variable that may be used in the file.  If 
;        more that one label is provided they should be given in a 
;        comma-delimited format, e.g. 'lon,longitude'.  The default is to use 
;        the list defined in netcdf_read_geo_varinfo.xml.
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
;        netcdf_read_geo_varinfo.xml.
;    LABEL_IN_VAR:  An optional scalar string containing possible labels for 
;        the field variable that may be used in the file.  If more that one 
;        label is provided they should be given in a comma-delimited format, 
;        e.g. 'tas,TREFHT,T2m'.  The default is to use the list defined in 
;        netcdf_read_geo_varinfo.xml, or if such a list is absent then the 
;        standard variable label as provided in VAR_LABEL.
;    LAT:  Returns a vector containing the latitude coordinate variable.  Of 
;        length N_LAT.  If the output field variable does not use a latitude 
;        dimension, then N_LAT=1 and LAT returns NaN.
;    LON:  Returns a vector containing the longitude coordinate variable.  Of 
;        length N_LON.  If the output field variable does not use a longitude 
;        dimension, then N_LON=1 and LON returns NaN.
;    NO_SHRINK:  If set, then extraction of regional data using the REGION 
;        keyword input is returned on the original grid, but with areas outside 
;        the region(s) specified in REGION being assigned NaN values.  The 
;        default is to restrict the spatial extent of the VAR_DATA output 
;        to include the smallest possible area fully including the regions(s) 
;        specified in REGION.  This has no effect if REGION is not input.
;    NO_UNITS_CONVERSION:  If set, then the function does not try to ensure 
;        that all variables, including dimensions, are in standard units.  The 
;        default is to check and try to convert if necessary.
;    ORIGIN_TIME:  An optional scalar string containing the time origin value 
;        to use when converting the values of the time dimension to the 
;        "days since ORIGIN_TIME" format.  The default is to obtain the origin 
;        value from the file's metadata.  This should be in one of the 
;        following formats:  "yyyy-mm-dd", "yyyy-mm-dd-hh", 
;        "yyyy-mm-dd-hh-mm", or "yyyy-mm-dd-hh-mm-ss".
;    REALIZATION:  Returns a vector containing the realization coordinate 
;        variable.  Of length N_REALIZATION.  If the output field variable does 
;        not use a realization dimension, then N_REALIZATION=1 and REALIZATION 
;        returns NaN.
;    QUIET:  If set then non-critical messages are not reported.  The default 
;        is to report all messages.
;    REGION:  An optional input describing regional extraction of data to 
;        be performed.  See extract_region.pro for acceptable input formats.
;    SHIFT_TIME:  A scalar number defining a forward (negative for backward) 
;        shift in time values.  This can be useful for instance if the time 
;        stamp at the end of the interval was used for mean values, instead of 
;        the standard middle of the interval.
;    TIME:  Returns a vector containing the time coordinate variable.  Of 
;        length N_TIME.  If the output field variable does not use a time 
;        dimension, then N_TIME=1 and TIME returns NaN.
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
;        variable in the requested file.  Of size 
;        N_LON*N_LAT*N_HEIGHT*N_TIME*N_REALIZATION.
;    ATTRIBUTE_HEIGHT_LABEL, ATTRIBUTE_HEIGHT_TYPE, ATTRIBUTE_HEIGHT_VALUE, 
;      ATTRIBUTE_LAT_LABEL, ATTRIBUTE_LAT_TYPE, ATTRIBUTE_LAT_VALUE, 
;      ATTRIBUTE_LON_LABEL, ATTRIBUTE_LON_TYPE, ATTRIBUTE_LON_VALUE, 
;      ATTRIBUTE_REALIZATION_LABEL, ATTRIBUTE_REALIZATION_TYPE, 
;      ATTRIBUTE_REALIZATION_VALUE, ATTRIBUTE_TIME_LABEL, ATTRIBUTE_TIME_TYPE, 
;      ATTRIBUTE_TIME_VALUE, ATTRIBUTE_VAR_LABEL, ATTRIBUTE_VAR_TYPE, 
;      ATTRIBUTE_VAR_VALUE, GLOBAL_TYPE, GLOBAL_VALUE, HEIGHT, LAT, LON, 
;      REALIZATON, TIME, UNITS_VAR, UNITS_HEIGHT, UNITS_LAT, UNITS_LON, 
;      UNITS_REALIZATION, UNITS_TIME
;
; USES:
;    convert_time_format.pro
;    extract_region.pro
;    markup_read.pro
;    netcdf_read_geo_varinfo.xml
;    ncdump
;    netcdf_read.pro
;    str.pro
;    var_type.pro
;
; PROCEDURE:
;    This function uses netcdf_read.pro to read data from a NetCDF file.  It 
;    assumes that the file contains some sort of geographical data and thus has 
;    certain dimensions to return and other properties, thus acting as a driver 
;    for a full read operation on files satisfying those assumptions.
;
; EXAMPLE:
;    var_data = netcdf_read_geo( 'rain_data.nc', 'pr', lon=lon, lat=lat, $
;        time=time )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-04-28, as 
;        climate_fileread.pro.
;    Modified:  DAS, 2017-10-10 (Branched from climate_fileread.pro;  $
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2017-10-27 (Added GLOBAL_LABEL and GLOBAL_VALUE keyword 
;        parameters;  fixed bug when working with Mac systems;  added
;        capability to convert between 'fraction' and '%' units)
;    Modified:  DAS, 2017-11-08 (Fixed reading of dimensions of variable when 
;        variable label is not in file name;  fixed was procedure deals with 
;        dimensions not existing in file)
;    Modified:  DAS, 2017-12-24 (Added GLOBAL_TYPE keyword output)
;    Modified:  DAS, 2018-02-15 (Improved search for variable dimensions in 
;        NetCDF output to avoid some situations of confusion;  corrected bug 
;        in variable size check for very large data sets)
;    Modified:  DAS, 2018-03-11 (Added capability to handle character 
;        dimensions;  added capabilty to have different labels for realization 
;        dimension and dimension variable, implemented through a 
;        comma-delimited list)
;    Modified:  DAS, 2018-08-30 (Fixed variable name error in region 
;        extraction;  Implemented automatic retrieval of vertical dimension 
;        variables with name other than plev)
;    Modified:  DAS, 2018-09-18 (Added "1" to list of possible variable units 
;        equivalent to fraction)
;    Modified:  DAS, 2018-10-09 (Removed possibly ambiguous variable 
;        identification in dimension determination;  Added another possible
;        abbreviation of the units attribute for Celsius)
;    Modified:  DAS, 2018-10-16 (Assumed all but first entry are to be ignored 
;        when multiple instances of variable label found in file)
;    Modified:  DAS, 2018-11-26 (Ensured load of calendar specification when 
;        not input;  Added ATTRIBUTE_* keywords)
;-

;***********************************************************************

FUNCTION NETCDF_READ_GEO, $
    FILE_NAME, $
    VAR_LABEL, $
    ATTRIBUTE_HEIGHT_LABEL=height_attribute_label, $
      ATTRIBUTE_HEIGHT_VALUE=height_attribute_value, $
      ATTRIBUTE_HEIGHT_TYPE=height_attribute_type, $
    ATTRIBUTE_LAT_LABEL=lat_attribute_label, $
      ATTRIBUTE_LAT_VALUE=lat_attribute_value, $
      ATTRIBUTE_LAT_TYPE=lat_attribute_type, $
    ATTRIBUTE_LON_LABEL=lon_attribute_label, $
      ATTRIBUTE_LON_VALUE=lon_attribute_value, $
      ATTRIBUTE_LON_TYPE=lon_attribute_type, $
    ATTRIBUTE_REALIZATION_LABEL=realization_attribute_label, $
      ATTRIBUTE_REALIZATION_VALUE=realization_attribute_value, $
      ATTRIBUTE_REALIZATION_TYPE=realization_attribute_type, $
    ATTRIBUTE_TIME_LABEL=time_attribute_label, $
      ATTRIBUTE_TIME_VALUE=time_attribute_value, $
      ATTRIBUTE_TIME_TYPE=time_attribute_type, $
    ATTRIBUTE_VAR_LABEL=var_attribute_label, $
      ATTRIBUTE_VAR_VALUE=var_attribute_value, $
      ATTRIBUTE_VAR_TYPE=var_attribute_type, $
    CALENDAR=time_calendar, ORIGIN_TIME=time_origin, $
    FILE_NETCDF_READ_GEO_VARINFO=file_netcdf_read_geo_varinfo, $
    GLOBAL_LABEL=global_label, GLOBAL_TYPE=global_type, $
      GLOBAL_VALUE=global_value, $
    LABEL_IN_HEIGHT=height_label_in, LABEL_IN_LAT=lat_label_in, $
      LABEL_IN_LON=lon_label_in, LABEL_IN_REALIZATION=realization_label_in, $
      LABEL_IN_TIME=time_label_in, LABEL_IN_VAR=var_label_in, $
    LABEL_IN_MISSING=missing_label_in, $
    SHIFT_TIME=shift_time, $
    LON=lon_data, LAT=lat_data, HEIGHT=height_data, TIME=time_data, $
      REALIZATION=realization_data, $
    REGION=region, $
    UNITS_VAR=var_units, UNITS_LON=lon_units, UNITS_LAT=lat_units, $
      UNITS_HEIGHT=height_units, UNITS_TIME=time_units, $
      UNITS_REALIZATION=units_realization, $
    KEEP_TIME_FORMAT=keep_time_format_opt, $
    NO_UNITS_CONVERSION=no_units_conversion_opt, $
    NO_SHRINK=no_shrink_opt, $
    QUIET=quiet_opt


;***********************************************************************
; Constants

; Quiet option (no reporting to terminal)
quiet_opt = keyword_set( quiet_opt )
if quiet_opt eq 1 then !except = 0

; Default missing data flag
nan = !values.f_nan
; An environment variable correction when running spawn on Mac systems
spawn, 'echo $OSTYPE', os_type
if strpos( strlowcase( os_type ), 'darwin' ) ge 0 then begin
  spawn_str = 'DYLD_LIBRARY_PATH="" ; '
endif else begin
  spawn_str = ''
endelse

; Default missing value label
if not( keyword_set( missing_label_in ) ) then missing_label_in = '_FillValue'

; The default file containing information on variables
if n_elements( file_netcdf_read_geo_varinfo ) eq 0 then begin
  file_netcdf_read_geo_varinfo = file_which( 'netcdf_read_geo_varinfo.xml' )
endif

; Build the list of possible labels for the variable used in the file
if not( keyword_set( var_label_in ) ) and keyword_set( var_label ) then begin
  if file_netcdf_read_geo_varinfo ne '' then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers='label', select_values='label='+var_label, $
        settings=var_label_in
    var_label_in = strsplit( var_label_in, ',', extract=1 )
  endif else begin
    var_label_in = var_label
  endelse
endif
; Build the list of possible labels for the longitude dimension and variable 
; used in the file
if not( keyword_set( lon_label_in ) ) then begin
  if file_netcdf_read_geo_varinfo ne '' then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers='label', select_values='label=lon', $
        settings=lon_label_in
    lon_label_in = strsplit( lon_label_in, ',', extract=1 )
  endif else begin
    lon_label_in = 'lon'
  endelse
endif
; Build the list of possible labels for the latitude dimension and variable 
; used in the file
if not( keyword_set( lat_label_in ) ) then begin
  if file_netcdf_read_geo_varinfo ne '' then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers='label', select_values='label=lat', $
        settings=lat_label_in
    lat_label_in = strsplit( lat_label_in, ',', extract=1 )
  endif else begin
    lat_label_in = 'lat'
  endelse
endif
; Build the list of possible labels for the vertical dimension and variable 
; used in the file
if not( keyword_set( height_label_in ) ) then begin
  if file_netcdf_read_geo_varinfo ne '' then begin
    temp_label = [ 'height', 'plev' ]
    for i_label = 0, n_elements( temp_label ) - 1 do begin
      markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
          select_headers='label', select_values='label='+temp_label[i_label], $
          settings=temp_height_label_in
      temp_height_label_in = strsplit( temp_height_label_in, ',', extract=1 )
      if i_label eq 0 then begin
        height_label_in = temporary( temp_height_label_in )
      endif else begin
        height_label_in = [ height_label_in, temporary( temp_height_label_in ) ]
      endelse
    endfor
  endif else begin
    height_label_in = 'height'
  endelse
endif
; Build the list of possible labels for the time dimension and variable used 
; in the file
if not( keyword_set( time_label_in ) ) then begin
  if file_netcdf_read_geo_varinfo ne '' then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers='label', select_values='label=time', $
        settings=time_label_in
    time_label_in = strsplit( time_label_in, ',', extract=1 )
  endif else begin
    time_label_in = 'time'
  endelse
endif

;***********************************************************************
; Load Data

; Load longitude data
if keyword_set( lon_label_in ) then begin
  lon_data = netcdf_read( file_name, lon_label_in, units=lon_units, $
      attribute_label=lon_attribute_label, $
      attribute_value=lon_attribute_value, attribute_type=lon_attribute_type )
endif
n_lon = n_elements( lon_data )
if n_lon eq 0 then begin
  lon_data = nan
  n_lon = 1
endif else if n_lon eq 1 then begin
  if finite( lon_data[0] ) eq 0 then lon_label_in = ''
endif
; Load latitude data
if keyword_set( lat_label_in ) then begin
  lat_data = netcdf_read( file_name, lat_label_in, units=lat_units, $
      attribute_label=lat_attribute_label, $
      attribute_value=lat_attribute_value, attribute_type=lat_attribute_type )
endif
n_lat = n_elements( lat_data )
if n_lat eq 0 then begin
  lat_data = nan
  n_lat = 1
endif else if n_lat eq 1 then begin
  if finite( lat_data[0] ) eq 0 then lat_label_in = ''
endif
; Load height (vertical) data
if keyword_set( height_label_in ) then begin
  height_data = netcdf_read( file_name, height_label_in, units=height_units, $
      attribute_label=height_attribute_label, $
      attribute_value=height_attribute_value, $
      attribute_type=height_attribute_type )
endif
n_height = n_elements( height_data )
if n_height eq 0 then begin
  height_data = nan
  n_height = 1
endif else if n_height eq 1 then begin
  if finite( height_data[0] ) eq 0 then height_label_in = ''
endif
; Load time data
if keyword_set( time_label_in ) then begin
  time_attribute_label_all = 1
  time_data = netcdf_read( file_name, time_label_in, $
      attribute_label=time_attribute_label_all, $
      attribute_value=time_attribute_value_all, $
      attribute_type=time_attribute_type_all, units=temp_units )
  if not( keyword_set( time_units ) ) then begin
    if keyword_set( temp_units ) then time_units = temp_units
  endif
  if keyword_set( time_attribute_label ) then begin
    if str( time_attribute_label[0] ) eq '1' then begin
      time_attribute_label = time_attribute_label_all
      time_attribute_value = time_attribute_value_all
      time_attribute_type = time_attribute_type_all
    endif else begin
      n_time_attribute = n_elements( time_attribute_label )
      time_attribute_value = strarr( n_time_attribute )
      time_attribute_type = 7 + intarr( n_time_attribute )
      for i_attribute = 0, n_time_attribute - 1 do begin
        id = where( $
            time_attribute_label_all eq time_attribute_value[i_attribute], $
            n_id )
        if n_id eq 0 then begin
          time_attribute_value[i_attribute] = ''
        endif else begin
          time_attribute_value[i_attribute] = time_attribute_value_all[id[0]]
          time_attribute_type[i_attribute] = time_attribute_type_all[id[0]]
        endelse
      endfor
    endelse
  endif
endif
n_time = n_elements( time_data )
if n_time eq 0 then begin
  time_data = nan
  n_time = 1
endif else if n_time eq 1 then begin
  if finite( time_data[0] ) eq 0 then time_label_in = ''
endif
; Shift time if requested
if keyword_set( shift_time ) then time = time + shift_time

; Load realization data
if keyword_set( realization_label_in ) then begin
  temp = strsplit( realization_label_in, ',', extract=1 )
  realization_data = netcdf_read( file_name, temp, units=realization_units, $
      attribute_label=realization_attribute_label, $
      attribute_value=realization_attribute_value, $
      attribute_type=realization_attribute_type )
  ; In case it seems that we have two labels, one for the dimension and one for 
  ; the dimension variable
  if ( n_elements( realization_data ) eq 1 ) $
      and ( strpos( realization_label_in, ',' ) gt 0 ) $
      and ( finite( realization_data[0] ) eq 0 ) then begin
    temp_realization_label_in = strsplit( realization_label_in, ',', $
        extract=1, count=temp )
    if temp ne 2 then stop
    realization_data = netcdf_read( file_name, temp_realization_label_in[0], $
        units=realization_units )
    if ( n_elements( realization_data ) eq 1 ) $
        and ( finite( realization_data[0] ) eq 0 ) then begin
      realization_data = netcdf_read( file_name, temp_realization_label_in[1], $
          units=realization_units, $
          attribute_label=realization_attribute_label, $
          attribute_value=realization_attribute_value, $
          attribute_type=realization_attribute_type )
    endif
  endif
  if var_type( realization_data ) eq 1 then begin
    n_realization = n_elements( realization_data[0,*] )
    temp_realization_data = strarr( n_realization )
    for i_realization = 0, n_realization - 1 do begin
      temp_realization_data[i_realization] $
          = string( realization_data[*,i_realization] )
    endfor
    realization_data = temporary( temp_realization_data )
  endif
endif
n_realization = n_elements( realization_data )
if n_realization eq 0 then begin
  realization_data = nan
  n_realization = 1
endif

; Load the main data field
if keyword_set( var_label_in ) then begin
  ; Read data from file
  var_data = netcdf_read( file_name, var_label_in, $
      missing_label=missing_label_in, units=var_units, $
      global_label=global_label, global_value=global_value, $
      global_type=global_type, attribute_label=var_attribute_label, $
      attribute_value=var_attribute_value, attribute_type=var_attribute_type )
  ; Determine which dimensions this variable uses, as it may not use all the 
  ; dimensions listed in the file
  spawn, spawn_str + 'ncdump -h ' + file_name + ' | grep " ' + var_label_in $
      + '("', $
      dim_list
  ;if strpos( dim_list[0], 'netcdf ' ) eq 0 then begin
  ;  dim_list = dim_list[1]
  ;endif else begin
  ;  dim_list = dim_list[0]
  ;endelse
  if n_elements( dim_list ) gt 1 then begin
    if quiet_opt eq 0 then begin
      print, 'Warning netcdf_read_geo.pro:  Multiple instances of ' $
          + var_label_in + ' within input file ' + file_name $
          + '.  Assuming all but first are in attributes only.'
    endif
    dim_list = dim_list[0]
  endif else if n_elements( dim_list ) ne 1 then begin
    stop, 'Error netcdf_read_geo.pro:  Cannot isolate requested variable ' $
        + var_label_in + ' within input file ' + file_name + '.'
  endif
  dim_list = dim_list[0]
  dim_list = strsplit( dim_list, '(,)', extract=1, count=n_dim_list )
  dim_list = strtrim( dim_list, 2 )
  n_dim_list = n_dim_list - 2
  if n_dim_list lt 1 then stop
  dim_list = reverse( dim_list[1:n_dim_list] )
  n_dim = intarr( 5 )
  ctr_dim = 0
  if not( keyword_set( lon_label_in ) ) then begin
    n_dim[0] = 1
  endif else if dim_list[ctr_dim] eq lon_label_in then begin
    n_dim[0] = n_lon
    ctr_dim = ctr_dim + 1
  endif else begin
    n_dim[0] = 1
  endelse
  if not( keyword_set( lat_label_in ) ) then begin
    n_dim[1] = 1
  endif else if dim_list[ctr_dim] eq lat_label_in then begin
    n_dim[1] = n_lat
    ctr_dim = ctr_dim + 1
  endif else begin
    n_dim[1] = 1
  endelse
  if not( keyword_set( height_label_in ) ) then begin
    n_dim[2] = 1
  endif else if dim_list[ctr_dim] eq height_label_in then begin
    n_dim[2] = n_height
    ctr_dim = ctr_dim + 1
  endif else begin
    n_dim[2] = 1
  endelse
  if not( keyword_set( time_label_in ) ) then begin
    n_dim[3] = 1
  endif else if dim_list[ctr_dim] eq time_label_in then begin
    n_dim[3] = n_time
    ctr_dim = ctr_dim + 1
  endif else begin
    n_dim[3] = 1
  endelse
  if not( keyword_set( realization_label_in ) ) then begin
    n_dim[4] = 1
  endif else if max( dim_list[ctr_dim] $
      eq strsplit( realization_label_in, ',', extract=1 ) ) eq 1 then begin
    n_dim[4] = n_realization
    ctr_dim = ctr_dim + 1
  endif else begin
    n_dim[4] = 1
  endelse
  if n_elements( var_data ) ne long64( product( long( n_dim ) ) ) then stop
  ; Reform data array to lon*lat*height*time*realization format for convenience
  var_data = reform( var_data, n_dim )
; Or just read global attributes
endif else if keyword_set( global_label ) then begin
  ; Read data from file
  var_data = netcdf_read( file_name, '', global_label=global_label, $
      global_value=global_value, global_type=global_type )
endif

;***********************************************************************
; Adjust data units to standard format

; Set units check flag
check_units = 0
; If units defined
if keyword_set( var_units ) and not( keyword_set( no_units_conversion_opt ) ) $
    then begin
  ; Determine the standard units for this data variable
  if file_netcdf_read_geo_varinfo ne '' then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers='units', select_values='label='+var_label, $
        settings=var_units_standard
  endif else begin
    var_units_standard = ''
  endelse
  ; Convert data variable to standard units if possible
  if var_units_standard eq '' then begin
    print, 'Warning netcdf_read_geo.pro:  No standard units for ' + var_label $
        + ' variable.  Assuming okay.'
  endif else begin
    ; Standardise description of data units
    if max( strlowcase( var_units ) eq [ 'kelvin', 'deg_k' ] ) eq 1 then begin
      var_units = 'K'
    endif else if max( strlowcase( var_units ) eq [ 'kg/m^2', 'kg m**-2' ] ) $
        eq 1 then begin
      var_units = 'kg m-2'
    endif else if max( strlowcase( var_units ) $
        eq [ 'kg/m^2/s', 'kg m**-2 s*-1' ] ) eq 1 then begin
      var_units = 'kg m-2 s-1'
    endif else if max( var_units eq [ 'm/s' ] ) eq 1 then begin
      var_units = 'm s-1'
    endif
    ; Check if we need to convert units
    if var_units ne var_units_standard then begin
      ; Convert deg_C to K
      if ( max( strlowcase( var_units ) $
          eq [ 'celsius', 'degrees c', 'deg_c', 'degc' ] ) eq 1 ) $
          and ( var_units_standard eq 'K' ) then begin
        var_data = var_data + 273.15
        var_units = var_units_standard
      endif
      ; Convert hPa to Pa
      if ( var_units eq 'hPa' ) and ( var_units_standard eq 'Pa' ) then begin
        var_data = var_data * 100.
        var_units = var_units_standard
      endif
      ; Convert 1 to fraction
      if ( var_units eq '1' ) and ( var_units_standard eq 'fraction' ) $
          then begin
        var_units = var_units_standard
      endif
      ; Convert fraction to %
      if ( max( var_units eq [ 'fraction', '1' ] ) eq 1 ) $
          and ( var_units_standard eq '%' ) then begin
        var_data = var_data * 100.
        var_units = var_units_standard
      endif
      ; Convert % to fraction
      if ( var_units eq '%' ) and ( var_units_standard eq 'fraction' ) $
          then begin
        var_data = var_data / 100.
        var_units = var_units_standard
      endif
      ; Check if we successfully converted
      if var_units ne var_units_standard then begin
        stop, 'ERROR netcdf_read_geo.pro:  Unable to convert data variable ' $
            + 'to standard units.'
      endif
    endif
  endelse
endif

;***********************************************************************
; Ensure standard format of coordinates

; Ensure increasing latitude
if n_lat gt 1 then begin
  if min( lat_data[1:n_lat-1] - lat_data[0:n_lat-2] ) lt 0. then begin
    id_lat = sort( lat_data )
    lat_data = lat_data[id_lat]
    if keyword_set( var_data ) then var_data = var_data[*,id_lat,*,*,*]
    id_lat = 0
  endif
endif  

; Ensure increasing levels
if n_height gt 1 then begin
  if min( height_data[1:n_height-1] - height_data[0:n_height-2] ) lt 0. $
      then begin
    id_height = sort( height_data )
    height_data = height_data[id_height]
    if keyword_set( var_data ) then var_data = var_data[*,*,id_height,*,*]
    id_height = 0
  endif
endif  

; Determine the calendar type
if ( max( finite( time_data ) ) eq 1 ) and not( keyword_set( time_calendar ) ) $
    then begin
  if keyword_set( time_attribute_value_all ) then begin
    id = where( strpos( time_attribute_label_all, 'calendar' ) ge 0, n_id )
    if n_id gt 1 then begin
      id = where( time_attribute_label_all eq 'calendar', n_id )
      if n_id ne 1 then stop
    endif
    if time_attribute_value_all[id[0]] eq '' then stop
    time_calendar = time_attribute_value_all[id[0]]
  endif
  time_calendar = strlowcase( time_calendar )
  ; Assume "standard" calendar is Gregorian
  if time_calendar eq 'standard' then begin
    if quiet_opt eq 0 then begin
      print, 'Warning netcdf_read_geo.pro:  ' $
          + 'Assuming "standard" calendar is Gregorian.'
    endif
    time_calendar = 'gregorian'
  endif
endif

; Ensure increasing time vector
if n_time gt 1 then begin
  temp = time_data[1:n_time-1] - time_data[0:n_time-2]
  if min( temp ) lt 0 then stop
endif

; Convert to "days since" format
if not( keyword_set( keep_time_format_opt ) ) then begin
  if max( finite( time_data ) ) eq 1 then begin
    ; If a time origin has not been specified in the input
    if not( keyword_set( time_origin ) ) then begin
      ; Determine if there is an origin specified in the file's time:units 
      ; attribute
      if strpos( time_units, ' since ' ) gt 0 then begin
        temp = strsplit( time_units, ' ', extract=1, count=n_temp )
        if n_temp gt 4 then stop
        if n_temp lt 3 then stop
        if n_temp eq 4 then temp[2] = temp[2] + ':' + temp[3]
        time_origin = temp[2]
      endif
    endif
    if not( keyword_set( time_origin ) ) then stop
    ; Convert format if necessary
    if time_units ne 'days since ' + time_origin then begin
      time_data = convert_time_format( time_data, time_units, $
          'days since '+time_origin, calendar=time_calendar )
    endif
  endif
endif


;***********************************************************************
; Perform regional extraction

; Select region if requested
if keyword_set( region ) then begin
  ; Extract data
  if n_elements( var_data ) gt 0 then begin
    var_data = reform( var_data, n_lon, n_lat, $
        n_height * n_time * n_realization )
    var_data = extract_region( var_data, lat=lat_data, lon=lon_data, $
        region=region, noshrink=noshrink_opt )
    n_lon = n_elements( lon_data )
    n_lat = n_elements( lat_data )
    var_data = reform( var_data, n_lon, n_lat, n_height, n_time, n_realization )
  endif
endif

;***********************************************************************
; Ensure output

; Create default no-data output
if n_elements( var_data ) eq 0 then var_data = nan

; Release quiet option
if quiet_opt eq 1 then !except = 1

;***********************************************************************
; The End

;stop
return, var_data
END
