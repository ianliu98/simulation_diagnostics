;+
; NAME:
;    c20c_dtos_v2_make_delta
;
; PURPOSE:
;    This procedure produces the deltaSST fields for use in the C20C+ 
;    Detection and Attribution Project experiments and the Weather Risk 
;    Attribution Forecast.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_make_delta, attrib_formula=attrib_formula, $
;        indep_file=indep_file, indep_period=indep_period, $
;        indep_var_label=indep_var_label, out_file=out_file
;
; INPUTS:
;    ATTRIB_FORMULA, GRID_FILE, INDEP_FILE, INDEP_PERIOD, INDEP_VAR_LABEL, 
;      OPTIONS, OUT_FILE, TIME_FILTER, OUT_FILE, OUT_GLOBAL_ATTRIBUTE_C20C, 
;      OUT_GLOBAL_ATTRIBUTE_LABEL, OUT_GLOBAL_ATTRIBUTE_TYPE, 
;      OUT_GLOBAL_ATTRIBUTE_VALUE, OUT_PERIOD, OUT_VAR_ATTRIBUTE_LABEL, 
;      OUT_VAR_ATTRIBUTE_TYPE, OUT_VAR_ATTRIBUTE_VALUE, OUT_VAR_LABEL
;
; KEYWORD PARAMETERS:
;    ATTRIB_FORMULA:  A required scalar string containing the formula to use in 
;      when using the indepependent scenarios to calculate the delta field.  
;      For two scenario, an example might be '+1-2', which would subtract the 
;      second scenario from the first scenario.  See algebra_to_matrix.pro for 
;      more information.
;    GRID_FILE:  An optional scalar string naming the file containing the 
;        target longitude-latitude grid  for interpolation of the output.
;    INDEP_FILE:  A required N_INDEP_REALIZATION,N_INDEP_SCENARIO string array 
;        containing the lists of files containing the data for the independent 
;        scenarios.  For each of the N_INDEP_REALIZATION realizations (e.g. 
;        simulations) of the N_INDEP_SCENARIO scenarios, the list of files 
;        should be comma-delimited and can include the "*" wildcard.  Of fewer 
;        than N_INDEP_REALIZATION realizations exist for a given scenario, the 
;        remaining entries should be left blank ('').
;    INDEP_PERIOD:  A required 2-element or 2,N_INDEP_SCENARIO string array 
;        describing the start and end dates of the period over which to 
;        estimate the attributable delta.  Of format [<start>,<end>] with the 
;        dates in the 'yyyymm' format.  If of size 2*N_INDEP_SCENARIO, where 
;        N_INDEP_SCENARIO is the number of independent scenarios to be used in 
;        estimating the attributable difference, this can be used for instance 
;        if the estimate is based on the change between two time periods in the 
;        observed record.
;    INDEP_VAR_LABEL:  A required scalar string or vector string of length 
;        N_INDEP_SCENARIO listing the variable(s) to load from the 
;        N_INDEP_SCENARIO scenarios.  If a vector is input, the variable labels 
;        can differ across scenarios.
;    OPTIONS:  An optional vector string listing various processing options to 
;        implement.  Possible entries are:
;        * 'extend last year':  If data from an input source ends before the 
;          end of INDEP_PERIOD[1,i_scenario], then if this option is input then 
;          those extra years will be taken as a repeat of the annual cycle of 
;          the last year with data.  Note that this option is implemented after 
;          application of a temporal filter, so if a boxcar filter of length 3 
;          is used, the repeat of the annual cycle of the final three years 
;          with data will be used.
;    OUT_FILE:  A required scalar string containing the name of the output 
;        NetCDF file containing the delta pattern.  This can include a 
;        directory path.  The file name must follow the C20C+ D&A conventions.
;    OUT_GLOBAL_ATTRIBUTE_C20C:  If set then a list of default global 
;        attributes for the C20C+ D&A project will be added to the 
;        OUT_GLOBAL_ATTRIBUTE_* lists.
;    OUT_GLOBAL_ATTRIBUTE_LABEL:  An optional string vector of length 
;        N_OUT_GLOBAL_ATTRIBUTE containing a list of global attribute labels to 
;        be included in the output NetCDF file.  The entries here correspond to 
;        the entries in OUT_GLOBAL_ATTRIBUTE_TYPE and 
;        OUT_GLOBAL_ATTRIBUTE_VALUE.
;    OUT_GLOBAL_ATTRIBUTE_TYPE:  An optional integer vector of length 
;        N_OUT_GLOBAL_ATTRIBUTE containing a list of the variable type of the 
;        global attribute values to be included in the output NetCDF file.  The 
;        entries here correspond to the entries in OUT_GLOBAL_ATTRIBUTE_VALUE.  
;        The default is a string array with all entries being 7 (type string).
;    OUT_GLOBAL_ATTRIBUTE_VALUE:  An optional string vector of length 
;        N_OUT_GLOBAL_ATTRIBUTE containing a list of global attribute values to 
;        be included in the output NetCDF file.  The entries here correspond to 
;        the entries in OUT_GLOBAL_ATTRIBUTE_LABEL and 
;        OUT_GLOBAL_ATTRIBUTE_TYPE.
;    OUT_PERIOD:  An optional two-element string vector containing the start 
;        and end dates respectively of the time period of the output data.  
;        Dates should be of the 'yyyymm' format.  The default is 
;        INDEP_PERIOD[*,0].
;    OUT_VAR_ATTRIBUTE_LABEL:  An optional string vector of length 
;        N_OUT_VAR_ATTRIBUTE containing a list of attribute labels for the data 
;        variable to be included in the output NetCDF file.  The entries here 
;        correspond to the entries in OUT_VAR_ATTRIBUTE_TYPE and 
;        OUT_VAR_ATTRIBUTE_VALUE.
;    OUT_VAR_ATTRIBUTE_TYPE:  An optional integer vector of length 
;        N_OUT_VAR_ATTRIBUTE containing a list of the variable type of the 
;        attribute values for the data variable to be included in the output 
;        NetCDF file.  The entries here correspond to the entries in 
;        OUT_VAR_ATTRIBUTE_VALUE.  The default is a string array with all 
;        entries being 7 (type string).
;    OUT_VAR_ATTRIBUTE_VALUE:  An optional string vector of length 
;        N_OUT_VAR_ATTRIBUTE containing a list of attribute values for the data 
;        variable to be included in the output NetCDF file.  The entries here 
;        correspond to the entries in OUT_VAR_ATTRIBUTE_LABEL and 
;        OUT_VAR_ATTRIBUTE_TYPE.
;    OUT_VAR_LABEL:  An optional scalar string specifying the label of the 
;        variable in the output data.  The default is INDEP_VAR_LABEL[0].
;    TIME_FILTER:  An optional two-element string vector with the first element 
;        describing the type of temporal filter to use and the second element 
;        describing the length of the filter.  For example ['boxcar','5'] would 
;        apply a boxcar filter of length 5.
;
; OUTPUTS:
;    The file specified in OUT_FILE.
;
; USES:
;    add_dim.pro
;    algebra_to_matrix.pro
;    convert_time_format.pro
;    filter.pro
;    month_name.pro
;    netcdf_read_geo.pro
;    netcdf_read_geo_multitime.pro
;    process_lonlatmonth.pro
;    str.pro
;    string_from_vector.pro
;    string_substitute.pro
;
; PROCEDURE:
;    This procedure loads data from various sources, performs the specified 
;    arithmetic on that data to produce a delta field, and saves the delta 
;    to a NetCDF file.
;
; EXAMPLES:
;    See c20c_dtos_v2_make_delta_driver.pro
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2011-01-05, as 
;        c20c_dsst_generate.pro.
;    Modified:  DAS, 2013-06-10
;    Modified:  DAS, 2016-05-05
;    Modified:  DAS, 2018-10-08 (Branched from c20c_dsst_generate.pro;  Switched from use of c20c_dsst_scenario_arithmetic.pro to algebra_to_matrix.pro)
;    Modified:  DAS, 2018-11-14 (Fixed bug in file name of IDL restore dump)
;    Modified:  DAS, 2018-12-12 (Completed documentation.  Added 
;        OUT_GLOBAL_ATTRIBUTE_C20C option.)
;-

PRO C20C_DTOS_V2_MAKE_DELTA, $
    ATTRIB_FORMULA=attrib_formula, $
    GRID_FILE=grid_file, $
    INDEP_FILE=indep_file, INDEP_PERIOD=indep_period, $
      INDEP_VAR_LABEL=indep_var_label, $
    OPTIONS=options, $
    OUT_FILE=out_file, OUT_PERIOD=out_period, OUT_VAR_LABEL=out_var_label, $
    OUT_GLOBAL_ATTRIBUTE_LABEL=out_global_attribute_label, $
      OUT_GLOBAL_ATTRIBUTE_TYPE=out_global_attribute_type, $
      OUT_GLOBAL_ATTRIBUTE_VALUE=out_global_attribute_value, $
      OUT_GLOBAL_ATTRIBUTE_C20C=out_global_attribute_c20c_opt, $
    OUT_VAR_ATTRIBUTE_LABEL=out_var_attribute_label, $
      OUT_VAR_ATTRIBUTE_TYPE=out_var_attribute_type, $
      OUT_VAR_ATTRIBUTE_VALUE=out_var_attribute_value, $
    TIME_FILTER=time_filter

;***********************************************************************
; Constants

; Not-a-number
nan = !values.f_nan
; Months in a year
mina = 12
; Hard carriage return character
hard_return = string( 10B )

;***********************************************************************
; Check for required input

; Required input data for calculating the attributable change.
; The list of input files
if not( keyword_set( indep_file ) ) then stop
n_indep_scenario = n_elements( indep_file[0,*] )
n_indep_realization = n_elements( indep_file[*,0] )
; The time periods to use
if not( keyword_set( indep_period ) ) then stop
if n_elements( indep_period[*,0] ) ne 2 then stop
; Ensure that the period array is consistent
if n_elements( indep_period[0,*] ) ne n_indep_scenario then begin
  if n_elements( indep_period[0,*] ) ne 1 then stop
  indep_period = add_dim( indep_period, 1, n_indep_scenario )
endif
; The labels of the variables to load
if not( keyword_set( indep_var_label ) ) then stop
if n_elements( indep_var_label ) ne n_indep_scenario then begin
  if n_elements( indep_var_label ) ne 1 then stop
  indep_var_label = indep_var_label[0] + strarr( n_indep_scenario )
endif

; The number of attributable estimates to calculate
n_attrib_formula = n_elements( attrib_formula )
if n_attrib_formula ne 1 then stop

; The file containing the target grid
if n_elements( grid_file ) gt 1 then stop

; The output file
if n_elements( out_file ) ne 1 then stop
; The output variable label
if not( keyword_set( out_var_label ) ) then out_var_label = indep_var_label[0]

; Default options (none)
if not( keyword_set( options ) ) then options = ''

; Ensure legal time filter request
if keyword_set( time_filter ) then begin
  if n_elements( time_filter ) ne 2 then stop
  if fix( time_filter[1] ) eq 0 then stop
  ; Stretch the period to load in order to apply filter on ends
  indep_period_load = [ $
      str( long( strmid( indep_period[0,*], 0, 4 ) ) $
      - ( long( time_filter[1] ) - 1 ) / 2, length=4, filler='0' ) $
      + strmid( indep_period[0,*], 4, 2 ), $
      str( long( strmid( indep_period[1,*], 0, 4 ) ) $
      + long( time_filter[1] ) / 2, length=4, filler='0' ) $
      + strmid( indep_period[1,*], 4, 2 ) ]
; Default dummy filter
endif else begin
  time_filter = [ '', '' ]
  indep_period_load = indep_period
endelse

; Default null global attribute in the output file
if not( keyword_set( out_global_attribute_label ) ) then begin
  out_global_attribute_label = ''
  out_global_attribute_value = ''
  out_global_attribute_type = 7
; Otherwise check global attribute inputs and adjust format
endif else begin
  n_out_global_attribute = n_elements( out_global_attribute_label )
  if n_elements( out_global_attribute_value ) ne n_out_global_attribute $
      then stop
  if n_elements( out_global_attribute_type ) eq 0 then begin
    out_global_attribute_type = 7 + intarr( n_out_global_attribute )
  endif else begin
    if n_elements( out_global_attribute_type ) ne n_out_global_attribute $
        then stop
  endelse
  for i_global = 0, n_out_global_attribute - 1 do begin
    if out_global_attribute_type[i_global] eq 7 then begin
      out_global_attribute_value[i_global] = string_substitute( $
          out_global_attribute_value[i_global], '\\', hard_return, regex=1 )
    endif
  endfor
endelse

;***********************************************************************
; Load data

; Load target grid dimension
if keyword_set( grid_file ) then begin
  temp = netcdf_read_geo( grid_file, lon=grid_lon, lat=grid_lat )
endif

; Initialise time metadata variables
indep_time_units = strarr( n_indep_scenario )
indep_time_calendar = strarr( n_indep_scenario )

; Load independent variable data.
; Iterate through independent scenarios
for i_indep = 0, n_indep_scenario - 1 do begin
  ; Iterate through realisations
  for i_realization = 0, n_indep_realization - 1 do begin
    ; Print status
    print, 'Loading independent realization files ' $
        + indep_file[i_realization,i_indep]
    ; If files are defined for this realization
    if indep_file[i_realization,i_indep] ne '' then begin
      ; Load data for this realization
      temp_indep_time_units = ''
      temp_indep_time_calendar = ''
      temp_file = strsplit( indep_file[i_realization,i_indep], ',', extract=1 )
      temp_data = netcdf_read_geo_multitime( temp_file, $
          indep_var_label[i_indep], mask_lon=grid_lon, mask_lat=grid_lat, $
          period_time=indep_period_load[*,i_indep], lon=indep_lon, $
          lat=indep_lat, time=temp_indep_time, $
          units_time=temp_indep_time_units, calendar=temp_indep_time_calendar, $
          quiet=1 )
      n_temp_indep_time = n_elements( temp_indep_time )
      ; Record dimensions and metadata
      if indep_time_units[i_indep] eq '' then begin
        ; Copy time metadata
        indep_time_units[i_indep] = temp_indep_time_units
        indep_time_calendar[i_indep] = temp_indep_time_calendar
        ; Extract time data variable
        if time_filter[1] ne '' then begin
          temp_time = convert_time_format( temp_indep_time, $
              temp_indep_time_units, 'yyyymm', $
              calendar=temp_indep_time_calendar )
          id = where( ( temp_time ge indep_period[0,i_indep] ) $
              and ( temp_time le indep_period[1,i_indep] ), temp_n_indep_time )
          if keyword_set( n_indep_time ) then begin
            if temp_n_indep_time ne n_indep_time then stop
          endif else begin
            n_indep_time = temp_n_indep_time
          endelse
          temp = n_temp_indep_time - ( long( time_filter[1] ) - 1 ) * mina
          if n_indep_time ne temp then stop
          if i_indep eq 0 then begin
            indep_time = nan + fltarr( n_indep_time, n_indep_scenario )
          endif
          indep_time[*,i_indep] = temp_indep_time[id]
        endif else begin
          if i_indep eq 0 then begin
            n_indep_time = n_elements( temp_indep_time )
            indep_time = nan + fltarr( n_indep_time, n_indep_scenario )
          endif
          indep_time[*,i_indep] = temp_indep_time
        endelse
        ; Record spatial dimension sizes
        if not( keyword_set( n_indep_lon ) ) then begin
          n_indep_lon = n_elements( indep_lon )
          n_indep_lat = n_elements( indep_lat )
        endif
      endif
      ; Reform data array to standard format
      temp_data = reform( temp_data, n_indep_lon, n_indep_lat, $
          n_temp_indep_time )
      ; Apply temporal smoothing
      if time_filter[0] ne '' then begin
        for i_month = 0, mina - 1 do begin
          index_month = ( n_temp_indep_time + mina - i_month - 1 ) / mina
          index_month = i_month + mina * indgen( index_month )
          temp_data_month = temp_data[*,*,index_month]
          for i_lon = 0, n_indep_lon - 1 do begin
            for i_lat = 0, n_indep_lat - 1 do begin
              temp_data_month[i_lon,i_lat,*] = filter( $
                  reform( temp_data_month[i_lon,i_lat,*] ), $
                  fix( time_filter[1] ), time_filter[0], nan=0 )
            endfor
          endfor
          temp_data[*,*,index_month] = temp_data_month
        endfor
        temp_data_month = 0
        ; Restrict to indep_period 
        temp_time = convert_time_format( temp_indep_time, $
            temp_indep_time_units, 'yyyymm', calendar=temp_indep_time_calendar )
        id = where( ( temp_time ge indep_period[0,i_indep] ) $
            and ( temp_time le indep_period[1,i_indep] ), n_id )
        if n_id ne n_indep_time then stop
        temp_data = temp_data[*,*,id]
      endif
      ; If we are extending the last year where needed (because of missing 
      ; values at the end)
      if max( options eq 'extend last year' ) eq 1 then begin
        temp_finite = total( total( finite( temp_data ), 1 ), 1 )
        if temp_finite[n_indep_time-1] lt 0.5 then begin
          id_good_last = max( where( temp_finite gt 0.5 ) )
          id_good = id_good_last - mina + 1 + indgen( mina )
          for i_month = 0, mina - 1 do begin
            n_id_miss = ( n_indep_time - id_good[i_month] - 1 ) / mina
            if n_id_miss gt 0 then begin
              id_miss = id_good[i_month] + mina + mina * indgen( n_id_miss )
              for i_miss = 0, n_id_miss - 1 do begin
                temp_data[*,*,id_miss[i_miss]] = temp_data[*,*,id_good[i_month]]
              endfor
            endif
          endfor
        endif
      endif
      ; Record data
      if not( keyword_set( indep_data ) ) then begin
        indep_data = fltarr( n_indep_lon, n_indep_lat, n_indep_time, $
            n_indep_scenario )
        indep_realization_ctr = intarr( n_indep_lon, n_indep_lat, $
            n_indep_time, n_indep_scenario )
      endif
      temp_data_finite = finite( temp_data )
      id = where( temp_data_finite eq 0, n_id )
      if n_id gt 0 then temp_data[id] = 0.
      indep_data[*,*,*,i_indep] = indep_data[*,*,*,i_indep] + temp_data
      indep_realization_ctr[*,*,*,i_indep] $
          = indep_realization_ctr[*,*,*,i_indep] + temp_data_finite
      ; Add Arctic- and domain-mean time series, for quality control
      if not( keyword_set( indep_series_global ) ) then begin
        indep_series_global = fltarr( n_indep_time, n_indep_scenario )
        indep_series_global_finite = intarr( n_indep_time, n_indep_scenario )
        indep_series_arctic = fltarr( n_indep_time, n_indep_scenario )
        indep_series_arctic_finite = intarr( n_indep_time, n_indep_scenario )
      endif
      temp = temp_data
      temp_lon = indep_lon
      temp_lat = indep_lat
      process_lonlatmonth, temp, lon=temp_lon, lat=temp_lat, $
          integrate='mean=1,2'
      temp = reform( temp )
      temp_finite = finite( temp )
      id = where( temp_finite eq 0, n_id )
      if n_id gt 0 then temp[id] = 0.
      indep_series_global[*,i_indep] = indep_series_global[*,i_indep] + temp
      indep_series_global_finite[*,i_indep] $
          = indep_series_global_finite[*,i_indep] + temp_finite
      id_lat = where( indep_lat gt 60., n_id_lat )
      if n_id_lat gt 0 then begin
        temp = temp_data[*,id_lat,*]
        temp_lon = indep_lon
        temp_lat = indep_lat[id_lat]
        process_lonlatmonth, temp, lon=temp_lon, lat=temp_lat, $
            integrate='mean=1,2'
        temp = reform( temp )
        temp_finite = finite( temp )
        id = where( temp_finite eq 0, n_id )
        if n_id gt 0 then temp[id] = 0.
        indep_series_arctic[*,i_indep] = indep_series_arctic[*,i_indep] + temp
        indep_series_arctic_finite[*,i_indep] $
            = indep_series_arctic_finite[*,i_indep] + temp_finite
      endif
      ; Clear memory
      temp_data = 0
      temp_data_finite = 0
    endif
  endfor
endfor
; Calculate averages across realisations
indep_data = indep_data / indep_realization_ctr
id = where( indep_realization_ctr eq 0, n_id )
if n_id gt 0 then indep_data[id] = nan

;***********************************************************************
; Calculate attributable change

; Generate the matrix algebra for calculating the attributable change
attrib_transform = algebra_to_matrix( attrib_formula, n_var=n_indep_scenario )
; Calculate the attributable change
attrib_data = reform( indep_data, n_indep_lon * n_indep_lat * n_indep_time, $
    n_indep_scenario ) # attrib_transform
attrib_data = reform( attrib_data, n_indep_lon, n_indep_lat, n_indep_time )
; Copy dimension
attrib_lon = indep_lon
attrib_lat = indep_lat
; Generate attrib_time if required
if keyword_set( out_period ) then begin
  ; Define time vector properties
  n_attrib_time = n_indep_time
  attrib_time_units =  'days since ' + strmid( out_period[0], 0, 4 ) + '-' $
      + strmid( out_period[0], 4, 2 ) + '-01-00-00-00'
  attrib_time_calendar = 'gregorian'
  ; Generate new time data from scratch
  attrib_time = long( strmid( out_period[0], 4, 2 ) ) - 1 $
      + lindgen( n_indep_time + 1 )
  attrib_time = str( long( strmid( out_period[0], 0, 4 ) ) $
      + attrib_time / mina, $
      length=4, filler='0' ) $
      + str( ( attrib_time mod mina ) + 1, length=2, filler='0' )
  if attrib_time[n_attrib_time-1] ne out_period[1] then stop
  attrib_time = convert_time_format( attrib_time + '01', 'yyyymmdd', $
        attrib_time_units, calendar=attrib_time_calendar )
  attrib_time = ( attrib_time[0:n_attrib_time-1] $
      + attrib_time[1:n_attrib_time] ) / 2.
; Otherwise just copy indep_time from the first scenario
endif else begin
  attrib_time = indep_time[*,0]
  attrib_time_units = indep_time_units[0]
  attrib_time_calendar = indep_time_calendar[0]
endelse

; Scale by regression coefficient
if max( options eq 'regress' ) eq 1 then begin
  stop, 'c20c_dtos_v2_make_delta.pro: "regress" option not yet implemented.'
  attrib_data = beta[0] * attrib_data
endif

;***********************************************************************
; Define attributes for output file

; Define time attributes
out_time_attribute_label = [ 'standard_name', 'long_name', 'axis', 'units', $
    'calendar' ]
out_time_attribute_value = [ '', '', '', attrib_time_units, $
    attrib_time_calendar ]
out_time_attribute_type = [ 7, 7, 7, 7, 7 ]

; If C20C+ D&A project global attributes are requested
if keyword_set( out_global_attribute_c20c_opt ) then begin
  ; Add project_id attribute
  if max( out_global_attribute_label eq 'project_id' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'project_id' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        'C20C+ Detection and Attribution Project' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add license attribute
  if max( out_global_attribute_label eq 'license' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'license' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        'Creative Commons License: ' $
        + 'http://creativecommons.org/licenses/by-nc-sa/2.0/' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add references attribute
  if max( out_global_attribute_label eq 'references' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'references' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        'http://portal.nersc.gov/c20c' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add parent_experiment_family attrbute
  if max( out_global_attribute_label eq 'parent_experiment_family' ) eq 0 $
      then begin
    out_global_attribute_label = [ out_global_attribute_label, $
        'parent_experiment_family' ]
    out_global_attribute_value = [ out_global_attribute_value, 'N/A' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add parent_experiment attrbute
  if max( out_global_attribute_label eq 'parent_experiment' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, $
        'parent_experiment' ]
    out_global_attribute_value = [ out_global_attribute_value, 'N/A' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add parent_subexperiment attrbute
  if max( out_global_attribute_label eq 'parent_subexperiment' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, $
        'parent_subexperiment' ]
    out_global_attribute_value = [ out_global_attribute_value, 'N/A' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Add parent_run_id attrbute
  if max( out_global_attribute_label eq 'parent_run_id' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'parent_run_id' ]
    out_global_attribute_value = [ out_global_attribute_value, 'N/A' ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Determine attributes based on file name
  out_file_parsed = strsplit( out_file, '/', extract=1, $
      count=n_out_file_parsed )
  out_file_parsed = strsplit( out_file_parsed[n_out_file_parsed-1], '_', $
      extract=1, count=n_out_file_parsed )
  if n_out_file_parsed ne 8 then stop
  if max( out_global_attribute_label eq 'realm' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'realm' ]
    if max( out_file_parsed[1] eq [ 'O3hr', 'Oday', 'Omon' ] ) eq 1 then begin
      out_global_attribute_value = [ out_global_attribute_value, 'ocean' ]
    endif else begin
      stop
    endelse
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'frequency' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'frequency' ]
    if max( out_file_parsed[1] eq [ 'Oday' ] ) eq 1 then begin
      out_global_attribute_value = [ out_global_attribute_value, 'day' ]
    endif else if max( out_file_parsed[1] eq [ 'Omon' ] ) eq 1 $
        then begin
      out_global_attribute_value = [ out_global_attribute_value, 'mon' ]
    endif else begin
      stop
    endelse
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'model_id' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'model_id' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        out_file_parsed[2] ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'experiment_family' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, $
        'experiment_family' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        out_file_parsed[3] ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'experiment' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'experiment' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        out_file_parsed[4] ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'subexperiment' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'subexperiment' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        out_file_parsed[5] ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'run_id' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'run_id' ]
    out_global_attribute_value = [ out_global_attribute_value, $
        out_file_parsed[6] ]
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  ; Warn if the following global attributes are missing
  temp_global = [ 'title', 'contact', 'institution', 'institute_id' ]
  id = where( isin( out_global_attribute_label, temp_global ) eq 0, n_id )
  if n_id gt 0 then begin
    print, 'WARNING (c20c_dtos_v2_make_delta.pro):  The following global ' $
        + 'attributes are missing but should be present:  ' $
        + string_from_vector( temp_global[id], addand=1 ) + '.'
  endif
endif

; Add creation date global attribute
creation_date = systime( utc=1 )
creation_date = strsplit( creation_date, ' ', extract=1, count=n_creation_date )
if n_creation_date ne 5 then stop
temp_year = creation_date[4]
temp_month = where( month_name( indgen( 12 ), abbreviate=1 ) $
    eq creation_date[1] )
temp_month = str( temp_month[0] + 1, length=2, filler='0' )
temp_day = creation_date[2]
temp_time = creation_date[3]
creation_date = temp_year + '-' + temp_month + '-' + temp_day + 'T' $
    + temp_time + 'Z'
out_global_attribute_label = [ out_global_attribute_label, 'creation_date' ]
out_global_attribute_value = [ out_global_attribute_value, creation_date ]
out_global_attribute_type = [ out_global_attribute_type, 7 ]
; Add or update history global attribute
temp_history = creation_date + ': Processed and written to ' + out_file $
    + ' by c20c_dtos_v2_make_tossic.pro.'
id = where( out_global_attribute_label eq 'history', n_id )
if n_id gt 1 then stop
if n_id eq 1 then begin
  if out_global_attribute_value[id[0]] eq '' then begin
    out_global_attribute_value[id[0]] = temp_history
  endif else begin
    out_global_attribute_value[id[0]] = out_global_attribute_value[id[0]] $
        + ' ' + hard_return + ' ' + temp_history
  endelse
endif else begin
  out_global_attribute_label = [ out_global_attribute_label, 'history' ]
  out_global_attribute_value = [ out_global_attribute_value, temp_history ]
endelse
out_global_attribute_type = [ out_global_attribute_type, 7 ]

;***********************************************************************
; Write the adjusted data to NetCDF file

; Write data to NetCDF file
netcdf_write, out_file, data_array=attrib_data, data_label=out_var_label, $
    data_attribute_label=out_var_attribute_label, $
    data_attribute_value=out_var_attribute_value, $
    data_attribute_type=out_var_attribute_type, $
    dim1_vector=attrib_lon, dim1_label='lon', $
    dim2_vector=attrib_lat, dim2_label='lat', $
    dim3_vector=attrib_time, dim3_label='time', $
    dim3_attribute_label=out_time_attribute_label, $
    dim3_attribute_value=out_time_attribute_value, $
    dim3_attribute_type=out_time_attribute_type, $
    global_attribute_label=out_global_attribute_label, $
    global_attribute_value=out_global_attribute_value, $
    global_attribute_type=out_global_attribute_type

; Print progress
print, 'Results of c20c_dtos_v2_make_delta.pro written to ' + out_file + '.'

;***********************************************************************
; The End

;stop
return
END
