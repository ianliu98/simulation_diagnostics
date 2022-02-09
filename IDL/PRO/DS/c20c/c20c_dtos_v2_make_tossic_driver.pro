;+
; NAME:
;    c20c_dtos_v2_make_tossic_driver
;
; PURPOSE:
;    This procedure runs c20c_dtos_v2_make_tossic.pro to generate standard 
;    versions of sea surface temperature and sea ice concentration data for use 
;    in the C20C+ Detection and Attribution project.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_make_tossic_driver, out_file
;
; INPUTS:
;    OUT_FILE:  A required scalar string containing the output NetCDF file 
;        name.  The format should be 
;        <variable>_<domain>_<model_id>_<experiment_family>_<experiment>_<subexperiment>_<realization>_<start-date>-<end-date>.nc.
;    CONTACT, DELTA_FILE, DELTA_PERIOD, GRID_FILE, IN_FILE, IN_PERIOD
;
; KEYWORD PARAMETERS:
;    CONTACT:  An optional scalar string containing the name and/or contact 
;        details of the person running the code.  If input, this gets added to 
;        any default value listed in c20c_dtos_v2_map_model.xml.
;    DELTA_EXTEND_CYCLE:  If set, then the delta data is assumed to cover only 
;        one annual cycle, and is extended to cover the full IN_PERIOD period.
;    DELTA_FILE:  An optional string vector containing the names of the files 
;        containing the sea surface temperature delta field to be 
;        added/subtracted to the sea surface temperature, if such a delta is 
;        being applied.  A default will be provided in 
;        c20c_dtos_v2_map_experiment.xml, if a delta is being applied, but the 
;        directory will be unknown unless the files are in the current 
;        directory or the NERSC option is set.
;    DELTA_PERIOD:  An optional 2*N_DELTA string array containing the time 
;        period for which to restrict data read from the DELTA_FILE files.  The 
;        format should be 'yyyymmdd', so for instance ['20010101','20101231'] 
;        would input data covering the 1 January 2001 through 31 December 2010 
;        period.  The default is IN_PERIOD.
;    GRID_FILE:  An optional string scalar containing the name of the file 
;        containing the grid dimensions for use in the output file.  A default 
;        will be provided in c20c_dtos_v2_map_experiment.xml, but the directory 
;        will be unknown unless the files are in the current directory or the 
;        NERSC option is set.
;    IN_FILE:  A string vector containing the list of files containing the data 
;        to prepare for C20C+ D&A use.  Typically the data will be observed sea 
;        surface temperature and sea ice concentration.  A file can contain 
;        both variables, and/or the variables can be randomly found in 
;        different files in the list.  The code will sort out which file(s) to 
;        use for which variable.  Entries can include the "*" wildcard, for 
;        instance '*.nc'. Required if the NERSC option is not set.
;    IN_PERIOD:  An optional two-element string vector defining the input 
;        start and end months of the period of time to load from the input 
;        data.  Each entry is of the format 'yyyymm'.  The default is the 
;        "<start-date>-<end-date>" component of OUT_FILE.
;    NERSC:  If set then input directories are assumed to be the standard C20C+ 
;        D&A archive directories on NERSC disk systems.  If not set then 
;        DELTA_FILE, GRID_FILE and IN_FILE must be input.
;    READY_TO_USE:  If set then the output file is modified such that it is 
;        ready for immediate use by the intended climate model.  The default is 
;        an output file that conforms to C20C+ D&A protocols.  For example, 
;        C20C+ D&A protocols dictate that the sea surface temperature variable 
;        should be labeled as "tos", but if the model is CAM5-1-1degree and 
;        this option is set, then the variable will be labeled as 
;        "SST_cpl_prediddle".  This option only works with certain climate 
;        models.
;
; OUTPUTS:
;    ---
;
; USES:
;    date
;    ncks
;    add_dim.pro
;    c20c_dtos_v2_download_obs.pro
;    c20c_dtos_v2_make_tossic.pro
;    c20c_dtos_v2_map_institute.xml
;    c20c_dtos_v2_map_model.xml
;    c20c_dtos_v2_map_scenario.xml
;    isin.pro
;    markup_read.pro
;    month_day.pro
;    netcdf_read_geo_multitime.pro
;    netcdf_read_metadata.pro
;    str.pro
;    string_from_vector.pro
;
; PROCEDURE:
;    This procedure uses information from a set of .xml library files in order 
;    to generate c20c_dtos_v2_make_tossic.pro calls following standard 
;    settings.  It also understands where data should be sitting on NERSC 
;    systems and so requires less guidance to find relevant files.
;
; EXAMPLES:
;    ; When on NERSC systems the following produces the ocean surface input 
;    ; file needed to run the LBNL/CAM5-1-1degree All-Hist/est1/v2-0 
;    ; simulations during 2015/01/01 through 2015/12/31.
;    out_file = 'tosbcs-sicbcs_Omon_CAM5-1-1degree_All-Hist_est1_v2-0_period201412-201601_201412-201601.nc'
;    c20c_dtos_v2_make_tossic_driver, out_file, nersc=1, ready_to_use=1
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-12-23.
;    Modified:  DAS, 2018-03-02 (Fixed error in default variable attribute 
;        settings;  added file search on NERSC systems;  modified interaction 
;        with .xml library files;  added CONTACT keyword input)
;    Modified:  DAS, 2018-05-06 (Added missing GRID_FILE keyword in call to 
;        c20c_dtos_v2_make_tossic.pro)
;    Modified:  DAS, 2018-08-23 (Completed header documentation;  Added missing 
;        out_time_label usage)
;    Modified:  DAS, 2018-11-16 (Added DELTA_PERIOD and IN_PERIOD keyword 
;        inputs, and DELTA_EXTEND_CYCLE keyord option;  Added more attribute 
;        defaults for the CAM series of models under READY_TO_USE=1;  Updated 
;        default input directories when NERSC=1;  Allowed usage of multiple 
;        delta patterns;  Added default DELTA_APPLY_VAR_LABEL in call to 
;        c20c_dtos_v2_make_tossic.pro)
;    Modified:  DAS, 2018-11-29 (Fixed bugs in new observations download 
;        section)
;-

PRO C20C_DTOS_V2_MAKE_TOSSIC_DRIVER, $
    OUT_FILE, $
    IN_FILE=in_file, IN_PERIOD=in_period, $
    CONTACT=attribute_contact, $
    DELTA_FILE=delta_file, DELTA_PERIOD=delta_period, $
      DELTA_EXTEND_CYCLE=delta_extend_cycle_opt, $
    GRID_FILE=grid_file, $
    NERSC=nersc_opt, $
    READY_TO_USE=ready_to_use_opt

;***********************************************************************
; Constants and options

; Set to run in "v1-0" mode
v1_opt = 1

; Determine the directory of the main code (and associated library files)
dir_c20c_dtos_v2 = file_which( 'c20c_dtos_v2_make_tossic.pro' )
if not( keyword_set( dir_c20c_dtos_v2 ) ) then stop
temp = strsplit( dir_c20c_dtos_v2, '/', extract=1, count=n_temp )
if n_temp eq 1 then begin
  dir_c20c_dtos_v2 = ''
endif else begin
  if strmid( dir_c20c_dtos_v2, 0, 1 ) eq '/' then temp[0] = '/' + temp[0]
  dir_c20c_dtos_v2 = string_from_vector( temp[0:n_temp-2], spacer='/', $
      nospace=1 ) + '/'
endelse
; Define experiment library file
file_experiment = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_experiment.xml'
; Define institution library file
file_institute = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_institute.xml'
; Define model library file
file_model = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_model.xml'
; Define scenario library file
file_scenario = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_scenario.xml'

; Create a mapping vector for domain to frequency and realm
domain_to_frequency_realm = [ [ 'fx', 'fx', '' ], [ 'Amon', 'mon', 'atmos' ], $
    [ 'Lmon', 'mon', 'land' ], [ 'Omon', 'mon', 'ocean' ], $
    [ 'OImon', 'mon', 'seaIce' ] ]

; Parse the output file name in order to determine the output properties
if n_elements( out_file ) ne 1 then stop
out_file_parsed = strsplit( out_file, '/', extract=1, count=n_out_file_parsed )
out_file_parsed = strsplit( out_file_parsed[n_out_file_parsed-1], '_', $
    extract=1, count=n_out_file_parsed )
if n_out_file_parsed ne 8 then stop

; Determine the output variables
out_var_label = strsplit( out_file_parsed[0], '-', extract=1, count=n_out_var )

; Determine the input variables
if not( keyword_set( in_var_label ) ) then begin
  in_var_label = ''
  for i_out_var = 0, n_out_var - 1 do begin
    temp = [ 'tos-sic', 'tos', 'sic', 'tosbcs-sicbcs', 'tosbcs', 'sicbcs' ]
    if ( max( out_var_label[i_out_var] eq temp ) eq 1 ) $
        and ( max( in_var_label eq 'tos' ) eq 0 ) then begin
      in_var_label = [ in_var_label, 'tos' ]
    endif
    temp = [ 'tos-sic', 'sic', 'tosbcs-sicbcs', 'sicbcs' ]
    if ( max( out_var_label[i_out_var] eq temp ) eq 1 ) $
        and ( max( in_var_label eq 'sic' ) eq 0 ) then begin
      in_var_label = [ in_var_label, 'sic' ]
    endif
  endfor
endif
id = where( in_var_label ne '', n_in_var )
if n_in_var eq 0 then stop
in_var_label = in_var_label[id]

; Determine the output realm and frequency
out_domain = out_file_parsed[1]
id = where( domain_to_frequency_realm[0,*] eq out_domain, n_id )
if n_id ne 1 then stop
out_frequency = domain_to_frequency_realm[1,id[0]]
out_realm = domain_to_frequency_realm[1,id[0]]

; Determine the output model_id
out_model_id = out_file_parsed[2]

; Determine the output experiment_family, experiment, and subexperiment
out_experiment_family = out_file_parsed[3]
out_experiment = out_file_parsed[4]
out_subexperiment = out_file_parsed[5]
out_scenario = out_experiment_family + '_' + out_experiment + '_' $
    + out_subexperiment

; Determine the output run_id
out_run_id = out_file_parsed[6]

; Determine the output period
out_period = strsplit( out_file_parsed[7], '.', extract=1, count=temp )
if temp ne 2 then stop
out_period = strsplit( out_period[0], '-', extract=1, count=temp )
if temp ne 2 then stop
if strlen( out_period[0] ) ne 6 then stop
if strlen( out_period[1] ) ne 6 then stop

; Determine the input period
if not( keyword_set( in_period ) ) then in_period = out_period

; Determine the minimum input period.
; If a *bcs variable is requested then this is minimum 2 months padded to 
; either end of in_period)
in_period_minrange = in_period
if max( strpos( out_var_label, 'bcs' ) ge 0 ) eq 1 then begin
  temp_month = fix( strmid( in_period[0], 4, 2 ) ) - 2
  temp_year = fix( strmid( in_period[0], 0, 4 ) )
  if temp_month lt 1 then begin
    temp_year = temp_year - 1
    temp_month = temp_month + 12
  endif
  in_period_minrange[0] = str( temp_year ) $
      + str( temp_month, length=2, filler='0' ) + '01'
  temp_month = fix( strmid( in_period[1], 4, 2 ) ) + 2
  temp_year = fix( strmid( in_period[1], 0, 4 ) )
  if temp_month gt 12 then begin
    temp_year = temp_year + 1
    temp_month = temp_month - 12
  endif
  in_period_minrange[1] = str( temp_year ) $
      + str( temp_month, length=2, filler='0' )
endif

; Determine the maximum preferred input period to attempt to have.
; If a *bcs variable is requested then one year is padded to either end of 
; out_period plus additional months such that the start/end are a 
; January/December )
in_period_maxrange = in_period
if max( strpos( out_var_label, 'bcs' ) ge 0 ) eq 1 then begin
  temp_year = fix( strmid( in_period[0], 0, 4 ) )
  in_period_maxrange[0] = str( temp_year - 1 ) + '01'
  temp_year = fix( strmid( in_period[1], 0, 4 ) )
  in_period_maxrange[1] = str( temp_year + 1 ) + '12'
endif

;***********************************************************************
; Settings for various source/model_id combinations

; Initialise the global attribute list
out_global_attribute_list = ''

; Determine model-specific attributes and properties
temp_select_headers = [ 'attribute_model_name', 'attribute_institute_id', $
    'attribute_contact', 'attribute_acknowledgement', 'file_sftlf', $
    'attribute_license', 'tossic_source' ]
temp_select_values = [ 'attribute_model_id=' + out_model_id ]
markup_read, file_model, comment_char=';', select_values=temp_select_values, $
    select_headers=temp_select_headers, settings=temp_out
; Copy the model name attribute
id = where( temp_select_headers eq 'attribute_model_name', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    out_global_attribute_list = [ out_global_attribute_list, $
        'model_name=' + temp_out[id[0]] + '=7' ]
  endif
endif
; Copy the institute_id attribute
id = where( temp_select_headers eq 'attribute_institute_id', n_id )
if n_id eq 0 then begin
  print, 'ERROR c20c_dtos_v2_make_tossic_driver:  For some reason I have not ' $
      + 'looked for the institute_id attribute'
  stop
endif else begin
  if temp_out[id[0]] eq '' then begin
    print, 'ERROR c20c_dtos_v2_make_tossic_driver:  Unable to determine ' $
        + 'institute_id attribute for model ' + out_model_id + '.'
    stop
  endif else begin
    out_institute_id = temp_out[id[0]]
    out_global_attribute_list = [ out_global_attribute_list, $
        'institute_id=' + out_institute_id  + '=7' ]
  endelse
endelse
; Copy the contact attribute
id = where( temp_select_headers eq 'attribute_contact', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    if keyword_set( attribute_contact ) then begin
      out_global_attribute_list = [ out_global_attribute_list, $
          'contact=' + temp_out[id[0]]  + ', ' + attribute_contact + '=7' ]
    endif else begin
      out_global_attribute_list = [ out_global_attribute_list, $
          'contact=' + temp_out[id[0]]  + '=7' ]
    endelse
  endif else begin
   if keyword_set( attribute_contact ) then begin
      out_global_attribute_list = [ out_global_attribute_list, $
            'contact=' + attribute_contact + '=7' ]
    endif
  endelse
endif else begin
  if keyword_set( attribute_contact ) then begin
    out_global_attribute_list = [ out_global_attribute_list, $
          'contact=' + attribute_contact + '=7' ]
  endif
endelse
; Copy the acknowledgement attribute
id = where( temp_select_headers eq 'attribute_acknowledgement', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    out_global_attribute_list = [ out_global_attribute_list, $
        'acknowledgement=' + temp_out[id[0]] + '=7' ]
  endif
endif
; Copy the license attribute
id = where( temp_select_headers eq 'attribute_license', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    out_global_attribute_list = [ out_global_attribute_list, $
        'license=' + temp_out[id[0]] + '=7' ]
  endif
endif
; Record the grid file
if not( keyword_set( grid_file ) ) then begin
  id = where( temp_select_headers eq 'file_sftlf', n_id )
  if n_id ge 1 then begin
    if temp_out[id[0]] ne '' then grid_file = temp_out[id[0]]
  endif
endif
; Record the source observational product
if not( keyword_set( in_obs_label ) ) then begin
  id = where( temp_select_headers eq 'tossic_source', n_id )
  if n_id ge 1 then begin
    if temp_out[id[0]] ne '' then in_obs_label = temp_out[id[0]]
  endif
endif

; Determine institute-specific attributes and properties
temp_select_headers = [ 'attribute_institution' ]
temp_select_values = [ 'attribute_institute_id=' + out_institute_id ]
markup_read, file_institute, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
; Copy the institution attribute
id = where( temp_select_headers eq 'attribute_institution', n_id )
if n_id eq 0 then begin
  print, 'ERROR c20c_dtos_v2_make_tossic_driver:  For some reason I have not ' $
      + 'looked for the institution attribute'
  stop
endif else begin
  if temp_out[id[0]] eq '' then begin
    print, 'ERROR c20c_dtos_v2_make_tossic_driver:  Unable to determine ' $
        + 'institution attribute for institute ' + out_institute_id + '.'
    stop
  endif else begin
    out_global_attribute_list = [ out_global_attribute_list, $
        'institution=' + temp_out[id[0]] + '=7' ]
  endelse
endelse

; Determine model-scenario-specific attributes and properties
temp_select_headers = [ 'tossic_source', 'attribute_title' ]
temp_select_values = [ 'attribute_model_id=' + out_model_id, $
    'attribute_experiment_family_experiment_subexperiment=' + out_scenario ]
markup_read, file_scenario, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
; Copy the source observational product
if not( keyword_set( in_obs_label ) ) then begin
  if n_id eq 0 then begin
    print, 'ERROR c20c_dtos_v2_make_tossic_driver:  For some reason I have ' $
      + 'not looked for the original source tos and sic data.'
    stop
  endif else begin
    if temp_out[id[0]] eq '' then begin
      print, 'ERROR c20c_dtos_v2_make_tossic_driver:  Unable to determine ' $
          + 'original source tos and sic data for mode ' + out_model_id $
          + 'run under scenario ' + out_scenario + '.'
      stop
    endif else begin
      in_obs_label = temp_out[id[0]]
    endelse
  endelse
endif
; Copy the title attribute
id = where( temp_select_headers eq 'attribute_title', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    out_global_attribute_list = [ out_global_attribute_list, $
        'title=' + temp_out[id[0]] + '=7' ]
  endif
endif

; Determine "experiment"-specific attributes and properties
temp_select_headers = [ 'attribute_acknowledgement', 'attribute_comment', $
    'delta_file', 'delta_factor', 'delta_var_label', 'tossic_fit_method', $
    'tossic_fit_period' ]
temp_select_values = [ $
    'attribute_experiment_family_experiment=' + out_experiment_family + '_' $
    + out_experiment ]
markup_read, file_experiment, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
; Copy the acknowledgement attribute
id = where( temp_select_headers eq 'attribute_acknowledgement', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    temp = strpos( out_global_attribute_list, 'acknowledgement' )
    id = where( temp eq 0, n_id )
    if n_id ge 1 then begin
      temp = strsplit( out_global_attribute_list[id[0]], '=', extract=1 )
      out_global_attribute_list[id[0]] = temp[0] + '=' + temp[1] + '  ' $
          + temp_out[id[0]] + '=7'
    endif else begin
      out_global_attribute_list = [ out_global_attribute_list, $
          'acknowledgement=' + temp_out[id[0]] + '=7' ]
    endelse
  endif
endif
; Copy the comment attribute
id = where( temp_select_headers eq 'attribute_comment', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    temp = strpos( out_global_attribute_list, 'comment=' )
    id = where( temp eq 0, n_id )
    if n_id ge 1 then begin
      temp = strsplit( out_global_attribute_list[id[0]], '=', extract=1 )
      out_global_attribute_list[id[0]] = temp[0] + '=' + temp[1] + '  ' $
          + temp_out[id[0]] + '=7'
    endif else begin
      out_global_attribute_list = [ out_global_attribute_list, $
          'comment=' + temp_out[id[0]] + '=7' ]
    endelse
  endif
endif
; Record the file containing the deltaSST data
if not( keyword_set( delta_file ) ) then begin
  id = where( temp_select_headers eq 'delta_file', n_id )
  if n_id ge 1 then begin
    if temp_out[id[0]] ne '' then begin
      delta_file = strsplit( temp_out[id[0]], '&', extract=1, count=n_delta )
    endif
  endif
endif
; Record the factor by which to multiply the deltaSST data when adding to the 
; input  data
id = where( temp_select_headers eq 'delta_factor', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    delta_factor = float( strsplit( temp_out[id[0]], '&', extract=1 ) )
    if n_elements( delta_factor ) ne n_delta then stop
  endif
endif
; Record the label of the variable to read from the deltaSST file
id = where( temp_select_headers eq 'delta_var_label', n_id )
if n_id ge 1 then begin
  if temp_out[id[0]] ne '' then begin
    delta_var_label = strsplit( temp_out[id[0]], '&', extract=1 )
    if n_elements( delta_var_label ) ne n_delta then stop
  endif
endif
; Record the method used for developing the sea ice adjustment function
if keyword_set( delta_file ) then begin
  id = where( temp_select_headers eq 'tossic_fit_method', n_id )
  if n_id ge 1 then begin
    if temp_out[id[0]] ne '' then begin
      tossic_fit_method = strsplit( temp_out[id[0]], ',', extract=1 )
    endif
  endif
endif
; Record the time period to use for fitting
if keyword_set( delta_file ) then begin
  id = where( temp_select_headers eq 'tossic_fit_period', n_id )
  if n_id ge 1 then begin
    if temp_out[id[0]] ne '' then begin
      tossic_fit_period = strsplit( temp_out[id[0]], ',', extract=1, $
          count=n_temp )
      if n_temp ne 2 then stop
    endif
  endif
endif
; Set the variable to which the delta should be applied
if keyword_set( delta_var_label ) then begin
  delta_apply_var_label = 'tos'
endif

; The default list of attributes to include for the output variables
out_var_attribute_list_default = [ 'long_name==7', 'standard_name==7', $
    'units==7', 'missing_value=1.e+20=4', '_FillValue=1.e+20=4' ]
if n_out_var eq 1 then begin
  out_var_attribute_list = out_var_attribute_list_default
endif else begin
  out_var_attribute_list = add_dim( out_var_attribute_list_default, 1, $
      n_out_var )
endelse

; Make modifications for producing a file that is ready-to-use by the intended 
; model, rather than one strictly following C20C+ D&A protocols
if keyword_set( ready_to_use_opt ) then begin
  ; For the CAM family of models
  temp = [ 'CAM5-1-2degree', 'CAM5-1-1degree', 'CAM5-1-025degree', $
      'CAM5-1-2-025degree', 'CAM4-2degree' ]
  if max( temp eq out_model_id ) eq 1 then begin
    ; Then include non-Sheng-Zwiers adjusted variables if appropriate
    if ( max( out_var_label eq 'tosbcs' ) eq 1 ) $
        and ( max( out_var_label eq 'tos' ) eq 0 ) then begin
      out_var_label = [ out_var_label, 'tos' ]
      n_out_var = n_out_var + 1
      out_var_attribute_list = [ [ out_var_attribute_list ], $
          [ out_var_attribute_list_default ] ]
    endif
    if ( max( out_var_label eq 'sicbcs' ) eq 1 ) $
        and ( max( out_var_label eq 'sic' ) eq 0 ) then begin
      out_var_label = [ out_var_label, 'sic' ]
      n_out_var = n_out_var + 1
      out_var_attribute_list = [ [ out_var_attribute_list ], $
          [ out_var_attribute_list_default ] ]
    endif
    ; The list of attributes to include for the output variables.
    ; Note that the "label" attribute is in fact not an attribute but a command 
    ; to rename the variable in the output NetCDF file.
    for i_out_var = 0, n_out_var - 1 do begin
      temp_attribute = ''
      if out_var_label[i_out_var] eq 'tos' then begin
        temp_attribute = [ 'label=SST_cpl_prediddle', 'units=deg_C' ]
      endif else if out_var_label[i_out_var] eq 'sic' then begin
        temp_attribute = [ 'label=ice_cov_prediddle', 'units=fraction' ]
      endif else if out_var_label[i_out_var] eq 'tosbcs' then begin
        temp_attribute = [ 'label=SST_cpl', 'units=deg_C', 'comment' ]
      endif else if out_var_label[i_out_var] eq 'sicbcs' then begin
        temp_attribute = [ 'label=ice_cov', 'units=fraction', 'comment' ]
      endif else begin
        stop
      endelse
      if temp_attribute[0] ne '' then begin
        for i_attribute = 0, n_elements( temp_attribute ) - 1 do begin
          temp = strsplit( temp_attribute[i_attribute], '=', extract=1 )
          temp = strpos( out_var_attribute_list[*,i_out_var], temp[0] + '=' )
          id = where( temp eq 0, n_id )
          if n_id gt 1 then stop
          if n_id eq 1 then begin
            out_var_attribute_list[id[0],i_out_var] $
                = temp_attribute[i_attribute]
          endif else begin
            id = where( out_var_attribute_list[*,i_out_var] eq '', n_id )
            if n_id eq 0 then begin
              out_var_attribute_list = [ out_var_attribute_list, $
                  strarr( 1, n_out_var ) ]
              id = n_elements( out_var_attribute_list[*,0] ) - 1
            endif
            out_var_attribute_list[id[0],i_out_var] $
                = temp_attribute[i_attribute]
          endelse
        endfor
      endif
    endfor
    ; Add the date and datesec time variables
    out_time_label = [ 'time', 'date', 'datesec' ]
    out_time_attribute_list = strarr( 1, 3 )
    out_time_attribute_list[0,0] = 'calendar=365_day'
  endif else begin
    stop
  endelse
endif

; Determine period of availability of NOAA-OI-v2
; (Note we are extending input4MIPs-NOAA-OI-v2 ourselves with NOAA-OI-v2, 
; rather than waiting for input4MIPs to do so.)
if ( in_obs_label eq 'NOAA-OI-v2' ) $
    or ( in_obs_label eq 'input4MIPs-NOAA-OI-v2' ) then begin
  ; Determine the current date
  date_today = bin_date( systime() )
  date_previous = date_today
  date_today = str( date_today[0] ) $
      + str( date_today[1], length=2, filler='0' ) $
      + str( date_today[2], length=2, filler='0' )
  ; Determine the last day of the previous month
  date_previous[1] = date_previous[1] - 1
  if date_previous[1] lt 1 then begin
    date_previous[0] = date_previous[0] - 1
    date_previous[1] = date_previous[1] + 12
  endif
  if ( date_previous[0] mod 4 eq 0 ) $
      and ( ( date_previous[1] mod 100 ne 0 ) $
      or ( date_previous[1] mod 400 eq 0 ) ) then begin
    temp = month_day( date_previous[1] - 1, leapyear=1 )  
  endif else begin
    temp = month_day( date_previous[1] - 1, leapyear=0 )  
  endelse
  date_previous[2] = temp[1] - temp[0] + 1
  date_previous = str( date_previous[0] ) $
      + str( date_previous[1], length=2, filler='0' ) $
      + str( date_previous[2], length=2, filler='0' )
  ; Determine the day of the week of the last day of the previous month
  spawn, 'date --date=' + date_previous, date_previous_dayofweek
  date_previous_dayofweek = strmid( date_previous_dayofweek[0], 0, 3 )
  ; Generate the advertised timetable for availability
  temp = [ [ 'Sun', '08' ], [ 'Mon', '07' ], [ 'Tue', '06' ], [ 'Wed', '05' ], $
      [ 'Thu', '11' ], [ 'Fri', '10' ], [ 'Sat', '09' ] ]
  ; Determine the day when that month's data will be available
  id = where( temp[0,*] eq date_previous_dayofweek, n_id )
  if n_id ne 1 then stop
  date_available = strmid( date_today, 0, 6 ) + temp[1,id[0]]
  ; If that day has arrived then take the previous month as available
  if long( date_today ) gt long( date_available ) then begin
    date_available = strmid( date_previous, 0, 6 )
  ; If that day has not yet arrived then interpret the "previous" month as 
  ; being the one before that
  endif else begin
    date_available = fix( [ strmid( date_available, 0, 4 ), $
        strmid( date_available, 4, 2 ) ] )
    date_available[1] = date_available[1] - 1
    if date_available[1] lt 1 then begin
      date_available[0] = date_available[0] - 1
      date_available[1] = date_available[1] + 12
    endif
    date_available = str( date_available[0] ) $
        + str( date_available[1], length=2, filler='0' )
  endelse
endif

; Restrict global attribute list to defined values
id = where( strlen( out_global_attribute_list ) gt 0 )
out_global_attribute_list = out_global_attribute_list[id]

;***********************************************************************
; Settings for the experiment_family/experiment/subexperiment combinations

; The default directories for data on NERSC systems
if keyword_set( nersc_opt ) then begin
  ; The CASCADE project directory
  dir_c20c = '/project/projectdirs/m1517/'
  ; The deltaSST directory
  if keyword_set( delta_file ) then begin
    if strpos( delta_file[0], '/' ) lt 0 then begin
      for i_delta = 0, n_delta - 1 do begin
        temp_delta_file = strsplit( delta_file[i_delta], ',', extract=1, $
            count=temp_n_delta_file )
        temp = strarr( temp_n_delta_file, 8 )
        for i_file = 0, temp_n_delta_file - 1 do begin
          temp_delta = strsplit( temp_delta_file[i_file], '_', extract=1, $
              count=n_temp )
          if n_temp ne 8 then stop
          temp[i_file,*] = temp_delta
        endfor
        temp_delta_dir = dir_c20c + 'C20C/C20C/' + temp[*,2] + '/' + temp[*,3] $
            + '/' + temp[*,4] + '/' + temp[*,5]
        for i_file = 0, temp_n_delta_file - 1 do begin
          id = where( domain_to_frequency_realm[0,*] eq temp[i_file,1], n_id )
          if n_id ne 1 then stop
          temp_delta_dir[i_file] = temp_delta_dir[i_file] + '/' $
              + domain_to_frequency_realm[1,id[0]] + '/' $
              + domain_to_frequency_realm[2,id[0]]
        endfor
        temp_delta_dir = temp_delta_dir + '/' + temp[*,0] + '/' + temp[*,6] $
            + '/'
        temp_delta_file = temp_delta_dir + temp_delta_file
        if temp_n_delta_file gt 1 then begin
          temp_delta_file = strjoin( temp_delta_file, ',' )
        endif
        delta_file[i_delta] = temp_delta_file
      endfor
    endif
  endif
  ; The grid directory if this is to be interpolated to a model grid
  if keyword_set( grid_file ) then begin
    if strpos( grid_file, '/' ) lt 0 then begin
      temp = strsplit( grid_file, '_', extract=1, count=n_temp )
      if n_temp ne 8 then stop
      grid_dir = dir_c20c + 'C20C/' + out_institute_id + '/' + temp[2] + '/' $
          + temp[3] + '/' + temp[4] + '/' + temp[5]
      if temp[1] eq 'fx' then begin
        grid_dir = grid_dir + '/fx'
        if max( temp[0] eq [ 'areacella', 'orog', 'sftlf' ]) eq 1 then begin
          grid_dir = grid_dir + '/atmos/' + temp[0]
        endif else begin
          stop
        endelse
      endif else begin
        id = where( domain_to_frequency_realm[0,*] eq temp[1], n_id )
        if n_id ne 1 then stop
        grid_dir = grid_dir + '/' + domain_to_frequency_realm[1,id[0]] $
            + '/' + domain_to_frequency_realm[2,id[0]] + temp[0]
      endelse
      grid_dir = grid_dir + '/' + temp[6] + '/'
      grid_file = grid_dir + grid_file
    endif
  endif
  ; The input observation data
  markup_read, file_model, comment_char=';', $
      select_values='attribute_model_id='+in_obs_label, $
      select_headers='attribute_institute_id', settings=in_institute_id
  in_institute_id = in_institute_id[0]
  if in_institute_id[0] eq '' then stop
  if max( in_var_label eq 'tos' ) eq 1 then begin
    in_file = dir_c20c + 'daithi/data/' + in_institute_id + '/' + in_obs_label $
        + '/All-Hist/est1/v1-0/mon/ocean/tos/observed/*.nc'
  endif
  if max( in_var_label eq 'sic' ) eq 1 then begin
    temp_in_file = dir_c20c + 'daithi/data/' + in_institute_id + '/' $
        + in_obs_label + '/All-Hist/est1/v1-0/mon/seaIce/sic/observed/*.nc'
    if keyword_set( in_file ) then begin
      in_file = [ in_file, temp_in_file ]
    endif else begin
      in_file = temp_in_file
    endelse
  endif
; Otherwise ensure the input files have been defined
endif else begin
  ; We need input data files
  if not( keyword_set( in_file ) ) then begin
    print, 'ERROR c20c_dtos_v2_make_tossic_driver:  No IN_FILE input data ' $
        + 'files defined.'
    stop
  endif
  ; We need a grid file if interpolation is to be performed
  if ( out_model_id ne in_obs_label ) and not( keyword_set( grid_file ) ) $
      then begin
    print, 'ERROR c20c_dtos_v2_make_tossic_driver:  Output model differs ' $
        + 'from input source, but no GRID_FILE defined for interpolation.'
    stop
  endif
endelse

;***********************************************************************
; Determine, and possibly improve, availability of input data

; Generate the full list of files consistent with in_file
in_file_list = file_search( in_file, count=n_in_file_list )
if n_in_file_list eq 0 then begin
  in_file_list = in_file + [ '/*.nc', '/*.nc4' ]
  in_file_list = file_search( in_file_list, count=n_in_file_list )
  if n_in_file_list eq 0 then stop
endif

; Initialise list of files for each variable
in_var_file = strarr( n_in_var )
; Iterate through input files
for i_in_file = 0, n_in_file_list - 1 do begin
  ; Get the list of variables available in this file
  temp_var_label = netcdf_read_metadata( in_file_list[i_in_file] )
  ; Determine which input variables are provided by this file
  index = isin( temp_var_label, in_var_label )
  id = where( index eq 1, n_id )
  ; Add this file to the lists for the available variables
  for i_id = 0, n_id - 1 do begin
    if in_var_file[id[i_id]] eq '' then begin
      in_var_file[id[i_id]] = in_file_list[i_in_file]
    endif else begin
      in_var_file[id[i_id]] = in_var_file[id[i_id]] + ',' $
          + in_file_list[i_in_file]
    endelse
  endfor
endfor
; Ensure we have files for all input variables
if min( strlen( in_var_file ) ) eq 0 then stop

; Generate list of files for sea ice adjustment calculations
if keyword_set( tossic_fit_method ) then begin
  if max( tossic_fit_method eq 'stonepall2018' ) eq 1 then begin
    id = where( in_var_label eq 'tos', n_id )
    if n_id ne 1 then stop
    tossic_fit_tos_file = in_var_file[id[0]]
    id = where( in_var_label eq 'sic', n_id )
    if n_id ne 1 then stop
    tossic_fit_sic_file = in_var_file[id[0]]
  endif
endif

; Iterate through input variables
for i_in_var = 0, n_in_var - 1 do begin
  ; Determine if the preferred period is covered by the files available for 
  ; this variable
  temp_period = ''
  temp = netcdf_read_geo_multitime( in_var_file[i_in_var], '', $
      period_time=temp_period, fix_time=1 )
  temp_period = strmid( temp_period, 0, 6 )
  ; Determine if we would like to try to retrieve more files for the end of 
  ; the period
  if temp_period[1] lt in_period_maxrange[1] then begin
    ; Determine the directories to put the downloaded data
    if not( keyword_set( dir_dest ) ) then begin
      dir_dest = strarr( n_in_var )
      for i_in_var_1 = 0, n_in_var - 1 do begin
        temp_dir = strsplit( in_var_file[i_in_var_1], '/', extract=1, $
            count=n_temp_dir )
        if n_temp_dir ge 2 then begin
          temp_dir = string_from_vector( temp_dir[0:n_temp_dir-2], spacer='/', $
              nospace=1 )
        endif else begin
          temp_dir = '.'
        endelse
        dir_dest[i_in_var_1] = temp_dir
      endfor
    endif
    ; If the source data is from HadISST1
    if in_obs_label eq 'HadISST1' then begin
      ; If we only need to update an amount less than the last 12 months
      temp_year = 12 * ( fix( strmid( in_period_maxrange[1], 0, 4 ) ) $
          - fix( strmid( temp_period[1], 0, 4 ) ) )
      temp_month = fix( strmid( in_period_maxrange[1], 4, 2 ) ) $
          - fix( strmid( temp_period[1], 4, 2 ) ) + 12 * temp_year
      flag_get_full = 1
      if temp_year lt 12 then begin
        ; Then we should try getting the data from the smaller "update" file
        url_list = 'https://www.metoffice.gov.uk/hadobs/hadisst/data/'
        if in_var_label[i_in_var] eq 'tos' then begin
          url_file = 'HadISST1_SST_update.nc.gz'
        endif else if in_var_label[i_in_var] eq 'sic' then begin
          url_file = 'HadISST1_ICE_update.nc.gz'
        endif else begin
          stop
        endelse
        url_list = url_list + url_file
        c20c_dtos_v2_download_obs, url_list, var_label=in_var_label, $
            dir_dest=dir_dest[0]+'/download'
        ; Determine the uncompressed file name
        pos = strpos( url_file, '.gz' )
        if pos ge 0 then url_file = strmid( url_file, 0, pos )
        ; Load the time vector
        temp = netcdf_read_geo( dir_dest[0] + '/download/' + url_file, '', $
            time=time_download )
        ; If the update does not start early enough
        if max( strmid( time_download, 0, 6 ) ) lt temp_period[1] then begin
          ; The "update" file is no use to us, so delete it
          spawn, 'rm dir_dest[0]/download/' + url_file
        ; If we have enough here
        endif else begin
          ; Flag that we do not need to retrieve the longer file
          flag_get_full = 0
        endelse
      endif
      ; If we need to download the complete file
      if flag_get_full eq 1 then begin
        ; Define the URL
        url_list = 'https://www.metoffice.gov.uk/hadobs/hadisst/data/'
        if in_var_label[i_in_var] eq 'tos' then begin
          url_file = 'HadISST_sst.nc.gz'
        endif else if in_var_label[i_in_var] eq 'sic' then begin
          url_file = 'HadISST_ice.nc.gz'
        endif else begin
          stop
        endelse
        url_list = url_list + url_file
        c20c_dtos_v2_download_obs, url_list, var_label=in_var_label, $
            dir_dest=dir_dest[0]+'/download'
        ; Determine the uncompressed file name
        pos = strpos( url_file, '.gz' )
        if pos ge 0 then url_file = strmid( url_file, 0, pos )
        ; Load the time vector
        temp = netcdf_read_geo( dir_dest[0] + '/download/' + url_file, '', $
            time=time_download )
      endif
      ; Identify the time steps we want to extract
      id_time = where( strmid( time_download, 0, 6 ) gt temp_period[1], $
          n_id_time )
      ; If we are to include the entire file
      if n_id_time eq n_elements( time_download ) then begin
        ; Move the entire file to the appropriate directory
        pos = strpos( url_file, '.nc' )
        temp_command = 'mv ' + dir_dest[0] + '/download/' + url_file + ' ' $
            + dir_dest[0] + '/' + strmid( url_file, 0, pos ) + '_' $
            + strmid( time_download[0], 0, 6 ) + '-' $
            + strmid( time_download[n_id_time-1], 0, 6 ) + '.nc'
        spawn, temp, error_status=temp_status
        if temp_status ne 0 then stop
      ; If we are to include only a portion of the file
      endif else if n_id_time gt 0 then begin
        ; Extract the desired months
        pos = strpos( url_file, '.nc' )
        temp_command = 'ncks -d time,' + str( id_time[0] ) + ',' $
            + str( id_time[n_id_time-1] ) + ' ' + dir_dest[0] + '/download/' $
            + url_file + ' ' + dir_dest[0] + '/' $
            + strmid( url_file, 0, pos ) + '_' $
            + strmid( time_download[id_time[0]], 0, 6 ) + '-' $
            + strmid( time_download[id_time[n_id_time-1]], 0, 6 ) + '.nc'
        spawn, temp, error_status=temp_status
        if temp_status ne 0 then stop
        ; Delete the file as we are done with it
        spawn, 'rm dir_dest[0]/download/' + url_file
      ; This file is no use to us
      endif else begin
        ; So delete it
        spawn, 'rm dir_dest[0]/download/' + url_file
      endelse
    ; If the source data is from NOAA-OI-v2 (including input4MIPs, which just 
    ; copies it)
    endif else if ( in_obs_label eq 'NOAA-OI-v2' ) or $
        ( in_obs_label eq 'input4MIPs-NOAA-OI-v2' ) then begin
      ; Check if more are available
      if long( temp_period[1] ) gt long( date_available  ) then stop
      ; Initialise list of months to retrieve
      url_list = ''
      ; Iterate through months to retrieve
      date_doing_year = fix( strmid( temp_period[1], 0, 4 ) )
      date_doing_month = fix( strmid( temp_period[1], 4, 2 ) )
      date_done = min( [ date_available, in_period_maxrange[1] ] )
      while ( date_doing_year ne fix( strmid( date_done, 0, 4 ) ) ) $
          or ( date_doing_month ne fix( strmid( date_done, 4, 2 ) ) ) $
          do begin
        ; Increment to the next month
        date_doing_month = date_doing_month + 1
        if date_doing_month gt 12 then begin
          date_doing_year = date_doing_year + 1
          date_doing_month = date_doing_month - 1
        endif
        ; Build URL list
        temp_url = 'oiv2mon. ' + str( date_doing_year ) $
            + str( date_doing_month, length=2, filler='0' ) + '.grb'
        temp_url = 'ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/GRIB/' $
            + temp_url
        if url_list[0] eq '' then begin
          url_list[0] = temp_url
        endif else begin
          url_list = [ url_list, temp_url ]
        endelse
      endwhile
      ; Download data
      if keyword_set( url_list ) then begin
        c20c_dtos_v2_download_obs, url_list, var_label=in_var_label, $
            dir_dest=dir_dest, cdo_shifttime='+15days'
      endif
    endif else begin
      stop
    endelse
  endif
  ; Determine if the minimum period is covered by the files available for this 
  ; variable
  temp_period = ''
  temp = netcdf_read_geo_multitime( in_var_file[i_in_var], '', $
      period_time=temp_period, fix_time=1 )
  temp_period = strmid( temp_period, 0, 6 )
  ; Determine if the period is adequate
  if long( temp_period[0] ) gt long( in_period_minrange[0] ) then stop
  if long( temp_period[1] ) lt long( in_period_minrange[1] ) then stop
  ; Determine the input period to use
  in_period = [ max( [ temp_period[0], in_period_maxrange[0] ] ), $
      min( [ temp_period[1], in_period_maxrange[1] ] ) ]
endfor

;***********************************************************************
; Generate the output NetCDF file

; Perform surface boundary condition calculations and create NetCDF file
in_var_file = reform( in_var_file, 1, n_in_var )
c20c_dtos_v2_make_tossic, in_file=in_var_file, in_var_label=in_var_label, $
    out_file=out_file, out_var_label=out_var_label, $
    out_var_attribute_list=out_var_attribute_list, $
    out_global_attribute_list=out_global_attribute_list, $
    out_global_attribute_c20c=1, cf_standard=1, out_period=out_period, $
    delta_file=delta_file, delta_factor=delta_factor, $
    delta_var_label=delta_var_label, in_period=in_period, $
    tossic_fit_method=tossic_fit_method, tossic_fit_period=tossic_fit_period, $
    tossic_fit_tos_file=tossic_fit_tos_file, $
    tossic_fit_sic_file=tossic_fit_sic_file, v1=v1_opt, $
    grid_file=grid_file, grid_var_label=grid_var_label, $
    out_time_label=out_time_label, delta_period=delta_period, $
    delta_extend_cycle=delta_extend_cycle_opt, $
    out_time_attribute_list=out_time_attribute_list, $
    delta_apply_var_label=delta_apply_var_label

;***********************************************************************
; The end

;stop
return
END


