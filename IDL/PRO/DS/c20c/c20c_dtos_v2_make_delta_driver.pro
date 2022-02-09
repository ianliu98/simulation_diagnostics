;+
; NAME:
;    c20c_dtos_v2_make_delta_driver
;
; PURPOSE:
;    This procedure runs the specific calls of c20c_dtos_v2_make_delta.pro 
;    required to produce deltaSST fields used by the C20C+ Detection and 
;    Attribution Project, the Weather Risk Attribution Forecast, and the HAPPI 
;    project..
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_make_delta_driver, scenario
;
; INPUTS:
;    SCENARIO:  A scalar string describing the delta sea surface temperature 
;        estimate to perform, in the format 
;        "<experiment_family>/<experiment>/<subexperiment>".  For instance, the 
;        benchmark estimate for C20C+ D&A described in Stone and Pall (2018) 
;        is defined by 'Nat-Hist/CMIP5-est1/v1-0'.
;    CONTACT, DEP_FILE, DEP_VAR, GRID_FILE, IN_DIR, INDEP_PERIOD, INSTITUTE_ID, 
;        INSTITUTION, OPTIONS, OUT_PERIOD, TIME_FILTER, TITLE
;
; KEYWORD PARAMETERS:
;    CONTACT:  An optional scalar string containing the name and/or contact 
;        details of the person running the code.
;    DEP_FILE:  A string vector listing the files containing the dependent data 
;        if a regression analysis is being used to calibrate the attributable 
;        warming estimate.  This may include the "*" wildcard.  It should 
;        include the directory as part of the file name.  This is required if 
;        the regression approach is used.  **NOT YET IMPLEMENTED**.
;    DEP_VAR:  A scalar string containing the label of the variable to read 
;        DEPENDENT_FILE.  This is required if the regression approach is used.
;        **NOT YET IMPLEMENTED**.
;    GRID_FILE:  An optional scalar string naming the file containing the 
;        target longitude-latitude grid  for interpolation of the output.
;    IN_DIR:  An optional string vector listing the directories in which to 
;        search for the files listed in c20c_dtos_v2_map_delta_source.xml 
;        which contain the input data for estimating the deltaSST scenario.  
;        The default is the present directory.
;    INDEP_PERIOD:  A required 2-element or 2,N_INDEP string array describing 
;        the start and end dates of the period over which to estimate the 
;        attributable delta.  Of format [<start>,<end>] with the dates in the 
;        'yyyymm' format.  If of size 2*N_INDEP, where N_INDEP is the number 
;        of independent scenarios to be used in estimating the attributable 
;        difference, this can be used for instance if the estimate is based on 
;        the change between two time periods in the observed record.
;    INSTITUTE_ID:  An optional scalar string providing the value for the 
;        global institute_id attribute in the output NetCDF file.  This 
;        overrides any default listed in c20c_dtos_v2_map_delta_scenario.xml.
;    INSTITUTION:  An optional scalar string providing the value for the global 
;        institution attribute in the output NetCDF file.  The default is 
;        the setting in c20c_dtos_v2_map_institute.xml for the INSTITUTE_ID 
;        value.
;    OPTIONS:  The OPTIONS input for c20c_dots_v2_make_delta.pro.  See that 
;        procedure for details.
;    OUT_PERIOD:  An optional two-element string vector containing the start 
;        and end dates respectively of the time period of the output data.  
;        Dates should be of the 'yyyymm' format.  The default is 
;        INDEP_PERIOD[*,0].
;    TIME_FILTER:  An optional two-element string vector with the first element 
;        describing the type of temporal filter to use and the second element 
;        describing the length of the filter.  For example ['boxcar','5'] would 
;        apply a boxcar filter of length 5.
;    TITLE:  An optional scalar string providing the global title attribute 
;        in the output NetCDF file.  A default can be constructed by the code 
;        instead.
;
; OUTPUTS:
;    ---
;
; USES:
;    c20c_dtos_v2_make_delta.pro
;    c20c_dtos_v2_map_delta_scenario.xml
;    c20c_dtos_v2_map_institute.xml'
;    markup_read.pro
;    str.pro
;    string_substitute.pro
;
; PROCEDURE:
;    This procedure uses c20c_dtos_v2_map_delta_scenario.xml to generate a call 
;    of c20c_dtos_v2_make_delta.pro in order to generate a standard estimate of 
;    attributable ocean warming.
;
; EXAMPLES:
;    ; The following script generates the C20C+ D&A Nat-Hist/CMIP5-est1/v1-0 benchmark deltaSST estimate on NERSC's Cori system.
;    dir_m1517 = '/project/projectdirs/m1517/'
;    c20c_dtos_v2_make_delta_driver, 'Nat-Hist/CMIP5-est1/v1-0', time_filter=['boxcar','10'], in_dir=dir_m1517+'daithi/data/CMIP5/*/'+['historical','rcp45','historicalNat']+'/est1/v1-0/mon/atmos/ts/r*', indep_period=['199001','201912'], contact='Daithi Stone (dastone@runbox.com)', institute_id='C20C', grid_file=dir_m1517+'daithi/data/NOAA-EMC-and-NCAR/NOAA-OI-v2/All-Hist/est1/v1-0/mon/ocean/tos/observed/oiv2mon.201101_tos.nc', out_period=['199001','201912'], options=['extend last year']
;    ; See c20c_dtos_v2_make_delta_driver_scripts.txt for further scripts.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2018-08-28
;    Modified:  DAS, 2018-11-07 (Set to output tos when ts input)
;    Modified:  DAS, 2018-11-15 (Included model "delta" in output file name)
;    Modified:  DAS, 2018-12-12 (Completed documentation.  Added 
;        out_global_attribute_c20c=1 in call to c20c_dtos_v2_make_delta.pro.)
;-

PRO C20C_DTOS_V2_MAKE_DELTA_DRIVER, $
    SCENARIO, $
    CONTACT=contact, $
    DEP_FILE=dep_file, DEP_VAR=dep_var, $
    GRID_FILE=grid_file, $
    IN_DIR=in_dir, $
    INDEP_PERIOD=indep_period, $
    INSTITUTE_ID=institute_id, INSTITUTION=institution, $
    OPTIONS=options, $
    OUT_PERIOD=out_period, $
    TIME_FILTER=time_filter, $
    TITLE=title

;***********************************************************************
; Check inputs

; Flag that DEP_FILE and DEP_VAR are not yet supported
if keyword_set( dep_file ) or keyword_set( dep_var ) then begin
  stop, 'c20c_dtos_v2_make_delta_driver.pro:  Abort.  DEP_FILE and DEP_VAR ' $
      + 'inputs are not implemented yet.'
endif

; Check that this is a supported scenario
if n_elements( scenario ) ne 1 then stop
temp = [ 'Nat-Hist/CMIP5-est1/v1-0', $
    'Nat-Hist/CCSM4-est1/v1-0', $
    'Nat-Hist/CESM1-CAM5-est1/v1-0' ]
if max( scenario eq temp ) ne 1 then stop

; The default input directory
if not( keyword_set( in_dir ) ) then in_dir = ''

;***********************************************************************
; Define constants and settings

; Determine the location of this code (and associated library files)
dir_c20c_dtos_v2 = file_which( 'c20c_dtos_v2_make_delta_driver.pro' )
if not( keyword_set( dir_c20c_dtos_v2 ) ) then stop
temp = strsplit( dir_c20c_dtos_v2, '/', extract=1, count=n_temp )
if n_temp eq 1 then begin
  dir_c20c_dtos_v2 = ''
endif else begin
  if strmid( dir_c20c_dtos_v2, 0, 1 ) eq '/' then temp[0] = '/' + temp[0]
  dir_c20c_dtos_v2 = strjoin( temp[0:n_temp-2], '/' ) + '/'
endelse

; Define the file name containing scenario definitions
file_scenario = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_delta_scenario.xml'

; Define institution library file
file_institute = dir_c20c_dtos_v2 + 'c20c_dtos_v2_map_institute.xml'

;***********************************************************************
; Build attribute list

; Initialise list of global and variable attributes
out_global_attribute_label = ''
out_global_attribute_value = ''
out_global_attribute_type = 7
out_var_attribute_label = ''
out_var_attribute_value = ''
out_var_attribute_type = 7
; Add the contact attribute
if keyword_set( contact ) then begin
  out_global_attribute_label = [ out_global_attribute_label, 'contact' ]
  out_global_attribute_value = [ out_global_attribute_value, contact ]
  out_global_attribute_type = [ out_global_attribute_type, 7 ]
endif
; Add the title attribute
if not( keyword_set( title ) ) then begin
  title = 'The ' + scenario + ' attributable warming estimate for the ' $
      + 'C20C+ Detection and Attribution project.'
endif
out_global_attribute_label = [ out_global_attribute_label, 'title' ]
out_global_attribute_value = [ out_global_attribute_value, title ]
out_global_attribute_type = [ out_global_attribute_type, 7 ]
; Add the institute_id attribute
if keyword_set( institute_id ) then begin
  out_global_attribute_label = [ out_global_attribute_label, 'institute_id' ]
  out_global_attribute_value = [ out_global_attribute_value, institute_id ]
  out_global_attribute_type = [ out_global_attribute_type, 7 ]
endif
; Add the institution attribute
if not( keyword_set( institution ) ) and keyword_set( institute_id ) then begin
  temp_select_headers = [ 'attribute_institution' ]
  temp_select_values = [ 'attribute_institute_id=' + institute_id ]
  markup_read, file_institute, comment_char=';', $
      select_values=temp_select_values, select_headers=temp_select_headers, $
      settings=temp_out
  id = where( temp_select_headers eq 'attribute_institution', n_id )
  if n_id gt 0 then institution = temp_out[0]
endif
if keyword_set( institution ) then begin
  out_global_attribute_label = [ out_global_attribute_label, 'institution' ]
  out_global_attribute_value = [ out_global_attribute_value, institution ]
  out_global_attribute_type = [ out_global_attribute_type, 7 ]
endif
; Add scenario-specific attributes
temp_select_headers = [ 'independent_variable', 'scenario_comment', $
    'scenario_domain' ]
temp_select_values = [ 'scenario_name=' + scenario ]
markup_read, file_scenario, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
id = where( temp_select_headers eq 'independent_variable', n_id )
if n_id ne 1 then stop
indep_var_label = temp_out[id[0]]
out_var_attribute_label = [ out_var_attribute_label, 'orig_name' ]
out_var_attribute_value = [ out_var_attribute_value, indep_var_label ]
out_var_attribute_type = [ out_var_attribute_type, 7 ]
id = where( temp_select_headers eq 'scenario_comment', n_id )
if n_id ne 1 then stop
out_global_attribute_label = [ out_global_attribute_label, 'comment' ]
out_global_attribute_value = [ out_global_attribute_value, temp_out[id[0]] ]
out_global_attribute_type = [ out_global_attribute_type, 7 ]

; Remove empty attribute entries
id = where( out_global_attribute_label ne '', n_id )
if n_id eq 0 then begin
  out_global_attribute_label = ''
  out_global_attribute_value = ''
  out_global_attribute_type = 0
endif else if n_id lt n_elements( out_global_attribute_label ) then begin
  out_global_attribute_label = out_global_attribute_label[id]
  out_global_attribute_value = out_global_attribute_value[id]
  out_global_attribute_type = out_global_attribute_type[id]
endif
id = where( out_var_attribute_label ne '', n_id )
if n_id eq 0 then begin
  out_var_attribute_label = ''
  out_var_attribute_value = ''
  out_var_attribute_type = 0
endif else if n_id lt n_elements( out_var_attribute_label ) then begin
  out_var_attribute_label = out_var_attribute_label[id]
  out_var_attribute_value = out_var_attribute_value[id]
  out_var_attribute_type = out_var_attribute_type[id]
endif

; Copy and modify output variable name
out_var_label = indep_var_label
if out_var_label eq 'ts' then out_var_label = 'tos'

; Generate file name
id = where( temp_select_headers eq 'scenario_domain', n_id )
if n_id ne 1 then stop
out_file = out_var_label + '_' + temp_out[id[0]] + '_delta_' $
    + string_substitute( scenario, '/', '_' ) + '_period' + indep_period[0] $
    + '-' + indep_period[1] + '_' + indep_period[0] + '-' + indep_period[1] $
    + '.nc'

;***********************************************************************
; Define the simulations to use

; Determine the number of source-scenarios to use
temp_select_headers = [ 'independent_algebra' ]
temp_select_values = [ 'scenario_name=' + scenario ]
markup_read, file_scenario, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
attrib_formula_label = temp_out[0]
if attrib_formula_label eq '' then stop
indep_label = strtrim( strsplit( attrib_formula_label, '-+', extract=1, $
    count=n_indep ), 2 )
; Determine the formula algebra
attrib_formula = ''
pos = [ strsplit( attrib_formula_label, '-+', count=n_pos ), $
    strlen( attrib_formula_label ) + 1 ]
for i_pos = 0, n_pos - 1 do begin
  temp = strmid( attrib_formula_label, pos[i_pos], $
      pos[i_pos+1] - pos[i_pos] - 1 )
  temp = strtrim( temp, 2 )
  id = where( indep_label eq temp, n_id )
  if n_id ne 1 then stop
  if pos[i_pos] eq 0 then begin
    ; Assume a lack of a sign means plus
    attrib_formula = attrib_formula + '+' + str( id[0] )
  endif else begin
    temp = strmid( attrib_formula_label, pos[i_pos] - 1, 1 )
    if temp eq '-' then begin
      attrib_formula = attrib_formula + '-' + str( id[0] )
    endif else if temp eq '+' then begin
      attrib_formula = attrib_formula + '+' + str( id[0] )
    endif else begin
      ; Assume a lack of a sign means plus
      attrib_formula = attrib_formula + '+' + str( id[0] )
    endelse
  endelse
endfor
; Load scenario definition
temp_select_headers = [ 'independent_' + indep_label + '_file' ]
markup_read, file_scenario, comment_char=';', $
    select_values=temp_select_values, select_headers=temp_select_headers, $
    settings=temp_out
for i_indep = 0, n_indep - 1 do begin
  id = where( temp_select_headers $
      eq 'independent_' + indep_label[i_indep] + '_file', n_id )
  if n_id ne 1 then stop
  temp_list = strtrim( strsplit( temp_out[id[0]], '&', extract=1, $
      count=n_temp_list ), 2 )
  if max( temp_list eq '' ) eq 1 then stop
  if i_indep eq 0 then begin
    n_indep_file = n_temp_list
    indep_file = strarr( n_indep_file, n_indep )
  endif
  if n_temp_list gt n_indep_file then begin
    indep_file = [ indep_file, strarr( n_temp_list - n_indep ) ]
    n_indep_file = n_temp_list
  endif
  indep_file[0:n_temp_list-1,i_indep] = temp_list
endfor
; Add directories to file lists
if keyword_set( in_dir ) then begin
  for i_indep = 0, n_indep - 1 do begin
    for i_file = 0, n_indep_file - 1 do begin
      temp_file = indep_file[i_file,i_indep]
      if temp_file ne '' then begin
        temp_file = strsplit( temp_file, ',', extract=1, count=n_temp_file )
        for i_temp_file = 0, n_temp_file - 1 do begin
          temp = file_search( in_dir + '/' + temp_file[i_temp_file], $
              count=n_temp )
          if n_temp eq 0 then stop
          if n_temp gt 1 then temp = strjoin( temp, ',' )
          temp_file[i_temp_file] = temp[0]
        endfor
        indep_file[i_file,i_indep] = strjoin( temp_file, ',' )
      endif
    endfor
  endfor
endif

; Build file list for target grid
if keyword_set( dependent_file ) then begin
  dependent_file_list = file_search( dependent_file, $
      count=n_dependent_file_list )
  if n_dependent_file_list eq 0 then stop
endif

;***********************************************************************
; Estimate attributable sea surface temperature difference

; Estimate delta
c20c_dtos_v2_make_delta, indep_var=indep_var_label, indep_file=indep_file, $
    attrib_formula=attrib_formula, $
    grid_file=grid_file, time_filter=time_filter, options=options, $
    indep_period=indep_period, out_file=out_file, out_var_label=out_var_label, $
    out_global_attribute_c20c=1, $
    out_global_attribute_label=out_global_attribute_label, $
    out_global_attribute_value=out_global_attribute_value, $
    out_global_attribute_type=out_global_attribute_type, $
    out_var_attribute_label=out_var_attribute_label, $
    out_var_attribute_value=out_var_attribute_value, $
    out_var_attribute_type=out_var_attribute_type, out_period=out_period
    ; These are not yet implemented in c20c_dtos_v2_make_delta.pro.
    ;dep_var=dep_var, dep_file=dep_file_list, 

;***********************************************************************
; The end

;stop
return
END
