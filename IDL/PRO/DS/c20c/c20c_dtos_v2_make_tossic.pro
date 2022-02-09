;+
; NAME:
;    c20c_dtos_v2_make_tossic
;
; PURPOSE:
;    This procedure generates sea surface temperature and sea ice concentration 
;    data for use in the C20C+ Detection and Attribution project.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_make_tossic, in_file=in_file, in_var_label=in_var_label, out_file=out_file, out_var_label=out_var_label
;
; INPUTS:
;    IN_FILE, IN_VAR_LABEL, OUT_FILE, OUT_VAR_LABEL
;    DELTA_APPLY_VAR_LABEL, DELTA_FACTOR, DELTA_FILE, DELTA_PERIOD, 
;      DELTA_VAR_LABEL, GRID_FILE, GRID_VAR_LABEL, IN_LAT_ATTRIBUTE_LIST, 
;      IN_LON_ATTRIBUTE_LIST, IN_PERIOD, IN_TIME_ATTRIBUTE_LIST, IN_TIME_LABEL, 
;      IN_VAR_ATTRIBUTE_LIST, OCEAN_MASK_DATA, OUT_GLOBAL_ATTRIBUTE_LIST, 
;      OUT_LAT_ATTRIBUTE_LIST, OUT_LON_ATTRIBUTE_LIST, OUT_PERIOD, 
;      OUT_TIME_ATTRIBUTE_LIST, OUT_TIME_LABEL, OUT_TIME_TYPE, 
;      OUT_VAR_ATTRIBUTE_LIST, OUT_VAR_TYPE, TOS_FREEZE, TOSSIC_FIT_METHOD, 
;      TOSSIC_FIT_PERIOD, TOSSIC_FIT_SIC_FILE, TOSSIC_FIT_TOS_FILE
;
; KEYWORD PARAMETERS:
;    CF_STANDARD:  If set then the output NetCDF file follows the CF (Climate 
;        Forecast) label and metadata standards.
;    DELTA_APPLY_VAR_LABEL:  An optional vector string containing the name of 
;        the input variable(s) to which the data from DELTA_VAR_LABEL in 
;        DELTA_FILE should be added.  The default is IN_VAR_LABEL.  This can be 
;        used if, for instance, DELTA_VAR_LABEL is 'ts' (skin temperature) and 
;        is being added to the input variable 'tos' (sea surface temperature) 
;        but not another requested input variable 'sic' (sea ice concentration).
;    DELTA_EXTEND_CYCLE:  If set, then the delta data is assumed to cover only 
;        one annual cycle, and is extended to cover the full IN_PERIOD period.
;    DELTA_FACTOR:  An optional vector number of length N_DELTA that lists the 
;        amount by which to multiply the data from each delta pattern in 
;        DELTA_FILE when adding that data to the input data.  For example, in 
;        order to subtract set this to "-1", as is done when using the 
;        Nat-Hist/CMIP5-est1 attributable warming estimate.  The default is 1.
;    DELTA_FILE:  An optional string vector of size N_DELTA containing the 
;        names of the files, including absolute or relative directory paths, 
;        containing data to be added to the/an input variable data when 
;        calculating the output.  N_DELTA is the number of different patterns 
;        to add in a linear combination (with amplitudes/signs defined in 
;        DELTA_FACTOR).  For each pattern, multiple files should be listed in 
;        comma-delimited form, with the inclusion of the "*" wildcard 
;        permitted.  For example, for the Nat-Hist/CMIP5-est1 sea surface 
;        estimate, DELTA_FILE should be a 1-element vector containing the name 
;        of the single file containing the Nat-Hist/CMIP5-est1 estimate of 
;        attributable warming.
;    DELTA_PERIOD:  An optional 2*N_DELTA string array containing the time 
;        period for which to restrict data read from the DELTA_FILE files.  The 
;        format should be 'yyyymmdd', so for instance ['20010101','20101231'] 
;        would input data covering the 1 January 2001 through 31 December 2010 
;        period.  The default is IN_PERIOD.
;    DELTA_VAR_LABEL:  An optional vector string of length N_DELTA containing 
;        the label of the variables for each delta pattern in DELTA_FILE which 
;        should be added to the input variable.  If a single entry is input 
;        then it is assumed to apply to all N_DELTA patterns.  The default is 
;        IN_VAR_LABEL[0].
;    DOUBLE:  If set then the Sheng and Zwiers (1998) calculation in 
;        shengzwiers.pro is performed in double numerical precision.
;    GRID_FILE:  An optional string scalar containing the name of a file, 
;        including the absolute or relative directory path, containing 
;        longitude and latitude vectors to which to interpolate the input data 
;        (and delta data) as it is being read from the file.
;    GRID_VAR_LABEL:  An optional string scalar containing the name of the 
;        data variable to read from GRID_FILE.  This may be needed if there are 
;        multiple longitude and latitude vectors in GRID_FILE, in order to 
;        identify those vectors corresponding to GRID_VAR_LABEL.  Also, if the 
;        value of 'sftlf' is input and OCEAN_MASK_DATA is not input, then the 
;        sftlf data (or rather one minus the data) is used to generate a 
;        land-ocean mask used by the sea ice concentration calculation code.
;    IN_FILE:  A required string array containing the names of the files, 
;        including absolute or relative directory paths, containing the input 
;        data to be used.  Of size N_IN_FILE*N_IN_VAR, with entry I in the 
;        N_IN_VAR dimension corresponding to entry I in the variables list 
;        IN_VAR_LABEL.  Entries can include the '*' wildcard and lists within 
;        '{' and '}' brackets, as for instance in 'tos_199{1,2,3}*.nc', and can 
;        include comma-delimited lists (thus "N_IN_FILE" above may not be the 
;        actual number of files being considered).
;    IN_LAT_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to read from or use when reading the input latitude 
;        dimension variable.  For use in reading the input (e.g. a missing 
;        'units' attribute in the input file) the format should be 
;        'name=value' to use attribute "name" with value "value".  To read a 
;        value for output to in OUT_FILE according to OUT_LAT_ATTRIBUTE_LIST, 
;        the format can be 'name=value' or just 'name', such that the value 
;        for attribute "name" is read from the file.
;    IN_LON_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to read from or use when reading the input longitude 
;        dimension variable.  For use in reading the input (e.g. a missing 
;        'units' attribute in the input file) the format should be 
;        'name=value' to use attribute "name" with value "value".  To read a 
;        value for output to in OUT_FILE according to OUT_LON_ATTRIBUTE_LIST, 
;        the format can be 'name=value' or just 'name', such that the value 
;        for attribute "name" is read from the file.
;    IN_PERIOD:  An optional two-element string vector containing the time 
;        period for which to restrict data read from the IN_FILE files.  The 
;        format should be 'yyyymmdd', so for instance ['20010101','20101231'] 
;        would input data covering the 1 January 2001 through 31 December 2010 
;        period.
;    IN_TIME_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to read from or use when reading the input time dimension 
;        variable.  For use in reading the input (e.g. a missing 'units' 
;        attribute in the input file) the format should be 'name=value' to use 
;        attribute "name" with value "value".  To read a value for output to 
;        in OUT_FILE according to OUT_LAT_ATTRIBUTE_LIST, the format can be 
;        'name=value' or just 'name', such that the value for attribute 
;        "name" is read from the file.
;    IN_TIME_LABEL:  An optional string scalar containing the label for the 
;        time dimension variable in the input files listed in IN_FILE.  The 
;        default is 'time'.
;    IN_VAR_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to read from or use when reading the input data variables.  
;        Of size N_IN_VAR_ATTRIBUTE_LIST*N_IN_VAR.  For use in reading the 
;        input (e.g. a missing 'units' attribute in the input file) the format 
;        should be 'name=value' to use attribute "name" with value "value".  
;        To read a value for output to in OUT_FILE according to 
;        OUT_VAR_ATTRIBUTE_LIST, the format can be 'name=value' or just 
;        'name', such that the value for attribute "name" is read from the 
;        file.
;    IN_VAR_LABEL:  A required string vector containing the labels of the input 
;        variables to read from IN_FILE and which to perform operations on in 
;        order to produce the output variables.  This may work with other 
;        variables, but thus far it has only been tested with 'sic' (sea ice 
;        concentration) and 'tos' (sea surface temperature).  Of length 
;        N_IN_VAR.
;    OCEAN_MASK_DATA:  An optional floating-point array defining the fraction 
;        the area of each grid cell that is ocean.  Of size N_IN_LON*N_IN_LAT, 
;        where these are the same sizes of the spatial dimensions of the 
;        variables read from the IN_FILE files.  Values should be 1 over ocean, 
;        0 over land, and between 0 and 1 for partially-oceanic cells.
;    OUT_FILE:  A required string scalar containing the name of the file, 
;        including the absolute or relative directory path, of the output 
;        file.  If the OUT_GLOBAL_ATTRIBUTE_C20C option is set then OUT_FILE 
;        must conform to the C20C+ D&A conventions.
;    OUT_GLOBAL_ATTRIBUTE_C20C:  If set then the standard C20C+ D&A project 
;        global attribute metadata is added to the output file.  OUT_FILE must 
;        conform to C20C+ D&A conventions in order for this option to be 
;        implemented correctly.
;    OUT_GLOBAL_ATTRIBUTE_LIST:  An optional string vector containing the list 
;        of global attributes to include in the output file.  The format should 
;        be 'name=value=type' to include attribute "name" with value "value" 
;        of data type "type" (e.g. "7" for string), or 'name=value' to include 
;        attribute "name" with value "value" of default type string, or 'name' 
;        tp include attribute "name" with an automatically determined value.  
;        Automatic determination may be from values in the input file or 
;        through other methods.  Note that 'creation_date' and 'history' are 
;        automatically included in this output list.  Note that if 
;        OUT_GLOBAL_ATTRIBUTE_C20C is set then the following attributes are 
;        automatically included with default C20C+ D&A project values but can 
;        be overridden with this keyword input:  'experiment', 
;        'experiment_family', 'frequency', 'license', 'model_id', 
;        'parent_experiment', 'parent_experiment_family', 'parent_run_id', 
;        'parent_subexperiment', 'project_id', 'realm', 'references', 'run_id', 
;        and 'subexperiment';  a default 'comment' attribute may be added too.
;    OUT_LAT_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to include for the output latitude dimension variable.  The 
;        format should be 'name=value=type' to include attribute "name" with 
;        value "value" of data type "type" (e.g. "7" for string), 'name=value' 
;        to include attribute "name" with value "value" of default type string, 
;        or 'name' to include attribute "name" with an automatically determined 
;        value.  Automatic determination may be from values in the input file 
;        (if specified in IN_LAT_ATTRIBUTE_LIST), from CF-standards, or through 
;        other methods.  Note that if CF_STANDARD is set then the following 
;        attributes are automatically included:  'axis', 'long_name', 
;        'standard_name', and 'units'.
;    OUT_LON_ATTRIBUTE_LIST:  An optional string vector containing the list of 
;        attributes to include for the output longitude dimension variable.  
;        The format should be 'name=value=type' to include attribute "name" 
;        with value "value" of data type "type" (e.g. "7" for string), 
;        'name=value' to include attribute "name" with value "value" of default 
;        type string, or 'name' to include attribute "name" with an 
;        automatically determined value.  Automatic determination may be from 
;        values in the input file (if specified in IN_LON_ATTRIBUTE_LIST), from 
;        CF-standards, or through other methods.  Note that if CF_STANDARD is 
;        set then the following attributes are automatically included:  'axis', 
;        'long_name', 'standard_name', and 'units'.
;    OUT_PERIOD:  An optional two-element string vector containing the time 
;        period for which to restrict data output to OUT_FILE file.  The format 
;        should be the same as for IN_PERIOD.  When producing 'bcs*' variables, 
;        the variability adjustments become stable for the first and last 
;        months when OUT_PERIOD exceeds IN_PERIOD by 2-3 months at either end.  
;        If IN_PERIOD is set then the default is IN_PERIOD.
;    OUT_TIME_ATTRIBUTE_LIST:  An optional string array containing the list of 
;        attributes to include for the output time dimension variables.  Of 
;        size N_OUT_TIME_ATTRIBUTE*N_OUT_TIME_VAR.  The format should be 
;        'name=value=type' to include attribute "name" with value "value" of 
;        data type "type" (e.g. "7" for string), 'name=value' to include 
;        attribute "name" with value "value" of default type string, or 'name' 
;        to include attribute "name" with an automatically determined value, 
;        or blank ('').  Automatic determination may be from values in the 
;        input file (if specified in IN_TIME_ATTRIBUTE_LIST), from 
;        CF-standards, or through other methods ("units").  The following 
;        attributes are automatically included:  'calendar' and 'units'.  Note 
;        that if CF_STANDARD is set then the following attributes are 
;        automatically included:  'axis', 'long_name', 'standard_name', and 
;        'units'.
;    OUT_TIME_LABEL:  An optional string vector containing the list of time 
;        dimension variables to output.  Acceptable values are 'time' (the 
;        regular 'days since...' time dimension variable), 'date' (containing 
;        the 'YYYYMMDD' date values), and 'datesec' (containing the number of 
;        seconds since the beginning of the date listed in 'date').  Of length 
;        The N_OUT_TIME_VAR.  The default is 'time'.
;    OUT_TIME_TYPE:  An optional integer vector containing the variable type 
;        identifiers for the time dimension variables listed in 
;        OUT_TIME_LABEL.  Of length N_OUT_TIME_VAR.  The default is to adopt 
;        the same type as for the IN_TIME_LABEL read from file, but that 
;        requires that the variable be included in IN_TIME_LABEL.  See 
;        var_type.pro for the variable type identifiers.
;    OUT_VAR_ATTRIBUTE_LIST:  An optional string array containing the list of 
;        attributes to include for the output data variables.  Of size 
;        N_OUT_VAR_ATTRIBUTE*N_OUT_VAR.  The format should be 
;        'name=value=type' to include attribute "name" with value "value" 
;        of data type "type" (e.g. "7" for string), 'name=value' to include 
;        attribute "name" with value "value" of default type string, or 'name' 
;        to include attribute "name" with an automatically determined value 
;        or blank ('').  Automatic determination may be from values in the 
;        input file (if specified in IN_VAR_ATTRIBUTE_LIST), from 
;        CF-standards, or through other methods.  Note that if CF_STANDARD is 
;        set then the following attributes are automatically included:  
;        'long_name', 'standard_name', and 'units'.  For a 'Xbcs' variable 
;        generated from an input 'X' variable the 'orig_name' attribute is 
;        automatically included.  Entering 'label=value' will not add a 
;        "label" attribute but rather will rename the variable to "value" in 
;        the output NetCDF file.
;    OUT_VAR_LABEL:  A required string vector containing the list of data 
;        variables to output to OUT_FILE.  Dimension variables should not be 
;        included in this list.  The code has been tested for 'sicbcs', 
;        'tosbcs', 'sic', and 'tos', but may work for other variables too.  A 
;        'Xbcs' variable is a version of the 'X' variable through application 
;        of the Sheng and Zwiers (1998) variability adjustment such that the 
;        monthly means of variable 'X' are recovered upon monthly averaging of 
;        the daily values that are linearly interpolated from the monthly 
;        'Xbcs' values.  If 'Xbcs' variable is included in OUT_VAR_LABEL but 
;        not in IN_VAR_LABEL (but variable 'X' is in IN_VAR_LABEL), then the 
;        'Xbcs' data is automatically calculated from variable 'X'.  Of length 
;        N_OUT_VAR.
;    OUT_VAR_TYPE:  An optional integer vector containing the variable type 
;        identifiers for the variables listed in OUT_VAR_LABEL.  Of length 
;        N_OUT_VAR.  The default is to adopt the same type as for the 
;        IN_VAR_LABEL read from file, but that requires that the variable be 
;        included in IN_VAR_LABEL.  See var_type.pro for the variable type 
;        identifiers.
;    SHENGZWIERS_EXTEND_CYCLE:  Passed as the EXTEND_CYCLE option to 
;        shengzwiers.pro when performing the Sheng and Zwiers (1998) 
;        calculation.
;    TOS_FREEZE:  An optional floating-point scalar defining the freezing point 
;        in Kelvin.  The default is 271.35K.
;    TOSSIC_FIT_METHOD:  The FIT_METHOD input required for use of the 
;        c20c_dtos_v2_adjust_sic.pro procedure for adjustment of sea ice 
;        concentration variables for consistency with adjusted sea surface 
;        temperature variables.
;    TOSSIC_FIT_PERIOD:  The FIT_PERIOD input optional for use of the 
;        c20c_dtos_v2_adjust_sic.pro procedure for adjustment of sea ice 
;        concentration variables for consistency with adjusted sea surface 
;        temperature variables.
;    TOSSIC_FIT_SIC_FILE:  The FIT_SIC_FILE input optional for use of the 
;        c20c_dtos_v2_adjust_sic.pro procedure for adjustment of sea ice 
;        concentration variables for consistency with adjusted sea surface 
;        temperature variables.
;    TOSSIC_FIT_TOS_FILE:  The FIT_TOS_FILE input optional for use of the 
;        c20c_dtos_v2_adjust_sic.pro procedure for adjustment of sea ice 
;        concentration variables for consistency with adjusted sea surface 
;        temperature variables.
;    V1: If set then the routine reproduces the output of the v1-0 code.  The 
;        difference is that the v1-0 set ice coverage to full as soon as the 
;        freezing point was passed, whereas the v2-0 code ignores the freezing 
;        point while shifting along the algorithm.  There is also a small 
;        difference in the selection of bins used to calculate the SST-SIC 
;        function, and v1-0 uses a time varying land-sea mask to for the 
;        calculation of the function whereas v2 uses a time-invariant mask.
;
; OUTPUTS:
;    -
;    (OUT_FILE)
;
; USES:
;    add_dim.pro
;    c20c_dtos_v2_adjust_sic.pro
;    c20c_dtos_v2_unnan.pro
;    convert_time_format.pro
;    isin.pro
;    month_day.pro
;    month_name.pro
;    netcdf_read.pro
;    netcdf_read_geo.pro
;    netcdf_read_geo_multitime.pro
;    netcdf_write.pro
;    mask_lonlattime.pro
;    shengzwiers.pro
;    str.pro
;    string_from_vector.pro
;    var_type.pro
;
; PROCEDURE:
;    The procedure:  reads input sea surface temperature and/or sea ice 
;    concentration data from file;  performs certain modifications such as 
;    spatial interpolation, addition of a delta field, and the Sheng and Zwiers 
;    (1998) adjustment;  and writes the result to a NetCDF file.
;
; EXAMPLE:
;    See c20c_dtos_v2_make_tossic_driver.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-06-18, as 
;        c20c_sheng_zwiers.pro.
;    Modified:  DAS, 2017-12-24 (Branched from c20c_sheng_zwiers.pro;  $
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2018-03-02 (Modified method of determining variable type 
;        of attributes;  added global attribute datatype information;  fixed 
;        issue with proper copying of variable attribute information;  modified 
;        interaction with .xml library files;  fixed error in passing input 
;        global attributes;  fixed metadata for SIC adjustment)
;    Modified:  DAS, 2018-05-06 (Added ability to extract land-sea mask from 
;        sftlf file with various possible units)
;    Modified:  DAS, 2018-08-29 (Added missing correction of sub-freezing 
;        adjusted sea surface temperatures;  Switched from usage of 
;        c20c_unnan.pro to c20c_dtos_v2_unnan.pro;  Fixed issues with 
;        implementation of multiple time dimension variables;  Set to default 
;        of no delta-consistent adjustment to sea ice concentration of 
;        TOSSIC_FIT_METHOD is not input)
;    Modified:  DAS, 2018-09-11 (Added capability for usage of multiple delta 
;        patterns in a linear combination;  Implemented DELTA_PERIOD)
;    Modified:  DAS, 2018-09-18 (Added informative statements on some stops;  
;        ensured input of OUT_PERIOD and OUT_VAR_ATTIBUTE_LIST is optional)
;    Modified:  DAS, 2018-10-11 (Added DELTA_VAR_ORIG_LABEL input;  Fixed error 
;        in call to c20c_unnan.pro for delta pattern;  Fixed misspelling of 
;        n_delta_time)
;    Modified:  DAS, 2018-11-12 (Added the DELTA_EXTEND_CYCLE keyword option;  
;        Added capability to modify time and variable data in response to 
;        prescribed calendar and units)
;    Modified:  DAS, 2018-11-28 (Corrected error in triggering data conversion 
;        when input and output units differ, shifted to later in code)
;-

;***********************************************************************

PRO C20C_DTOS_V2_MAKE_TOSSIC, $
    IN_FILE=in_file, $
    IN_VAR_LABEL=in_var_label, $
    IN_TIME_LABEL=in_time_label, $
    IN_VAR_ATTRIBUTE_LIST=in_var_attribute_list, $
      IN_LON_ATTRIBUTE_LIST=in_lon_attribute_list, $
      IN_LAT_ATTRIBUTE_LIST=in_lat_attribute_list, $
      IN_TIME_ATTRIBUTE_LIST=in_time_attribute_list, $
    IN_PERIOD=in_period, $
    OUT_FILE=out_file, $
    OUT_VAR_LABEL=out_var_label, OUT_VAR_TYPE=out_var_type, $
    OUT_TIME_LABEL=out_time_label, OUT_TIME_TYPE=out_time_type, $
    OUT_VAR_ATTRIBUTE_LIST=out_var_attribute_list, $
      OUT_LON_ATTRIBUTE_LIST=out_lon_attribute_list, $
      OUT_LAT_ATTRIBUTE_LIST=out_lat_attribute_list, $
      OUT_TIME_ATTRIBUTE_LIST=out_time_attribute_list, $
    OUT_PERIOD=out_period, $          
    OUT_GLOBAL_ATTRIBUTE_LIST=out_global_attribute_list, $
    DELTA_FILE=delta_file, DELTA_FACTOR=delta_factor, $
      DELTA_VAR_LABEL=delta_var_label, $
      DELTA_VAR_ORIG_LABEL=delta_var_orig_label, $
      DELTA_APPLY_VAR_LABEL=delta_apply_var_label, $
    DELTA_PERIOD=delta_period, DELTA_EXTEND_CYCLE=delta_extend_cycle_opt, $
    GRID_FILE=grid_file, GRID_VAR_LABEL=grid_var_label, $
    OCEAN_MASK_DATA=ocean_mask_data, $
    TOSSIC_FIT_METHOD=tossic_fit_method, TOSSIC_FIT_PERIOD=tossic_fit_period, $
      TOSSIC_FIT_SIC_FILE=tossic_fit_sic_file, $
      TOSSIC_FIT_TOS_FILE=tossic_fit_tos_file, $
    SHENGZWIERS_EXTEND_CYCLE=shengzwiers_extend_cycle, $
    TOS_FREEZE=tos_freeze, $
    CF_STANDARD=cf_standard_opt, $
    DOUBLE=double_opt, $
    OUT_GLOBAL_ATTRIBUTE_C20C=out_global_attribute_c20c_opt, $
    V1=v1_opt

;***********************************************************************
; Constants

; Missing data flag
nan = !values.f_nan
; Hard carriage return character
hard_return = string( 10B )
; The default freezing point
if n_elements( tos_freeze ) eq 0 then tos_freeze = 273.15 - 1.8

; Ensure at least one input variable has been requested
n_in_var = n_elements( in_var_label )
if n_in_var eq 0 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'No input variables defined (IN_VAR_LABEL).'
endif

; The input period
if keyword_set( in_period ) then begin
  if n_elements( in_period[*,0] ) ne 2 then stop
  if n_elements( in_period[0,*] ) ne n_in_var then begin
    if n_elements( in_period[0,*] ) eq 1 then begin
      in_period = add_dim( in_period, 1, n_in_var )
    endif else begin
      stop
    endelse
  endif
endif

; The default output time variable label
if not( keyword_set( out_time_label ) ) then out_time_label = 'time'

; The default time label
if keyword_set( out_time_label ) or keyword_set( out_time_attribute_list ) $
    and not( keyword_set( in_time_label ) ) then begin
  in_time_label = 'time'
endif

; Ensure at least one output variable has been requested
n_out_var = n_elements( out_var_label )
if n_out_var eq 0 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'No output variables defined (OUT_VAR_LABEL).'
endif
; Create a list of variables to be potentially prepared for output.
; This needs to be done to ensure that non-bcs variables are included in the 
; processing (even if they are not output) for bcs variables
act_var_label = out_var_label
for i_var = 0, n_out_var - 1 do begin
  if strpos( out_var_label[i_var], 'bcs' ) $
      eq strlen( out_var_label[i_var] ) - 3 then begin
    if ( max( in_var_label eq out_var_label[i_var] ) eq 0 ) $
        and ( max( out_var_label eq out_var_label[i_var] ) eq 0 ) then begin
      act_var_label = [ act_var_label, $
          strmid( out_var_label[i_var], 0, $
          strlen( out_var_label[i_var] ) - 3 ) ]
    endif
  endif
endfor
n_act_var = n_elements( act_var_label )
; Ensure the output variables are in this list
id = where( isin( act_Var_label, out_var_label ) eq 0, n_id )
if n_id gt 0 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'Some output variables in OUT_VAR_LABEL not supported (' $
      + strjoin( out_var_label[id], ', ' ) + ')'
endif
; Count the number of time-related variables to output
n_out_time_var = n_elements( out_time_label )
if ( n_elements( out_time_attribute_list ) gt 0 ) $
    and ( n_out_time_var eq 0 ) then begin
  out_time_label = 'time'
  n_out_time_var = 1
endif

; Ensure an input file or list of input files has been provided
if not( keyword_set( in_file ) ) then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  No input files defined (IN_FILE).'
endif
n_in_file = n_elements( in_file[*,0] )
if n_elements( in_file[0,*] ) ne n_in_var then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'The number of variables in the second dimension of IN_FILE (' $
      + str( n_elements( in_file[0,*] ) ) $
      + ') does not equal the number of input variables in IN_VAR_LABEL (' $
      + str( n_in_var ) + ').'
endif

; Ensure an output file has been provided
n_out_file = n_elements( out_file )
if n_out_file ne 1 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'No output file defined (OUT_FILE).'
endif

; Count the number of attributes for input data variables
n_in_var_attribute = n_elements( in_var_attribute_list )
if n_in_var_attribute ne 0 then begin
  n_in_var_attribute = n_elements( in_var_attribute_list[*,0] )
  if n_elements( in_var_attribute_list[0,*] ) ne n_in_var then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'The number of variables in the second dimension of ' $
        + 'IN_VAR_ATTRIBUTE_LIST (' $
        + str( n_elements( in_var_attribute_list[0,*] ) ) $
        + ') does not equal then number of input variables in IN_VAR_LABEL (' $
        + str( n_in_var ) + ').'
  endif
endif
; Count the number of attributes for input dimensions variables
n_in_lon_attribute = n_elements( in_lon_attribute_list )
n_in_lat_attribute = n_elements( in_lat_attribute_list )
n_in_time_attribute = n_elements( in_time_attribute_list )

; Count the number of attributes for output data variables
n_out_var_attribute = n_elements( out_var_attribute_list )
if n_out_var_attribute ne 0 then begin
  n_out_var_attribute = n_elements( out_var_attribute_list[*,0] )
  if n_elements( out_var_attribute_list[0,*] ) ne n_out_var then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'The number of variables in the second dimension of ' $
        + 'OUT_VAR_ATTRIBUTE_LIST (' $
        + str( n_elements( out_var_attribute_list[0,*] ) ) $
        + ') does not equal then number of output variables in OUT_VAR_LABEL ' $
        + '(' + str( n_in_var ) + ').'
  endif
endif
; Count the number of attributes for output dimensions variables
n_out_lon_attribute = n_elements( out_lon_attribute_list )
n_out_lat_attribute = n_elements( out_lat_attribute_list )
if n_out_time_var eq 0 then begin
  n_out_time_attribute = 0
endif else begin
  n_out_time_attribute = n_elements( out_time_attribute_list ) / n_out_time_var
endelse

; Ensure consistent output variable type vectors
if n_elements( out_var_type ) eq 0 then begin
  out_var_type = -1 + intarr( n_out_var )
endif else begin
  if n_elements( out_var_type ) ne n_out_var then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'The number of entries in OUT_VAR_TYPE should equal the number of ' $
        + 'output variables in OUT_VAR_LABEL (' + str( n_out_var ) + ').'
  endif
endelse
if n_out_time_var gt 0 then begin
  if n_elements( out_time_type ) eq 0 then begin
    out_time_type = -1 + intarr( n_out_time_var )
  endif else begin
    if n_elements( out_time_type ) ne n_out_time_var then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The number of entries in OUT_TIME_TYPE should equal the number ' $
        + 'of output time variables in OUT_TIME_LABEL (' $
        + str( n_out_time_var ) + ').'
  endif
  endelse
endif

; The number of global attributes requested
n_out_global_attribute = n_elements( out_global_attribute_list )

; The number of delta patterns to add
n_delta = n_elements( delta_file )
; If we are working with delta patterns
if n_delta gt 0 then begin
  ; The default factor by which to multiply the difference field before
  ; subtracting it (e.g. use -1 for addition)
  n_temp = n_elements( delta_factor )
  if n_temp eq 0 then begin
    delta_factor = 1 + intarr( n_delta )
  endif else if n_temp ne n_delta then begin
    if n_temp eq 1 then begin
      if n_delta ne 1 then delta_factor = delta_factor[0] + intarr( n_delta )
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The number of multiplication factors for the delta patterns as ' $
          + 'defined in DELTA_FACTOR (' + str( n_temp ) $
          + ') does not equal the number of delta patterns defined in ' $
          + 'DELTA_FILE (' + str( n_delta ) + ').'
    endelse
  endif
  ; The default delta variable
  n_temp = n_elements( delta_var_label )
  if n_temp eq 0 then begin
    delta_var_label = in_var_label[0] + strarr( n_delta )
  endif else if n_temp ne n_delta then begin
    if n_temp eq 1 then begin
      if n_delta ne 1 then begin
        delta_var_label = delta_var_label[0] + strarr( n_delta )
      endif
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The number of input variables for the delta patterns as ' $
          + 'defined in DELTA_VAR_LABEL (' + str( n_temp ) $
          + ') does not equal the number of delta patterns defined in ' $
          + 'DELTA_FILE (' + str( n_delta ) + ').'
    endelse
  endif
  ; The default variable to interpret the delta difference field as applying to
  n_temp = n_elements( delta_apply_var_label )
  if n_temp eq 0 then begin
    delta_apply_var_label = in_var_label
  endif else if n_temp ne n_in_var then begin
    if n_temp eq 1 then begin
      if n_in_var ne 1 then begin
        delta_apply_var_label = delta_apply_var_label[0] + strarr( n_in_var )
      endif
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The number of application variables for the delta patterns as ' $
          + 'defined in DELTA_APPLY_VAR_LABEL (' + str( n_temp ) $
          + ') does not equal the number of delta patterns defined in ' $
          + 'DELTA_FILE (' + str( n_delta ) + ').'
    endelse
  endif
  ; The default period for which to load data for the delta field
  n_temp = n_elements( delta_period ) / 2
  if n_temp eq 0 then begin
    ; The default is IN_PERIOD
    if keyword_set( in_period ) then begin
      delta_period = add_dim( in_period[*,0], 1, n_delta )
    endif
  endif else if n_temp ne n_delta then begin
    if n_temp eq 1 then begin
      delta_period = add_dim( delta_period, 1, n_delta )
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The number of input periods for the delta patterns as ' $
          + 'defined in DELTA_PERIOD (' + str( n_temp ) $
          + ') does not equal the number of delta patterns defined in ' $
          + 'DELTA_FILE (' + str( n_delta ) + ').'
    endelse
  endif
endif

; The default output period
if not( keyword_set( out_period ) ) and keyword_set( in_period ) then begin
  out_period = in_period[*,0]
endif
; Ensure the output period is of full yyyymmdd format
if keyword_set( out_period ) then begin
  if strlen( out_period[0] ) ne 8 then begin
    ; If it has six numbers then assume yyyymm format
    if strlen( out_period[0] ) eq 6 then begin
      out_period[0] = out_period[0] + '01'
    ; Otherwise the format is unrecognisable
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Unrecognised format for OUT_PERIOD[0] ' $
          + '(should be yyyymmdd or yyyymm).'
    endelse
  endif
  if strlen( out_period[1] ) ne 8 then begin
    ; If it has six numbers then assume yyyymm format
    if strlen( out_period[1] ) eq 6 then begin
      temp_year = fix( strmid( out_period[1], 0, 4 ) )
      if ( temp_year mod 4 eq 0 ) $
          and ( ( temp_year mod 400 eq 0 ) or ( temp_year mod 100 ne 0 ) ) $
          then begin
        temp_leap = 1
      endif else begin
        temp_leap = 0
      endelse
      temp = month_day( fix( strmid( out_period[1], 4, 2 ) ) - 1, $
          leap=temp_leap )
      out_period[1] = out_period[1] $
          + str( temp[1] - temp[0] + 1, length=2, filler='0' )
    ; Otherwise the format is unrecognisable
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Unrecognised format for OUT_PERIOD[1] ' $
          + '(should be yyyymmdd or yyyymm).'
    endelse
  endif
endif

; The default fit period
if not( keyword_set( fit_period ) ) then fit_period = [ 2001., 2011. ]

; Option to ensure CF-standard metadata is included
cf_standard_opt = keyword_set( cf_standard_opt )

;***********************************************************************
; Parse variable and attribute names

; Parse list of input variable attributes
if n_in_var_attribute gt 0 then begin
  ; Initialise arrays containing attribute labels and attribute values
  in_var_attribute_label = strarr( n_in_var_attribute, n_in_var )
  in_var_attribute_value = strarr( n_in_var_attribute, n_in_var )
  in_var_attribute_type = 7 + intarr( n_in_var_attribute, n_in_var )
  ; Iterate through variables
  for i_var = 0, n_in_var - 1 do begin
    ; Extract attribute labels and, if provided, values
    for i_attribute = 0, n_in_var_attribute - 1 do begin
      if in_var_attribute_list[i_attribute,i_var] ne '' then begin
        temp_attribute = strsplit( in_var_attribute_list[i_attribute,i_var], $
            '=', extract=1, count=n_temp_attribute )
        if n_temp_attribute gt 3 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
              + 'There should only be three components of ' $
              + 'IN_VAR_ATTRIBUTE_LIST[' + str( i_attribute ) + ',' $
              + str( i_var ) + '] separated by "=" (' $
              + in_var_attribute_list[i_attribute,i_var] + ').'
        endif
        in_var_attribute_label[i_attribute,i_var] = temp_attribute[0]
        if n_temp_attribute eq 2 then begin
          in_var_attribute_value[i_attribute,i_var] = temp_attribute[1]
          if n_temp_attribute eq 3 then begin
            in_var_attribute_type[i_attribute,i_var] = fix( temp_attribute[2] )
          endif
        endif
      endif
    endfor
  endfor
endif

; Parse list of output variable attributes
if n_out_var_attribute gt 0 then begin
  ; Initialise arrays containing attribute labels and attribute values for the 
  ; data variables
  out_var_attribute_label = strarr( n_out_var_attribute, n_out_var )
  out_var_attribute_value = strarr( n_out_var_attribute, n_out_var )
  out_var_attribute_type = 7 + intarr( n_out_var_attribute, n_out_var )
  ; Iterate through variables
  for i_var = 0, n_out_var - 1 do begin
    ; Iterate through requested attributes
    for i_attribute = 0, n_out_var_attribute - 1 do begin
      if out_var_attribute_list[i_attribute,i_var] ne '' then begin
        ; Determine output attribute label
        temp_attribute = strsplit( out_var_attribute_list[i_attribute,i_var], $
            '=', extract=1, count=n_temp_attribute, preserve_null=1 )
        if n_temp_attribute gt 3 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
              + 'There should only be three components of ' $
              + 'OUT_VAR_ATTRIBUTE_LIST[' + str( i_attribute ) + ',' $
              + str( i_var ) + '] separated by "=" (' $
              + out_var_attribute_list[i_attribute,i_var] + ').'
        endif
        out_var_attribute_label[i_attribute,i_var] = temp_attribute[0]
        ; Copy attribute value if provided
        if n_temp_attribute ge 2 then begin
          out_var_attribute_value[i_attribute,i_var] = temp_attribute[1]
          if n_temp_attribute eq 3 then begin
            out_var_attribute_type[i_attribute,i_var] = fix( temp_attribute[2] )
          endif
        endif
        ; If there is no value, then see if we can adopt an input value
        if ( out_var_attribute_value[i_attribute,i_var] eq '' ) $
            and keyword_set( in_var_attribute_value ) then begin
          ; Determine if the same variable exists in the input list
          id_var = where( in_var_label eq out_var_label[i_var], n_id_var )
          if n_id_var gt 1 then begin
            stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
                + 'There are duplicate variables defined in IN_VAR_LABEL ' $
                + 'with the label "' + out_var_label[i_var] + '".'
          endif
          if n_id_var eq 1 then begin
            ; Determine if the same attribute exists in the input list
            id_attribute = where( in_var_attribute_label[*,id_var[0]] $
                eq out_var_attribute_label[i_attribute,i_var], n_id_attribute )
            if n_id_attribute gt 1 then begin
              stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
                  + 'There are duplicate ' $
                  + out_var_attribute_label[i_attribute,i_var] $
                  + ' attributes for the variable ' + in_var_label[id_var[0]] $
                  + ' defined in IN_VAR_ATTRIBUTE_LABEL.'
            endif
            if n_id_attribute eq 1 then begin
              ; Adopt the value
              out_var_attribute_value[i_attribute,i_var] $
                  = in_var_attribute_value[id_attribute[0],id_var[0]]
              out_var_attribute_type[i_attribute,i_var] $
                  = in_var_attribute_type[id_attribute[0],id_var[0]]
            endif
          ; Otherwise determine if a non-bcs version of the variable exists in 
          ; the input list
          endif else begin
            ; But do not do this for long_name or standard_name
            if max( out_var_attribute_label[i_attribute,i_var] $
                eq [ 'long_name', 'standard_name' ] ) $
                eq 0 then begin
              id_var = where( in_var_label + 'bcs' eq out_var_label[i_var], $
                  n_id_var )
              if n_id_var gt 1 then begin
                stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
                    + 'Could not find ' + out_var_label[i_var] $
                    + ' or a bcs version of it in IN_VAR_LABEL for use in ' $
                    + 'determine long_name and standard_name attributes.'
              endif
              if n_id_var eq 1 then begin
                ; Determine if the same attribute exists in the input list
                id_attribute = where( in_var_attribute_label[*,id_var[0]] $
                    eq out_var_attribute_label[i_attribute,i_var], $
                    n_id_attribute )
                if n_id_attribute gt 1 then begin
                  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
                      + 'There are duplicate ' $
                      + out_var_attribute_label[i_attribute,i_var] $
                      + ' attributes for the variable ' $
                      + in_var_label[id_var[0]] $
                      + ' defined in IN_VAR_ATTRIBUTE_LABEL.'
                endif
                if n_id_attribute eq 1 then begin
                  ; Adopt the value
                  out_var_attribute_value[i_attribute,i_var] $
                      = in_var_attribute_value[id_attribute[0],id_var[0]]
                  out_var_attribute_type[i_attribute,i_var] $
                      = in_var_attribute_type[id_attribute[0],id_var[0]]
                endif
              endif
            endif
          endelse
        endif
      endif
    endfor
  endfor
endif

; If we have output longitude attributes defined
if n_out_lon_attribute gt 0 then begin
  ; Initialise arrays containing attribute labels and attribute values
  out_lon_attribute_label = strarr( n_out_lon_attribute )
  out_lon_attribute_value = strarr( n_out_lon_attribute )
  out_lon_attribute_type = 7 + intarr( n_out_lon_attribute )
  ; Iterate through longitude attributes
  for i_attribute = 0, n_out_lon_attribute - 1 do begin
    ; Parse this attribute
    temp_attribute = strsplit( out_lon_attribute_list[i_attribute], '=', $
        extract=1, count=n_temp_attribute )
    if n_temp_attribute gt 3 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'There should only be three components of OUT_LON_ATTRIBUTE_LIST[' $
          + str( i_attribute ) + '] separated by "=" (' $
          + out_lon_attribute_list[i_attribute] + ').'
    endif
    out_lon_attribute_label[i_attribute] = temp_attribute[0]
    if n_temp_attribute eq 2 then begin
      out_lon_attribute_value[i_attribute] = temp_attribute[1]
      if n_temp_attribute eq 3 then begin
        out_lon_attribute_type[i_attribute] = fix( temp_attribute[2] )
      endif
    endif
  endfor
endif
; Ensure CF-standard metadata is included for longitude, if requested
if cf_standard_opt eq 1 then begin
  temp_label = [ 'units', 'long_name', 'standard_name', 'axis' ]
  if n_out_lon_attribute eq 0 then begin
    out_lon_attribute_label = temp_label
    out_lon_attribute_value = strarr( n_elements( temp_label ) )
    out_lon_attribute_type = 7 + intarr( n_elements( temp_label ) )
  endif else begin
    temp = isin( out_lon_attribute_label, temp_label )
    id = where( temp eq 0, n_id )
    if n_id gt 0 then begin
      out_lon_attribute_label = [ out_lon_attribute_label, temp_label[id] ]
      out_lon_attribute_value = [ out_lon_attribute_value, strarr( n_id ) ]
      out_lon_attribute_type = [ out_lon_attribute_type, 7 + intarr( n_id ) ]
    endif
  endelse
  n_out_lon_attribute = n_elements( out_lon_attribute_label )
endif
; Adopt input values for empty output longitude attribute values, if provided
if n_in_lon_attribute gt 0 then begin
  for i_attribute = 0, n_out_lon_attribute - 1 do begin
    if out_lon_attribute_value[i_attribute] eq '' then begin
      ; Do not do this for the units attribute
      if out_lon_attribute_label[i_attribute] ne 'units' then begin
        id_in = where( in_lon_attribute_label $
            eq out_lon_attribute_label[i_attribute], n_id_in )
        if n_id_in gt 1 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
              + 'Multiple instances of ' $
              + out_lon_attribute_label[i_attribute] $
              + ' attribute in IN_LON_ATTRIBUTE_LABEL.'
        endif
        if n_id_in eq 1 then begin
          out_lon_attribute_value[i_attribute] $
              = in_lon_attribute_value[id_in[0]]
          out_lon_attribute_type[i_attribute] = in_lon_attribute_type[id_in[0]]
        endif
      endif
    endif
  endfor
endif

; If we have output latitude attributes defined
if n_out_lat_attribute gt 0 then begin
  ; Initialise arrays containing attribute labels and attribute values
  out_lat_attribute_label = strarr( n_out_lat_attribute )
  out_lat_attribute_value = strarr( n_out_lat_attribute )
  out_lat_attribute_type = 7 + intarr( n_out_lat_attribute )
  ; Iterate through latitude attributes
  for i_attribute = 0, n_out_lat_attribute - 1 do begin
    ; Parse this attribute
    temp_attribute = strsplit( out_lat_attribute_list[i_attribute], '=', $
        extract=1, count=n_temp_attribute )
    if n_temp_attribute gt 3 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'There should only be three components of OUT_LAT_ATTRIBUTE_LIST[' $
          + str( i_attribute ) + '] separated by "=" (' $
          + out_lat_attribute_list[i_attribute] + ').'
    endif
    out_lat_attribute_label[i_attribute] = temp_attribute[0]
    if n_temp_attribute eq 2 then begin
      out_lat_attribute_value[i_attribute] = temp_attribute[1]
      if n_temp_attribute eq 3 then begin
        out_lat_attribute_type[i_attribute] = fix( temp_attribute[2] )
      endif
    endif
  endfor
endif
; Ensure CF-standard metadata is included for latitude, if requested
if cf_standard_opt eq 1 then begin
  temp_label = [ 'units', 'long_name', 'standard_name', 'axis' ]
  if n_out_lat_attribute eq 0 then begin
    out_lat_attribute_label = temp_label
    out_lat_attribute_value = strarr( n_elements( temp_label ) )
    out_lat_attribute_type = 7 + intarr( n_elements( temp_label ) )
  endif else begin
    temp = isin( out_lat_attribute_label, temp_label )
    id = where( temp eq 0, n_id )
    if n_id gt 0 then begin
      out_lat_attribute_label = [ out_lat_attribute_label, temp_label[id] ]
      out_lat_attribute_value = [ out_lat_attribute_value, strarr( n_id ) ]
      out_lat_attribute_type = [ out_lat_attribute_type, 7 + intarr( n_id ) ]
    endif
  endelse
  n_out_lat_attribute = n_elements( out_lat_attribute_label )
endif
; Adopt input values for empty output latitude attribute values, if provided
if n_in_lat_attribute gt 0 then begin
  for i_attribute = 0, n_out_lat_attribute - 1 do begin
    if out_lat_attribute_value[i_attribute] eq '' then begin
      ; Do not do this for the units attribute
      if out_lat_attribute_label[i_attribute] ne 'units' then begin
        id_in = where( in_lat_attribute_label $
            eq out_lat_attribute_label[i_attribute], n_id_in )
        if n_id_in gt 1 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
              + 'Multiple instances of ' $
              + out_lat_attribute_label[i_attribute] $
              + ' attribute in IN_LAT_ATTRIBUTE_LABEL.'
        endif
        if n_id_in eq 1 then begin
          out_lat_attribute_value[i_attribute] $
              = in_lat_attribute_value[id_in[0]]
          out_lat_attribute_type[i_attribute] = in_lat_attribute_type[id_in[0]]
        endif
      endif
    endif
  endfor
endif

; If we have output time attributes defined
if n_out_time_attribute gt 0 then begin
  ; Initialise arrays containing attribute labels and attribute values
  out_time_attribute_label = strarr( n_out_time_attribute, n_out_time_var )
  out_time_attribute_value = strarr( n_out_time_attribute, n_out_time_var )
  out_time_attribute_type = 7 + intarr( n_out_time_attribute, n_out_time_var )
  ; Iterate through time variables and attributes
  for i_var = 0, n_out_time_var - 1 do begin
    for i_attribute = 0, n_out_time_attribute - 1 do begin
      ; Parse this attribute
      temp_attribute = strsplit( out_time_attribute_list[i_attribute,i_var], $
          '=', extract=1, count=n_temp_attribute )
      if n_temp_attribute gt 3 then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'There should only be three components of ' $
            + 'OUT_TIME_ATTRIBUTE_LIST[' + str( i_attribute ) + ',' $
            + str( i_var ) + '] separated by "=" (' $
            + out_time_attribute_list[i_attribute,i_var] + ').'
      endif
      out_time_attribute_label[i_attribute,i_var] = temp_attribute[0]
      if n_temp_attribute eq 2 then begin
        out_time_attribute_value[i_attribute,i_var] = temp_attribute[1]
        if n_temp_attribute eq 3 then begin
          out_time_attribute_type[i_attribute,i_var] = fix( temp_attribute[2] )
        endif
      endif
    endfor
  endfor
endif
; Ensure CF-standard metadata is included for time variables, if requested
if cf_standard_opt eq 1 then begin
  temp_label = [ 'units', 'long_name', 'standard_name', 'axis' ]
  if n_out_time_attribute eq 0 then begin
    if n_out_time_var gt 1 then begin
      out_time_attribute_label = add_dim( temp_label, 1, n_out_time_var )
    endif else begin
      out_time_attribute_label = reform( temp_label, n_elements( temp_label ), $
          1 )
    endelse
    out_time_attribute_value = strarr( n_elements( temp_label ), $
        n_out_time_var )
    out_time_attribute_type = 7 $
        + intarr( n_elements( temp_label ), n_out_time_var )
  endif else begin
    for i_var = 0, n_out_time_var - 1 do begin
      temp = isin( out_time_attribute_label[*,i_var], temp_label )
      id = where( temp eq 0, n_id )
      if n_id gt 0 then begin
        temp = strarr( n_id, n_out_time_var )
        temp[*,i_var] = temp_label[id]
        out_time_attribute_label = [ out_time_attribute_label, temp ]
        out_time_attribute_value = [ out_time_attribute_value, $
            strarr( n_id, n_out_time_var ) ]
        out_time_attribute_type = [ out_time_attribute_type, $
            7 + intarr( n_id, n_out_time_var ) ]
      endif
    endfor
  endelse
  n_out_time_attribute = n_elements( out_time_attribute_label[*,0] )
endif
; Adopt input values for empty output time attribute values, if provided
if n_in_time_attribute gt 0 then begin
  for i_var = 0, n_out_time_var - 1 do begin
    for i_attribute = 0, n_out_time_attribute - 1 do begin
      if out_time_attribute_value[i_attribute,i_var] eq '' then begin
        ; Do not do this for the units attribute
        if out_time_attribute_label[i_attribute,i_var] ne 'units' then begin
          id_in_var = where( in_time_label eq out_time_label[i_var], $
              n_id_in_var )
          if n_id_in_var gt 1 then begin
            stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
               + 'Multiple instance of ' + out_time_label[i_var] $
               + ' time variable in IN_TIME_LABEL.'
          endif
          if n_id_in_var eq 1 then begin
            id_in = where( in_time_attribute_label[*,id_in_var[0]] $
                eq out_time_attribute_label[i_attribute,i_attribute], n_id_in )
            if n_id_in gt 1 then begin
              stop, 'Error c20c_dtos_v2_make_tossic.pro:  Mulitple ' $
                  + out_time_attribute_label[i_attribute,i_attribute] $
                  + ' attributes defined in IN_TIME_ATTRIBUTE_LABEL for time ' $
                  + 'variable ' + in_time_label[id_in_var[0]] + '.'
            endif
            if n_id_in eq 1 then begin
              out_time_attribute_value[i_attribute,i_var] $
                  = in_time_attribute_value[id_in[0],id_in_var[0]]
              out_time_attribute_type[i_attribute,i_var] $
                  = in_time_attribute_type[id_in[0],id_in_var[0]]
            endif
          endif
        endif
      endif
    endfor
  endfor
endif

; Parse global attribute names and values
if n_out_global_attribute eq 0 then begin
  out_global_attribute_label = ''
  out_global_attribute_value = ''
  out_global_attribute_type = 7
endif else begin
  out_global_attribute_label = strarr( n_out_global_attribute )
  out_global_attribute_value = strarr( n_out_global_attribute )
  out_global_attribute_type = 7 + intarr( n_out_global_attribute )
  for i_attribute = 0, n_out_global_attribute - 1 do begin
    temp_attribute = strsplit( out_global_attribute_list[i_attribute], '=', $
        extract=1, count=n_temp_attribute )
    if n_temp_attribute gt 3 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'There should only be three components of ' $
          + 'OUT_GLOBAL_ATTRIBUTE_LIST[' + str( i_attribute ) $
          + ' separated by "=" (' $
          + out_global_attribute_list[i_attribute] + ').'
    endif
    out_global_attribute_label[i_attribute] = temp_attribute[0]
    if n_temp_attribute ge 2 then begin
      out_global_attribute_value[i_attribute] = temp_attribute[1]
      if n_temp_attribute eq 3 then begin
        out_global_attribute_type[i_attribute] = fix( temp_attribute[2] )
      endif
    endif
  endfor
endelse

;***********************************************************************
; Load data

; Load data describing target grid
if keyword_set( grid_file ) then begin
  ; Prepare to extract ocean mask data if possible
  if not( keyword_set( grid_var_label ) ) then begin
    if strpos( grid_file, 'sftlf' ) ge 0 then grid_var_label = 'sftlf'
  endif
  ; Load target grid data
  temp = netcdf_read_geo( grid_file, grid_var_label, lon=grid_lon_data, $
      lat=grid_lat_data, quiet=1, units_var=temp_units )
  n_grid_lon = n_elements( grid_lon_data )
  n_grid_lat = n_elements( grid_lat_data )
  ; Extract an ocean mask if possible
  if not( keyword_set( ocean_mask_data ) ) $
      and ( keyword_set( grid_var_label ) ) then begin
    ocean_mask_data = reform( temp[*,*,0,0,0] )
    if temp_units eq '%' then begin
      ocean_mask_data = ocean_mask_data / 100.
    endif else if temp_units ne 'fraction' then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Units for the mask data defined in ' + grid_file $
          + ' are not recognised (should be % or fraction).'
    endif
    if strpos( grid_file, 'sftlf' ) ge 0 then begin
      ocean_mask_data = 1. - ocean_mask_data
    endif else begin
      ocean_mask_data = 0
    endelse
  endif
endif else begin
  n_grid_lon = 0
  n_grid_lat = 0
endelse

; Iterate through input variables
for i_var = 0, n_in_var - 1 do begin
  ; Copy any time attributes required for reading
  temp_in_time_calendar = ''
  temp_in_time_units = ''
  if keyword_set( in_time_attribute_list ) then begin
    ; Copy the time calendar attribute if provided
    id = where( strpos( in_time_attribute_list, 'calendar=' ) eq 0, n_id )
    if n_id gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  Multiple instances of ' $
          + 'time calendar attribute in IN_TIME_ATTRIBUTE_LIST.'
    endif
    if n_id eq 1 then begin
      temp = strsplit( in_time_attribute_list[id[0]], '=', extract=1, $
          count=n_temp )
      if ( n_temp ne 2 ) and ( n_temp ne 3 ) then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  Incorrect format for ' $
            + 'time calendar attribute in IN_TIME_ATTRIBUTE_LIST.'
      temp_in_time_calendar = temp[1]
    endif
    ; Copy the time units attribute if provided
    id = where( strpos( in_time_attribute_list, 'units=' ) eq 0, n_id )
    if n_id gt 1 then stop
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  Multiple instances of ' $
          + 'time units attribute in IN_TIME_ATTRIBUTE_LIST.'
    endif
    if n_id eq 1 then begin
      temp = strsplit( in_time_attribute_list[id[0]], '=', extract=1, $
          count=n_temp )
      if ( n_temp ne 2 ) and ( n_temp ne 3 ) then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  Incorrect format for ' $
            + 'time units attribute in IN_TIME_ATTRIBUTE_LIST.'
      endif
      temp_in_time_units = temp[1]
    endif
  endif
  ; Copy any data variable attributes required for reading
  temp_in_var_orig_label = ''
  temp_in_var_units = ''
  if keyword_set( in_var_attribute_list ) then begin
    ; Copy the data variable original label attribute if provided
    ; (Note that the program will assume that this is the label of the variable 
    ; within the input files)
    id = where( strpos( in_var_attribute_list[*,i_var], 'orig_name=' ) eq 0, $
        n_id )
    if n_id gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Multiple instance of orig_name attribute for input variable ' $
          + in_var_label[i_var] + '.'
    endif
    if n_id eq 1 then begin
      temp = strsplit( in_var_attribute_list[id[0],i_var], '=', extract=1, $
          count=n_temp )
      if ( n_temp ne 2 ) and ( n_temp ne 3 ) then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'Incorrect format for orig_name attribute for input variable ' $
            + in_var_label[i_var] + ' in IN_VAR_ATTRIBUTE_LIST.'
      endif
      temp_in_orig_label_units = temp[1]
    endif
    ; Copy the data variable units attribute if provided
    id = where( strpos( in_var_attribute_list[*,i_var], 'units=' ) eq 0, n_id )
    if n_id gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'No units attribute found for input variable ' $
          + in_var_label[i_var] + '.'
    endif
    if n_id eq 1 then begin
      temp = strsplit( in_var_attribute_list[id[0],i_var], '=', extract=1, $
          count=n_temp )
      if ( n_temp ne 2 ) and ( n_temp ne 3 ) then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'Incorrect format for units attribute for input variable ' $
            + in_var_label[i_var] + ' in IN_VAR_ATTRIBUTE_LIST.'
      endif
      temp_in_var_units = temp[1]
    endif
  endif
  ; Copy any global attributes that might need to be retrieved
  if n_out_global_attribute gt 0 then begin
    id_global = where( out_global_attribute_value eq '', n_id_global )
    if n_id_global eq 0 then begin
      temp_global_label = ''
      temp_global_value = ''
    endif else begin
      temp_global_label = out_global_attribute_label[id_global]
      temp_global_value = out_global_attribute_value[id_global]
    endelse
  endif else begin
    n_id_global = 0
  endelse
  ; Read data from files
  id = where( in_file[*,i_var] ne '', n_id )
  if n_id eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No input files defined for input variable ' + in_var_label[i_var] $
        + '.'
  endif
  temp_file = string_from_vector( in_file[id,i_var], spacer=',', nospace=1 )
  temp_in_lon_units = ''
  temp_in_lat_units = ''
  if keyword_set( in_period ) then begin
    temp_in_period = in_period[*,i_var]
  endif else begin
    temp_in_period = ''
  endelse
  temp_in_var_data = netcdf_read_geo_multitime( temp_file, $
      in_var_label[i_var], period_time=temp_in_period, lon=temp_in_lon_data, $
      lat=temp_in_lat_data, time=temp_in_time_data, $
      units_lon=temp_in_lon_units, units_lat=temp_in_lat_units, $
      units_time=temp_in_time_units, calendar=temp_in_time_calendar, $
      units_var=temp_in_var_units, label_in_var=temp_in_var_orig_label, $
      global_label=temp_global_label, global_value=temp_global_value, quiet=1, $
      fix_time=1 )
  n_temp_in_lon = n_elements( temp_in_lon_data )
  n_temp_in_lat = n_elements( temp_in_lat_data )
  n_temp_in_time = n_elements( temp_in_time_data )
  temp_in_var_data = reform( temp_in_var_data, n_temp_in_lon, n_temp_in_lat, $
      n_temp_in_time )
  ; Remove highly negative flag values in sea surface temperature data from 
  ; HadISST1
  if in_var_label[i_var] eq 'tos' then begin
    id = where( temp_in_var_data lt -500, n_id )
    if n_id gt 1 then temp_in_var_data[id] = tos_freeze
  endif
  ; Ensure there is data
  if max( finite( temp_in_var_data ) ) eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No non-NaN data in input data from ' + temp_file + '.'
  endif
  ; Ensure prescribed and input calendars agree
  if keyword_set( out_time_attribute_list ) then begin
    id_time = where( out_time_label eq 'time', n_id_time )
    if n_id_time gt 1 then stop
    if n_id_time eq 1 then begin
      id = where( strpos( $
          out_time_attribute_list[*,id_time[0]], 'calendar=' ) eq 0, n_id )
      if n_id gt 1 then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  Multiple instances of ' $
            + 'time calendar attribute in OUT_TIME_ATTRIBUTE_LIST.'
      endif
      if n_id eq 1 then begin
        temp_out_time_calendar = strsplit( $
            out_time_attribute_list[id[0],id_time[0]], '=', extract=1 )
        temp_out_time_calendar = temp_out_time_calendar[1]
        if temp_out_time_calendar ne temp_in_time_calendar then begin
          temp_in_time_data = convert_time_format( temp_in_time_data, $
              temp_in_time_units, 'yyyymmddhhmmss', $
              calendar=temp_in_time_calendar )
          temp_in_time_data = convert_time_format( temp_in_time_data, $
              'yyyymmddhhmmss', temp_in_time_units, $
              calendar=temp_out_time_calendar )
          temp_in_time_calendar = temp_out_time_calendar
        endif
      endif
    endif
  endif
  ; Convert to standard units
  if temp_in_var_units eq 'deg_C' then begin
    temp_in_var_data = temp_in_var_data + 273.15
    temp_in_var_units = 'K'
  endif else if temp_in_var_units eq '%' then begin
    temp_in_var_data = temp_in_var_data / 100.
    temp_in_var_units = 'fraction'
  endif
  ; Ensure that sea surface temperature data is within expected range
  if in_var_label[i_var] eq 'tos' then begin
    if min( temp_in_var_data, nan=1 ) lt tos_freeze - 0.0001 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'SST in input variable ' + str( i_var ) $
          + ' is less than the freezing point (' $
          + str( min( temp_in_var_data, nan=1 ) ) + ').'
    endif
    if max( temp_in_var_data, nan=1 ) gt 323.15 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'SST in input variable ' + str( i_var ) + ' is greater than 50C (' $
          + str( max( temp_in_var_data, nan=1 ) ) + ').'
    endif
  endif
  ; Ensure that sea ice concentration data is within expected range
  if in_var_label[i_var] eq 'sic' then begin
    if min( temp_in_var_data ) lt -0.01 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'SIC in input variable ' + str( i_var ) + ' is less than zero (' $
          + str( min( temp_in_var_data, nan=1 ) ) + ').'
    endif
    if max( temp_in_var_data ) gt 1.01 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'SIC in input variable ' + str( i_var ) $
          + ' is greater than 100% (' $
          + str( 100 * max( temp_in_var_data, nan=1 ) ) + ').'
    endif
  endif
  ; Perform interpolation
  if n_grid_lon + n_grid_lat gt 0 then begin
    mask_lonlattime, temp_in_var_data, lon=temp_in_lon_data, $
        lat=temp_in_lat_data, mask_lon=grid_lon_data, mask_lat=grid_lat_data
    n_temp_in_lon = n_elements( temp_in_lon_data )
    n_temp_in_lat = n_elements( temp_in_lat_data )
  endif
  ; If this is the first variable to load
  if i_var eq 0 then begin
    ; Copy the dimension variables
    in_lon_data = temp_in_lon_data
    n_in_lon = n_elements( in_lon_data )
    in_lat_data = temp_in_lat_data
    n_in_lat = n_elements( in_lat_data )
    in_time_data = temp_in_time_data
    n_in_time = n_elements( in_time_data )
    in_time_units = temp_in_time_units
    in_time_calendar = temp_in_time_calendar
    ; Initialise the data array and copy this variable's data in
    in_var_data = nan * fltarr( n_in_lon, n_in_lat, n_in_time, n_in_var )
  ; If this is not the first variable to load
  endif else begin
    ; Ensure consistency of dimensions
    if n_temp_in_lon ne n_in_lon then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The length of the longitude vector for input variable ' $
          + in_var_label[i_var] + ' (' + str( n_temp_in_lon ) $
          + ' ) differs from the length for other variables (' $
          + str( n_in_lon ) + ').'
    endif
    if max( abs( temp_in_lon_data - in_lon_data ) ) gt 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Values in the longitude vector for input variable ' $
          + in_var_label[i_var] + ' differ from values for other variables.'
    endif
    if n_temp_in_lat ne n_in_lat then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The length of the latitude vector for input variable ' $
          + in_var_label[i_var] + ' (' + str( n_temp_in_lat ) $
          + ' ) differs from the length for other variables (' $
          + str( n_in_lat ) + ').'
    endif
    if max( abs( temp_in_lat_data - in_lat_data ) ) gt 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Values in the latitude vector for input variable ' $
          + in_var_label[i_var] + ' differ from values for other variables.'
    endif
    if n_temp_in_time ne n_in_time then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'The length of the time vector for input variable ' $
          + in_var_label[i_var] + ' (' + str( n_temp_in_time ) $
          + ' ) differs from the length for other variables (' $
          + str( n_in_time ) + ').'
    endif
    if max( abs( temp_in_time_data - in_time_data ) ) gt 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Values in the time vector for input variable ' $
          + in_var_label[i_var] + ' differ from values for other variables.'
    endif
  endelse
  ; Determine if there are values masked to NaN
  id_nan = where( finite( temp_in_var_data ) eq 0, n_id_nan )
  ; Generate the ocean mask if we do not have one and this looks useful
  if not( keyword_set( ocean_mask_data ) ) then begin
    temp = float( n_id_nan ) / n_elements( temp_in_var_data )
    if ( temp gt 0.25 ) and ( temp lt 0.35 ) then begin
      ocean_mask_data = float( finite( temp_in_var_data[*,*,0] ) )
    endif
  endif
  ; Fill in values masked to NaN
  if n_id_nan gt 0 then c20c_dtos_v2_unnan, temp_in_var_data, wrap_x=1
  ; Copy this variable's data into the data array
  in_var_data[*,*,*,i_var] = temporary( temp_in_var_data )
  ; Copy global attribute
  if n_id_global gt 0 then begin
    out_global_attribute_value[id_global] = temp_global_value
    out_global_attribute_type[id_global] = var_type( temp_global_value )
  endif
endfor

; Copy output dimensions
n_out_lon = n_in_lon
out_lon_data = in_lon_data
n_out_lat = n_in_lat
out_lat_data = in_lat_data
n_out_time = n_in_time
out_time_data = double( in_time_data )
temp_label = [ 'units', 'calendar' ]
temp_value = [ in_time_units, in_time_calendar ]
if n_out_time_attribute eq 0 then begin
  out_time_attribute_label = temp_label
  out_time_attribute_value = temp_value
endif else begin
  id_time = where( out_time_label eq 'time', n_id_time )
  if n_id_time gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple instances of time variable in output time variable list.'
  endif
  if n_id_time eq 1 then begin
    for i_attribute = 0, n_elements( temp_label ) - 1 do begin
      id = where( out_time_attribute_label[*,id_time[0]] $
          eq temp_label[i_attribute], n_id )
      if n_id gt 1 then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'Multiple instances of ' + temp_label[i_attribute] $
            + ' attribute for time variable ' + out_time_label[id_time[0]] + '.'
      endif
      if n_id eq 1 then begin
        if out_time_attribute_value[id[0],id_time[0]] eq '' then begin
          out_time_attribute_value[id[0],id_time[0]] = temp_value[i_attribute]
        endif
      endif else begin
        temp = strarr( 1, n_out_time_var )
        temp[0,id_time[0]] = temp_label[i_attribute]
        out_time_attribute_label = [ out_time_attribute_label, temp ]
        temp[0,id_time[0]] = temp_value[i_attribute]
        out_time_attribute_value = [ out_time_attribute_value, temp ]
      endelse
    endfor
  endif
endelse
; Copy input data to active (potential) output data array
out_var_data = nan * fltarr( n_out_lon, n_out_lat, n_out_time, n_out_var )
for i_var = 0, n_act_var - 1 do begin
  id = where( in_var_label eq act_var_label[i_var], n_id )
  if n_id gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple instances of ' + act_var_label[i_var] $
        + ' variable in input variable list.'
  endif
  if n_id eq 1 then begin
    out_var_data[*,*,*,i_var] = in_var_data[*,*,*,id[0]]
    id_out = where( out_var_label eq act_var_label[i_var], n_id_out )
    if n_id_out gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Multiple instances of ' + act_var_label[i_var] $
          + ' variable in output variable list.'
    endif
    if n_id_out eq 1 then begin
      if out_var_type[id_out[0]] eq -1 then begin
        out_var_type[id_out[0]] = var_type( in_var_data )
      endif
    endif
  endif
endfor

; Fill in any missing output variable attribute values with those from the 
; input file.
if n_out_var_attribute gt 0 then begin
  ; Iterate through data variables
  for i_var = 0, n_out_var - 1 do begin
    ; Determine if we have any missing attribute values, omitting some if we 
    ; will be taking the CF standards
    if cf_standard_opt eq 1 then begin
      id_attribute = where( ( out_var_attribute_label[*,i_var] ne '' ) $
          and ( out_var_attribute_value[*,i_var] eq '' ) $
          and ( out_var_attribute_label[*,i_var] ne 'standard_name' ) $
          and ( out_var_attribute_label[*,i_var] ne 'long_name' ), $
          n_id_attribute )
    endif else begin
      id_attribute = where( ( out_var_attribute_label[*,i_var] ne '' ) $
          and ( out_var_attribute_value[*,i_var] eq '' ), n_id_attribute )
    endelse
    if n_id_attribute gt 0 then begin
      ; If the output variable is present in the input variable list then adopt 
      ; it
      if max( in_var_label eq out_var_label[i_var] ) eq 1 then begin
        temp_var_label = out_var_label[i_var]
      ; If a non-bcs version of the output variable is present in the input 
      ; list then adopt the non-bcs variable but without a few attributes
      endif else if max( in_var_label + 'bcs' eq out_var_label[i_var] ) eq 1 $
          then begin
        temp_var_label = strmid( out_var_label[i_var], 0, $
            strlen( out_var_label[i_var] ) - 3 )
        temp = [ 'comment', 'long_name', 'orig_name', 'standard_name' ]
        id = where( isin( temp, out_var_attribute_label[id_attribute,i_var] ) $
            eq 1, n_id_attribute )
        if n_id_attribute gt 0 then id_attribute = id_attribute[id]
      ; Otherwise we have no variable to work with
      endif else begin
        temp_var_label = ''
      endelse
      ; Proceed if we have an appropriate input variable and attribute list
      if ( temp_var_label ne '' ) and ( n_id_attribute gt 0 ) then begin
        ; Get the first file in the input file list (we will only check this 
        ; file)
        id_var = where( in_var_label eq out_var_label[i_var], n_id_var )
        if n_id_var gt 1 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
              + 'Multiple instances of ' + out_var_label[i_var] $
              + ' in input variable list.'
        endif
        if n_id_var eq 1 then begin
          temp_file = strsplit( in_file[0,id_var[0]], ',', extract=1 )
          temp_file = file_search( temp_file[0], count=n_temp_file )
          if n_temp_file eq 0 then begin
            stop, 'Error c20c_dtos_v2_make_tossic.pro:  Cannot find file (' $
                + temp_file + ') for use in copying ' + temp_var_label $
                + ' variable attributes.'
          endif
          temp_file = temp_file[0]
          ; Read the requested attributes from the input file
          temp_attribute_value = ''
          temp_attribute_type = -1
          temp = netcdf_read( temp_file[0], temp_var_label, no_data=1, $
              attribute_label=out_var_attribute_label[id_attribute,i_var], $
              attribute_value=temp_attribute_value, $
              attribute_type=temp_attribute_type )
          out_var_attribute_value[id_attribute,i_var] = temp_attribute_value
          out_var_attribute_type[id_attribute,i_var] = temp_attribute_type
        endif
      endif
    endif
  endfor
endif
; Fill in any missing longitude variable attribute values with those from the 
; input file.
if n_out_lon_attribute gt 0 then begin
  ; Determine if we have any missing attribute values, omitting some if we 
  ; will be taking the CF standards
  if cf_standard_opt eq 1 then begin
    id_attribute = where( ( out_lon_attribute_label ne '' ) $
        and ( out_lon_attribute_value eq '' ) $
          and ( out_lon_attribute_label ne 'standard_name' ) $
          and ( out_lon_attribute_label ne 'long_name' ), n_id_attribute )
  endif else begin
    id_attribute = where( ( out_lon_attribute_label ne '' ) $
        and ( out_lon_attribute_value eq '' ), n_id_attribute )
  endelse
  if n_id_attribute gt 0 then begin
    ; Get the first file in the input file list (we will only check this file)
    temp_file = strsplit( in_file[0,0], ',', extract=1 )
    temp_file = file_search( temp_file, count=n_temp_file )
    if n_temp_file eq 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  Cannot find file (' $
          + temp_file + ') for use in copying longitude variable attributes.'
    endif
    temp_file = temp_file[0]
    ; Read the requested attributes from the input file
    temp_attribute_value = ''
    temp = netcdf_read( temp_file[0], 'lon', no_data=1, $
        attribute_label=out_lon_attribute_label[id_attribute], $
        attribute_value=temp_attribute_value )
    out_lon_attribute_value[id_attribute] = temp_attribute_value
  endif
endif
; Fill in any missing latitude variable attribute values with those from the 
; input file.
if n_out_lat_attribute gt 0 then begin
  ; Determine if we have any missing attribute values, omitting some if we 
  ; will be taking the CF standards
  if cf_standard_opt eq 1 then begin
    id_attribute = where( ( out_lat_attribute_label ne '' ) $
        and ( out_lat_attribute_value eq '' ) $
          and ( out_lat_attribute_label ne 'standard_name' ) $
          and ( out_lat_attribute_label ne 'long_name' ), n_id_attribute )
  endif else begin
    id_attribute = where( ( out_lat_attribute_label ne '' ) $
        and ( out_lat_attribute_value eq '' ), n_id_attribute )
  endelse
  if n_id_attribute gt 0 then begin
    ; Get the first file in the input file list (we will only check this file)
    temp_file = strsplit( in_file[0,0], ',', extract=1 )
    temp_file = file_search( temp_file, count=n_temp_file )
    if n_temp_file eq 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  Cannot find file (' $
          + temp_file + ') for use in copying latitude variable attributes.'
    endif
    temp_file = temp_file[0]
    ; Read the requested attributes from the input file
    temp_attribute_value = ''
    temp = netcdf_read( temp_file[0], 'lat', no_data=1, $
        attribute_label=out_lat_attribute_label[id_attribute], $
        attribute_value=temp_attribute_value )
    out_lat_attribute_value[id_attribute] = temp_attribute_value
  endif
endif
; Fill in any missing time variable attribute values with those from the input 
; file.
if n_out_time_attribute gt 0 then begin
  ; Iterate through time variables
  for i_var = 0, n_out_time_var - 1 do begin
    ; Determine if we have any missing attribute values, omitting some if we 
    ; will be taking the CF standards
    if cf_standard_opt eq 1 then begin
      id_attribute = where( ( out_time_attribute_label[*,i_var] ne '' ) $
          and ( out_time_attribute_value[*,i_var] eq '' ) $
          and ( out_time_attribute_label[*,i_var] ne 'standard_name' ) $
          and ( out_time_attribute_label[*,i_var] ne 'long_name' ), $
          n_id_attribute )
    endif else begin
      id_attribute = where( ( out_time_attribute_label[*,i_var] ne '' ) $
          and ( out_time_attribute_value[*,i_var] eq '' ), n_id_attribute )
    endelse
    if n_id_attribute gt 0 then begin
      ; If the output variable is present in the input variable list then 
      ; proceed
      if max( in_time_label eq out_time_label[i_var] ) eq 1 then begin
        ; Get the first file in the input file list (we will only check this 
        ; file)
        temp_file = strsplit( in_file[0,0], ',', extract=1 )
        temp_file = file_search( temp_file, count=n_temp_file )
        if n_temp_file eq 0 then begin
          stop, 'Error c20c_dtos_v2_make_tossic.pro:  Cannot find file (' $
              + temp_file + ') for use in copying time variable attributes.'
        endif
        temp_file = temp_file[0]
        ; Read the requested attributes from the input file
        temp_attribute_value = ''
        temp = netcdf_read( temp_file[0], out_time_label[i_var], no_data=1, $
            attribute_label=out_time_attribute_label[id_attribute,i_var], $
            attribute_value=temp_attribute_value )
        out_time_attribute_value[id_attribute,i_var] = temp_attribute_value
      endif
    endif
  endfor
endif

;***********************************************************************
; Adjust according to attributable difference field

; If subtraction of a difference field has been requested
if n_delta gt 0 then begin
  ; Iterate through delta patterns
  for i_delta = 0, n_delta - 1 do begin
    ; Load difference data, including interpolation to target grid
    temp_period = 0
    if keyword_set( delta_period ) then temp_period = delta_period[*,i_delta]
    temp_file = strsplit( delta_file[*,i_delta], ',', extract=1 )
    if keyword_set( delta_var_orig_label ) then begin
      temp_delta_var_orig_label = delta_var_orig_label[i_delta]
    endif
    temp_delta_time_units = 0
    temp_delta_time_calendar = 0
    temp_delta_var_data = netcdf_read_geo_multitime( temp_file, $
        delta_var_label[i_delta], period_time=temp_period, $
        lon=temp_delta_lon_data, lat=temp_delta_lat_data, $
        time=temp_delta_time_data, units_var=temp_delta_var_units, quiet=1, $
        fix_time=1, label_in_var=temp_delta_var_orig_label, $
        no_units_conversion=1, units_time=temp_delta_time_units, $
        calendar=temp_delta_time_calendar )
    ;if temp_delta_var_units ne temp_in_var_units then begin
    ;  check = 1
    ;  temp = [ 'celsius', 'degrees c', 'deg_c', 'degc' ]
    ;  if ( max( strlowcase( temp_delta_var_units ) eq temp ) eq 1 ) $
    ;      and ( temp_in_var_units[0] eq 'K' ) then begin
    ;    check = 0
    ;  endif else if ( max( strlowcase( temp_in_var_units ) eq temp ) eq 1 ) $
    ;      and ( temp_delta_var_units[0] eq 'K' ) then begin
    ;    check = 0
    ;  endif
    ;  if check eq 1 then stop
    ;endif
    n_temp_delta_lon = n_elements( temp_delta_lon_data )
    n_temp_delta_lat = n_elements( temp_delta_lat_data )
    n_temp_delta_time = n_elements( temp_delta_time_data )
    temp_delta_var_data = reform( temp_delta_var_data, n_temp_delta_lon, $
        n_temp_delta_lat, n_temp_delta_time )
    ; Confirm consistencey with previously loaded deltas
    if i_delta gt 0 then begin
      if n_temp_delta_time ne n_delta_time then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'New delta data (' + str( i_delta ) $
            + ') has a different time length (' + str( n_temp_delta_time ) $
            + ') than the previous delta data (' + n_delta_time + ').'
      endif
      if temp_delta_var_units ne delta_var_units then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'New delta data units (' + temp_delta_var_units $
            + ') differ from previous delta data units (' + delta_var_units $
            + ').'
      endif
    endif
    ; Interpolate to output grid
    check_regrid = 0
    if n_temp_delta_lon ne n_in_lon then begin
      check_regrid = 1
    endif else if n_temp_delta_lat ne n_in_lat then begin
      check_regrid = 1
    endif else if max( abs( temp_delta_lon_data - in_lon_data ) ) gt 0 $
        then begin
      check_regrid = 1
    endif else if max( abs( temp_delta_lat_data - in_lat_data ) ) gt 0 $
        then begin
      check_regrid = 1
    endif
    if check_regrid eq 1 then begin
      mask_lonlattime, temp_delta_var_data, lon=temp_delta_lon_data, $
          lat=temp_delta_lat_data, mask_lon=in_lon_data, mask_lat=in_lat_data
    ; Otherwise confirm the grid is consistent with existing data
    endif
    ; Fill in values masked to NaN
    id_nan = where( finite( temp_delta_var_data ) eq 0, n_id_nan )
    if n_id_nan gt 0 then c20c_dtos_v2_unnan, temp_delta_var_data, wrap_x=1
    ; Apply scaling factor to delta pattern
    temp_delta_var_data = temp_delta_var_data * delta_factor[i_delta]
    ; Extend annual cycle if requested
    if keyword_set( delta_extend_cycle_opt ) then begin
      ; This currently only works with monthly data over one year
      if n_temp_delta_time ne 12 then stop
      ; Convert time data to yyyymmdd format
      temp_delta_time_data_ymd = convert_time_format( temp_delta_time_data, $
          temp_delta_time_units, 'yyyymmdd', calendar=temp_delta_time_calendar )
      out_time_data_ymd = convert_time_format( out_time_data, in_time_units, $
          'yyyymmdd', calendar=in_time_calendar )
      ; Initialise new delta data array
      temp_data = !values.f_nan * fltarr( n_out_lon, n_out_lat, n_out_time )
      ; Iterate through months in the new delta data array
      for i_time = 0, n_out_time - 1 do begin
        ; Identify month in new delta data array
        temp_month = strmid( out_time_data_ymd[i_time], 4, 2 )
        ; Identify month in annual delta data
        id_delta = where( $
            strmid( temp_delta_time_data_ymd, 4, 2 ) eq temp_month, n_id_delta )
        if n_id_delta ne 1 then stop
        ; Copy data
        temp_data[*,*,i_time] = temp_delta_var_data[*,*,id_delta[0]]
      endfor
      ; Adopt new delta data
      temp_delta_time_data = out_time_data
      temp_delta_var_data = temporary( temp_data )
    endif
    ; Record delta
    if i_delta eq 0 then begin
      delta_var_data = temporary( temp_delta_var_data )
      n_delta_time = n_temp_delta_time
      delta_var_units = temporary( temp_delta_var_units )
    endif else begin
      delta_var_data = delta_var_data + temp_delta_var_data
    endelse
  endfor
  ; Add the delta from the specified input variable
  id_var = where( isin( delta_apply_var_label, act_var_label ) eq 1, n_id_var )
  for i_var = 0, n_id_var - 1 do begin
    out_var_data[*,*,*,id_var[i_var]] = out_var_data[*,*,*,id_var[i_var]] $
        + delta_var_data
  endfor
  ; If we need to alter sea ice concentration data for consistency with 
  ; altered sea surface temperature data
  if ( max( delta_apply_var_label eq 'tos' ) eq 1 ) $
      and ( max( delta_apply_var_label eq 'sic' ) eq 0 ) $
      and ( max( out_var_label eq 'sic' ) eq 1 ) $
      and keyword_set( tossic_fit_method ) then begin
    ; Identify the sea ice concentration variable
    id_sic = where( in_var_label eq 'sic', n_id_sic )
    if n_id_sic ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot input sic data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    ; Determine sea ice units
    id_sic_units = where( out_var_attribute_label[*,id_sic[0]] eq 'units', $
        n_id )
    if n_id ne 1 then stop
    out_sic_units = out_var_attribute_value[id_sic_units[0],id_sic[0]]
    ; Identify the sea surface temperature variable
    id_tos = where( in_var_label eq 'tos', n_id_tos )
    if n_id_tos ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot input tos data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    ; Calculate adjusted SIC
    in_time_data_yyyymmdd = convert_time_format( in_time_data, in_time_units, $
        'yyyymmdd', calendar=in_time_calendar )
    c20c_dtos_v2_adjust_sic, fit_sic_file=tossic_fit_sic_file, $
        fit_tos_file=tossic_fit_tos_file, fit_method=tossic_fit_method, $
        fit_period=tossic_fit_period, $
        in_tos_data=in_var_data[*,*,*,id_tos[0]], $
        in_sic_data=in_var_data[*,*,*,id_sic[0]], in_lon=in_lon_data, $
        in_lat=in_lat_data, in_time=in_time_data_yyyymmdd, $
        delta_tos_data=delta_var_data, ocean_mask_data=ocean_mask_data, $
        out_sic_data=out_sic_data, v1=v1_opt
    ; Convert back to % (output above is in fraction)
    if out_sic_units ne 'fraction' then begin
      if out_sic_units eq '%' then begin
        temp = max( out_sic_data, nan=1 )
        if ( temp gt 0.5 ) and ( temp lt 2. ) then begin
          out_sic_data = out_sic_data * 100.
        endif
      endif else begin
        stop
      endelse
    endif
    ; Fill in values masked to NaN
    id_nan = where( finite( out_sic_data ) eq 0, n_id_nan )
    if n_id_nan gt 0 then c20c_dtos_v2_unnan, out_sic_data, wrap_x=1
    ; Copy data to active output array
    id_sic = where( act_var_label eq 'sic', n_id_sic )
    if n_id_sic ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot record sic data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    out_var_data[*,*,*,id_sic[0]] = out_sic_data
  endif
  if ( max( delta_apply_var_label eq 'tosbcs' ) eq 1 ) $
      and ( max( delta_apply_var_label eq 'sicbcs' ) eq 0 ) $
      and ( max( out_var_label eq 'sicbcs' ) eq 1 ) $
      and keyword_set( tossic_fit_method ) then begin
    ; Identify the sea ice concentration variable
    id_sic = where( in_var_label eq 'sicbcs', n_id_sic )
    if n_id_sic ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot input sicbcs data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    ; Determine sea ice units
    id_sic_units = where( out_var_attribute_label[*,id_sic[0]] eq 'units', $
        n_id )
    if n_id ne 1 then stop
    out_sicbcs_units = out_var_attribute_value[id_sic_units[0],id_sic[0]]
    ; Identify the sea surface temperature variable
    id_tos = where( in_var_label eq 'tosbcs', n_id_tos )
    if n_id_tos ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot input tosbcs data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    ; Calculate adjusted SIC
    in_time_data_yyyymmdd = convert_time_format( in_time_data, in_time_units, $
        'yyyymmdd', calendar=in_time_calendar ) 
    c20c_dtos_v2_adjust_sic, fit_sic_file=tossic_fit_sic_file, $
        fit_tos_file=tossic_fit_tos_file, fit_method=tossic_fit_method, $
        fit_period=tossic_fit_period, $
        in_tos_data=in_var_data[*,*,*,id_tos[0]], $
        in_sic_data=in_var_data[*,*,*,id_sic[0]], in_lon=in_lon_data, $
        in_lat=in_lat_data, in_time=in_time_data_yyyymmdd, $
        delta_tos_data=delta_var_data, ocean_mask_data=ocean_mask_data, $
        out_sic_data=out_sic_data, v1=v1_opt
    ;; Convert back to % (output above is in fraction)
    ;if out_sicbcs_units ne 'fraction' then begin
    ;  if out_sicbcs_units eq '%' then begin
    ;    temp = max( out_sic_data, nan=1 )
    ;    if ( temp gt 0.5 ) and ( temp lt 2. ) then begin
    ;      out_sic_data = out_sic_data * 100.
    ;    endif
    ;  endif else begin
    ;    stop
    ;  endelse
    ;endif
    ; Fill in values masked to NaN
    id_nan = where( finite( out_sic_data ) eq 0, n_id_nan )
    if n_id_nan gt 0 then c20c_dtos_v2_unnan, out_sic_data, wrap_x=1
    ; Copy data to active output array
    id_sic = where( act_var_label eq 'sicbcs', n_id_sic )
    if n_id_sic ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot record sicbcs data because no such variable is defined ' $
          + 'or multiple variables are defined with that label.'
    endif
    out_var_data[*,*,*,id_sic[0]] = out_sic_data
  endif
endif

;***********************************************************************
; Perform mean-preserving interpolation

; Iterate through output variables
for i_var = 0, n_out_var - 1 do begin
  ; If this is a variable requiring the Sheng-Zwiers adjustment
  if ( strpos( out_var_label[i_var], 'bcs', reverse_search=1 ) $
      eq strlen( out_var_label[i_var] ) - 3 ) $
      and ( max( in_var_label eq out_var_label[i_var] ) eq 0 ) then begin
    ; Find the source variable
    id_source = where( act_var_label + 'bcs' eq out_var_label[i_var], n_id )
    if n_id ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'OUT_VAR_LABEL[' + str( i_var ) $
          + '] defines a "bcs" variable to be calculated, but no non-bcs ' $
          + 'version exists in the current list of variables.'
    endif
    ; Perform adjustment
    temp_data = reform( out_var_data[*,*,*,id_source[0]], $
        n_out_lon * n_out_lat, n_out_time )
    temp_out_data = shengzwiers( temp_data, out_time_data, $
        extend_cycle=shengzwiers_extend_cycle, double=double_opt )
    temp_data = 0
    out_var_data[*,*,*,i_var] = reform( temporary( temp_out_data ), n_out_lon, $
        n_out_lat, n_out_time )
    ; Note variable type if necessary
    if out_var_type[i_var] eq -1 then begin
      out_var_type[i_var] = var_type( out_var_data )
    endif
    ;; Add comment variable attribute
    ;temp_comment = 'Adjusted to return original values upon averaging of ' $
    ;    + 'linearly interpolated values, following Sheng and Zwiers (1998, ' $
    ;    + 'Climate Dynamics, 14, 609-613).'
    ;id = where( out_var_attribute_label[*,i_var] eq 'comment', n_id )
    ;if n_id gt 1 then stop
    ;if n_id eq 1 then begin
    ;  if out_var_attribute_value[id[0],i_var] eq '' then begin
    ;    out_var_attribute_value[id[0],i_var] = temp_comment
    ;  endif else begin
    ;    out_var_attribute_value[id[0],i_var] $
    ;        = out_var_attribute_value[id[0],i_var] + hard_return + ' ' $
    ;        + temp_comment
    ;  endelse
    ;endif else begin
    ;  temp = strarr( 1, n_out_var )
    ;  temp[0,i_var] = 'comment'
    ;  out_var_attribute_label = [ out_var_attribute_label, temp ]
    ;  temp = strarr( 1, n_out_var )
    ;  temp[0,i_var] = temp_comment
    ;  out_var_attribute_value = [ out_var_attribute_value, temp ]
    ;endelse
  endif
endfor

; Confirm adjusted freezing temperatures are numerically at the defined 
; freezing temperature
id_tos = where( isin( [ 'tos', 'tosbcs' ], out_var_label ) eq 1, n_id_tos )
if n_id_tos ge 1 then begin
  out_var_data = reform( out_var_data, n_out_lon * n_out_lat * n_out_time, $
      n_out_var )
  for i_id_tos = 0, n_id_tos - 1 do begin
    id_units = where( out_var_attribute_label[*,id_tos[i_id_tos]] eq 'units', $
        n_id )
    if n_id ne 1 then stop
    temp_value = out_var_attribute_value[id_units[0],id_tos[i_id_tos]]
    if temp_value eq 'K' then begin
      temp_tos_freeze = tos_freeze
    endif else if temp_value eq 'deg_C' then begin
      temp_tos_freeze = tos_freeze - 273.15
    endif else begin
      stop
    endelse
    id = where( out_var_data[*,id_tos[i_id_tos]] lt temp_tos_freeze + 0.0001, $
        n_id )
    if n_id gt 0 then out_var_data[id,id_tos[i_id_tos]] = temp_tos_freeze
  endfor
  out_var_data = reform( out_var_data, n_out_lon, n_out_lat, n_out_time, $
      n_out_var )
endif

;***********************************************************************
; Prepare output variables and metadata

; Restrict output data array to requested output variables
id_out = intarr( n_out_var )
for i_var = 0, n_out_var - 1 do begin
  id_act = where( act_var_label eq out_var_label[i_var], n_id_act )
  if n_id_act ne 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Variable ' + out_var_label[i_var] + ' does not exist in list of ' $
        + 'possible output variables.'
  endif
  id_out[i_var] = id_act[0]
endfor
out_var_data = reform( out_var_data[*,*,*,id_out] )

; Ensure output data variable data is in desired format
if keyword_set( out_var_attribute_label ) then begin
  for i_var = 0, n_out_var - 1 do begin
    ; Copy the data variable units attribute if provided
    id = where( out_var_attribute_label[*,i_var] eq 'units', n_id )
    if n_id gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Multiple units attribute found for output variable ' $
          + out_var_label[i_var] + '.'
    endif
    if n_id eq 1 then begin
      temp_out_var_units = out_var_attribute_value[id[0],i_var]
      ; Kelvin to other
      if max( out_var_label[i_var] eq [ 'tos', 'tosbcs' ] ) eq 1 then begin
        if temp_out_var_units ne 'K' then begin
          if temp_out_var_units eq 'deg_C' then begin
            out_var_data[*,*,*,i_var] = out_var_data[*,*,*,i_var] - 273.15
          endif else begin
            stop
          endelse
        endif
      ; fraction to other
      endif else if max( out_var_label[i_var] eq [ 'sic', 'sicbcs' ] ) eq 1 $
          then begin
        if temp_out_var_units ne 'fraction' then begin
          if temp_out_var_units eq '%' then begin
            out_var_data[*,*,*,i_var] = out_var_data[*,*,*,i_var] / 100.
          endif else begin
            stop
          endelse
        endif
      ; Not yet implemented
      endif else begin
        stop
      endelse
    endif
  endfor
endif

; Restrict output data to requested period
if n_elements( out_period ) eq 2 then begin
  check = 0
  if n_elements( in_period ) eq 0 then begin
    check = 1
  endif else begin
    if ( out_period[0] ne in_period[0,0] ) $
        and ( out_period[1] ne in_period[1,0] ) then begin
      check = 1
    endif
  endelse
  if check eq 1 then begin
    id_time = where( out_time_label eq 'time', n_id_time )
    if n_id_time ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot restrict to OUT_PERIOD because no time variable defined.'
    endif
    id_units = where( out_time_attribute_label[*,id_time[0]] eq 'units', $
        n_id_units )
    if n_id_units ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot restrict to OUT_PERIOD because no time units attribute ' $
          + 'defined.'
    endif
    id_calendar = where( out_time_attribute_label[*,id_time[0]] eq 'calendar', $
        n_id_calendar )
    if n_id_calendar ne 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Cannot restrict to OUT_PERIOD because no time calendar ' $
          + 'attribute defined.'
    endif
    out_period_dayssince = convert_time_format( out_period, 'yyyymmdd', $
        out_time_attribute_value[id_units[0],id_time[0]], $
        calendar=out_time_attribute_value[id_calendar[0],id_time[0]] )
    id_time = where( ( out_time_data ge out_period_dayssince[0] ) $
        and ( out_time_data le out_period_dayssince[1] ), n_out_time )
    if n_out_time eq 0 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Output time values are all outside range of OUT_PERIOD.'
    endif
    out_var_data = out_var_data[*,*,id_time,*]
    out_time_data = out_time_data[id_time]
    id_time = 0
  endif
endif

; If there are more time-related variables to output
if n_out_time_var gt 1 then begin
  ; Modify output time vector into array for all time-related variables
  out_time_data = add_dim( out_time_data, 1, n_out_time_var )
  id = where( out_time_label ne 'time', n_id )
  if n_id eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple output time variables defined.'
  endif
  if n_id eq n_out_time_var then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No output time variable specified.'
  endif
  out_time_data[*,id] = nan
  id_time = where( out_time_label eq 'time', n_id_time )
  if n_id_time eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No output time variable specified.'
  endif else if n_id_time gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple output time variables defined.'
  endif
  id_units = where( out_time_attribute_label[*,id_time[0]] eq 'units', $
      n_id_units )
  if n_id_units eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No time units attribute defined or found.'
  endif else if n_id_units gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  Multiple time units attributes.'
  endif
  id_calendar = where( out_time_attribute_label[*,id_time[0]] eq 'calendar', $
      n_id_calendar )
  if n_id_calendar eq 0 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'No time calendar attribute defined or found.'
  endif else if n_id_calendar gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple time calendar attributes.'
  endif
  ; Create date and/or datesec output variables if requested.
  id_date = where ( out_time_label eq 'date', n_id_date )
  if n_id_date gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple date time variables defined in OUT_TIME_LABEL.'
  endif
  id_datesec = where ( out_time_label eq 'datesec', n_id_datesec )
  if n_id_datesec gt 1 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'Multiple datesec time variables defined in OUT_TIME_LABEL.'
  endif
  if ( n_id_date eq 1 ) or ( n_id_datesec eq 1 ) then begin
    ; Determine date and time
    date_time = convert_time_format( out_time_data[*,id_time[0]], $
        out_time_attribute_value[id_units[0],id_time[0]], 'yyyymmddhhmmss', $
        calendar=out_time_attribute_value[id_calendar[0],id_time[0]] )
    ; Extract date
    if n_id_date eq 1 then begin
      out_time_data[*,id_date[0]] = long( strmid( date_time, 0, 8 ) )
      ; Define variable type if necessary
      if out_time_type[id_date[0]] eq -1 then out_time_type[id_date[0]] = 3
    endif
    ; Extract datesec
    if n_id_datesec eq 1 then begin
      out_time_data[*,id_datesec[0]] $
          = long( strmid( date_time, 8, 2 ) ) * 60l * 60l $
          + long( strmid( date_time, 10, 2 ) ) * 60l $
          + long( strmid( date_time, 12, 2 ) )
      ; Define variable type if necessary
      if out_time_type[id_datesec[0]] eq -1 then begin
        out_time_type[id_datesec[0]] = 3
      endif
    endif
  endif
endif

; Ensure CF-standard metadata is included for the data variables, if requested
if cf_standard_opt eq 1 then begin
  ; Add any missing CF-standard attributes, assuming they will be filled in 
  ; with defaults by netcdf_write.pro
  cf_label = [ 'units', 'long_name', 'standard_name' ]
  for i_var = 0, n_out_var - 1 do begin
    temp = isin( out_var_attribute_label[*,i_var], cf_label )
    id = where( temp eq 0, n_id )
    if n_id gt 0 then begin
      temp = strarr( n_id, n_out_var )
      temp[*,i_var] = cv_label[id]
      out_var_attribute_label = [ out_var_attribute_label, temp ]
      out_var_attribute_value = [ out_var_attribute_value, $
          strarr( n_id, n_out_var ) ]
    endif
  endfor
endif
; Compress empty attributes in attribute arrays
if n_out_var_attribute gt 0 then begin
  for i_var = 0, n_out_var - 1 do begin
    id = where( out_var_attribute_label[*,i_var] ne '', n_id )
    if n_id gt 0 then begin
      out_var_attribute_label[0:n_id-1,i_var] $
          = out_var_attribute_label[id,i_var]
      out_var_attribute_value[0:n_id-1,i_var] $
          = out_var_attribute_value[id,i_var]
    endif
  endfor
  if n_out_var eq 1 then begin
    id = where( strlen( out_var_attribute_label ) ne 0, n_id )
  endif else begin
    id = where( total( strlen( out_var_attribute_label ), 2, integer=1 ) ne 0, $
        n_id )
  endelse
  if n_id gt 0 then begin
    out_var_attribute_label = out_var_attribute_label[id,*]
    out_var_attribute_value = out_var_attribute_value[id,*]
  endif
endif

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
  if n_out_file_parsed ne 8 then begin
    stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
        + 'OUT_FILE should have 8 components separated by "_" (not ' $
        + str( n_out_file_parsed ) + ').'
  endif
  if max( out_global_attribute_label eq 'realm' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'realm' ]
    if max( out_file_parsed[1] eq [ 'O3hr', 'Oday', 'Omon' ] ) eq 1 then begin
      out_global_attribute_value = [ out_global_attribute_value, 'ocean' ]
    endif else if max( out_file_parsed[1] eq [ 'OI3hr', 'OIday', 'OImon' ] ) $
        eq 1 then begin
      out_global_attribute_value = [ out_global_attribute_value, 'seaIce' ]
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Unrecognised realm defined through OUT_FILE (' $
          + out_file_parsed[1] + ').'
    endelse
    out_global_attribute_type = [ out_global_attribute_type, 7 ]
  endif
  if max( out_global_attribute_label eq 'frequency' ) eq 0 then begin
    out_global_attribute_label = [ out_global_attribute_label, 'frequency' ]
    if max( out_file_parsed[1] eq [ 'Oday', 'OIday' ] ) eq 1 then begin
      out_global_attribute_value = [ out_global_attribute_value, 'day' ]
    endif else if max( out_file_parsed[1] eq [ 'Omon', 'OImon' ] ) eq 1 $
        then begin
      out_global_attribute_value = [ out_global_attribute_value, 'mon' ]
    endif else begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Unrecognised frequency defined through OUT_FILE (' $
          + out_file_parsed[1] + ').'
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
    print, 'WARNING (c20c_dtos_v2_make_tossic.pro):  The following global ' $
        + 'attributes are missing but should be present:  ' $
        + string_from_vector( temp_global[id], addand=1 ) + '.'
  endif
  ; Add comment on sea ice concentration estimation
  if keyword_set( tossic_fit_method ) then begin
    temp_comment = ''
    if max( tossic_fit_method eq 'pall2007' ) eq 1 then begin
      temp_comment = 'Attributable sea ice concentration change estimated ' $
         + 'according to: ' + hard_return $
         + ' Pall, P.  2007.  Constraints on, and attribution of, changes in ' $
         + 'extreme precipitation under climate change.  Ph.D. thesis, St. ' $
         + 'Cross College, University of Oxford, Oxford, U.K., 187pp. ' $
         + hard_return $
         + ' Pall, P., T. Aina, D. A. Stone, P. A. Stott, T. Nozawa, ' $
         + 'A. G. J. Hilberts, D. Lohmann, and M. R. Allen.  2011.  ' $
         + 'Anthropogenic greenhouse gas contribution to flood risk in ' $
         + 'England and Wales in Autumn 2000.  Nature, 470, 382-385.'
    endif else if max( tossic_fit_method eq 'stonepall2018' ) eq 1 then begin
      temp_comment = 'Attributable sea ice concentration change estimated ' $
         + 'according to: ' + hard_return $
         + 'Stone, D. A., and P. Pall.  2018.  A benchmark estimate of the ' $
         + 'effect of anthropogenic emissions on the ocean surface.  ' $
         + 'Submitted to Geoscientific Model Development.'
    endif
    if temp_comment ne '' then begin
      id = where( out_global_attribute_label eq 'comment', n_id )
      if n_id gt 1 then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'Multiple global comment attributes defined.'
      endif
      if n_id eq 1 then begin
        out_global_attribute_value[id[0]] = out_global_attribute_value[id[0]] $
            + ' ' + hard_return + ' ' + temp_comment
        out_global_attribute_type[id[0]] = 7
      endif else begin
        out_global_attribute_label = [ out_global_attribute_label, 'comment' ]
        out_global_attribute_value = [ out_global_attribute_value, $
            temp_comment ]
        out_global_attribute_type = [ out_global_attribute_type, 7 ]
      endelse
    endif
  endif
endif

; Add creation date global attribute
creation_date = systime( utc=1 )
creation_date = strsplit( creation_date, ' ', extract=1, count=n_creation_date )
if n_creation_date ne 5 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'Something is non-standard about the output from systime().'
endif
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
if n_id gt 1 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
      + 'Multiple global history attributes defined.'
endif
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
; Write output to file

; Check if output directory exists, if one is specified
out_dir = strsplit( out_file, '/', extract=1, count=n_temp )
if n_temp gt 1 then begin
  out_dir = string_from_vector( out_dir[0:temp-2], nospace=1, spacer='/' )
  if strmid( out_file, 0, 1 ) eq '/' then out_dir = '/' + out_dir
  temp = file_search( out_dir, count=n_temp )
  ; If the output directory does not exist then create it
  if n_temp eq 0 then spawn, 'mkdir -p ' + out_dir
endif

; Ensure we are writing all requested variables
temp = total( total( total( finite( out_var_data ), $
    1, integer=1 ), 1, integer=1 ), 1, integer=1 )
if min( temp ) eq 0 then begin
  stop, 'Error c20c_dtos_v2_make_tossic.pro:  Missing output data....'
endif

; The default output time variable
if n_elements( out_time_data ) gt 0 then begin
  if not( keyword_set( out_time_label ) ) then begin
    temp_out_time_label = 'time'
  endif else begin
    temp_out_time_label = out_time_label
  endelse
endif

; Modify output variable labels if requested
out_var_label_use = out_var_label
if n_out_var_attribute gt 0 then begin
  out_var_attribute_label_use = out_var_attribute_label
  out_var_attribute_value_use = out_var_attribute_value
  out_var_attribute_type_use = out_var_attribute_type
  for i_var = 0, n_out_var - 1 do begin
    id = where( out_var_attribute_label_use[*,i_var] eq 'label', n_id )
    if n_id gt 1 then begin
      stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
          + 'Multiple "label" attributes defined for output variable ' $
          + out_var_label[i_var] + '.'
    endif
    if n_id eq 1 then begin
      if out_var_attribute_value_use[id[0],i_var] eq '' then begin
        stop, 'Error c20c_dtos_v2_make_tossic.pro:  ' $
            + 'Empty "label" attribute for output variable ' $
            + out_var_label[i_var] + '.'
      endif
      temp = out_var_label_use[i_var]
      out_var_label_use[i_var] = out_var_attribute_value_use[id[0],i_var]
      out_var_attribute_label_use[id[0],i_var] = 'name'
      out_var_attribute_value_use[id[0],i_var] = temp
    endif
  endfor
endif

; Write data to NetCDF file
netcdf_write, out_file, data_array=out_var_data, data_label=out_var_label_use, $
    data_attribute_label=out_var_attribute_label_use, $
    data_attribute_value=out_var_attribute_value_use, $
    data_attribute_type=out_var_attribute_type_use, data_type=out_var_type, $
    dim1_vector=out_lon_data, dim1_label='lon', $
    dim1_attribute_label=out_lon_attribute_label, $
    dim1_attribute_value=out_lon_attribute_value, dim2_vector=out_lat_data, $
    dim2_label='lat', dim2_attribute_label=out_lat_attribute_label, $
    dim2_attribute_value=out_lat_attribute_value, dim3_vector=out_time_data, $
    dim3_label=temp_out_time_label, dim3_type=out_time_type, $
    dim3_attribute_label=out_time_attribute_label, $
    dim3_attribute_value=out_time_attribute_value, $
    global_attribute_label=out_global_attribute_label, $
    global_attribute_value=out_global_attribute_value, $
    global_attribute_type=out_global_attribute_type

; Print progress
print, 'Results of c20c_dtos_v2_make_tossic.pro written to ' + out_file + '.'

;***********************************************************************
; The end

;stop
return
END
