;+
; NAME:
;    netcdf_read.pro
;
; PURPOSE:
;    This function reads and returns the desired data variable from the desired 
;    NetCDF file.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    data = netcdf_read( file_name, var_label )
;
; INPUTS:
;    FILE_NAME:  A required scalar string containing the name of the NetCDF 
;        file.
;    VAR_LABEL:  A required string containing the label of the variable to read 
;        from the NetCDF file.  This may be a vector of strings, but then only 
;        the first entry that is available in the file will be returned.  This 
;        is useful for instance if there are a number of different possible 
;        labels for a variable, e.g. [ 'lon', 'longitude' ] for longitude, and 
;        only one of them will be available in the file anyway.  Also returns 
;        the label of the variable that was actually read.
;    ATTRIBUTE_LABEL, GLOBAL_LABEL, MISSING_LABEL, MISSING_REPLACE
;
; KEYWORD PARAMETERS:
;    ATTRIBUTE_LABEL:  An optional vector string containing the labels of the 
;        attributes for the requested variable which should have their values 
;        returned in ATTRIBUTE_VALUE.  If set to 1 then all attributes will be 
;        returned.
;    ATTRIBUTE_TYPE:  If ATTRIBUTE_LABEL is input, then this returns the 
;        variable types of the attribute values returned in ATTRIBUTE_VALUE.  
;        For instance if two attributes are requested and one is a string and 
;        the other a short integer, then this would return [7,2].  Of same size 
;        as ATTRIBUTE_LABEL.
;    ATTRIBUTE_VALUE:  If ATTRIBUTE_LABEL is input, then this returns the 
;        values of the attributes requested in ATTRIBUTE_LABEL of the requested 
;        variable.  Of same size as ATTRIBUTE_LABEL.
;    GLOBAL_LABEL:  An optional vector string containing the labels of the 
;        global attributes which should have their values returned in 
;        GLOBAL_VALUE.  If set to 1 then all global attributes will be returned.
;    GLOBAL_TYPE:  If GLOBAL_LABEL is input, then this returns the variable 
;        type of the values returned in GLOBAL_VALUE, for instance "7" for 
;        a string.  Of same size as GLOBAL_LABEL, and of type integer.
;    GLOBAL_VALUE:  If GLOBAL_LABEL is input, then this returns the values of 
;        the global attributes requested in GLOBAL_LABEL.  Of same size as 
;        GLOBAL_LABEL.
;    MISSING_LABEL:  An optional string vector containing possible labels for 
;        the variable-specific missing_value attribute.  If the attribute is 
;        found, then the missing value is replaced with the value specified in 
;        MISSING_REPLACE if provided or else a default specified in the code.
;    MISSING_REPLACE:  If MISSING_LABEL is input, then this is an optional 
;        scalar containing the value with which to replace flagged missing 
;        values in the data array.
;    NO_DATA:  If set, then the variable's data is not returned, with only the 
;        requested metadata returned in the relevant keyword parameter outputs. 
;        The default is to return the data if a variable specified in VAR_LABEL 
;        exists in the file.
;    UNITS:  Returns the value of the units attribute for the variable.
;
; OUTPUTS:
;    DATA:  Returns an array containing the data from the requested variable 
;        in the requested file.
;    VAR_LABEL:  As well as being an input, this returns a scalar string 
;        containing the label of the variable that was actually read.
;    ATTRIBUTE_TYPE, ATTRIBUTE_VALUE, GLOBAL_TYEP, GLOBAL_VALUE, UNITS
;
; USES:
;    dimension.pro
;    var_type.pro
;
; PROCEDURE:
;    This function scans through variables in the specified file and returns 
;    data and/or information as requested for the requested variables.
;
; EXAMPLE:
;    lon = netcdf_read( 'file.nc', 'lon' )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2016-04-25, as 
;        ncdf_read.pro.
;    Modified:  DAS, 2017-10-09 (Branched from ncdf_read.pro;  standardised 
;        documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2017-11-08 (Changed no-data output to NaN instead of zero)
;    Modified:  DAS, 2017-12-24 (Added GLOBAL_TYPE keyword output)
;    Modified:  DAS, 2018-09-18 (Added ATTRIBUTE_TYPE keyword output)
;    Modified:  DAS, 2018-11-23 (Added ATTRIBUTE_LABEL=1 and GLOBAL_LABEL=1 
;        options)
;    Modified:  DAS, 2018-11-28 (Set GLOBAL_TYPE to return 7 for byte type)
;-

;***********************************************************************

FUNCTION NETCDF_READ, $
    FILE_NAME, $
    VAR_LABEL, $
    MISSING_LABEL=missing_label, MISSING_REPLACE=missing_replace, $
    ATTRIBUTE_LABEL=attribute_label, ATTRIBUTE_TYPE=attribute_type, $
      ATTRIBUTE_VALUE=attribute_value, $
    UNITS=units, $
    GLOBAL_LABEL=global_label, GLOBAL_TYPE=global_type, $
      GLOBAL_VALUE=global_value, GLOBAL_ALL=global_all_opt, $
    NO_DATA=no_data_opt

;***********************************************************************
; Constants

; Option not to return data values
no_data_opt = keyword_set( no_data_opt )

; The list of default replacements for the missing value flag for each variable 
; type
if n_elements( missing_replace ) eq 0 then begin
  missing_replace_list = [ [ '4', 'NaN' ], [ '5', 'NaN' ], [ '6', 'NaN' ], $
      [ '7', '' ], [ '9', 'NaN' ] ]
  missing_replace_list = transpose( missing_replace_list )
endif

;***********************************************************************
; Open the File

; Open the file
id_file = ncdf_open( file_name )

; Find out what is in the file
info = ncdf_inquire( id_file )
; The number of input dimensions
n_dim = info.ndims
; The number of variables
n_var = info.nvars

;; Find the dimension labels
;dim_label = strarr( n_dim )
;dim_len = intarr( n_dim )
;for i_dim = 0, n_dim - 1 do begin
;  ncdf_diminq, id_file, i_dim, temp_1, temp_2
;  dim_label[i_dim] = temp_1
;  dim_len[i_dim] = reverse( temp_2 )
;endfor

; Find the variable labels
var_label_file = strarr( n_var )
for i_var = 0, n_var - 1 do begin
  var_info = ncdf_varinq( id_file, i_var )
  var_label_file[i_var] = var_info.name
endfor

; Find the requested global attributes
n_global = n_elements( global_label )
if n_global gt 0 then begin
  ; Determine what global attributes are in the file
  temp = ncdf_inquire( id_file )
  n_global_file = temp.ngatts
  global_label_file = strarr( n_global_file )
  for i_global = 0, n_global_file - 1 do begin
    global_label_file[i_global] = ncdf_attname( id_file, i_global, global=1 )
  endfor
  ; Copy all global attributes if requested
  if var_type( global_label ) eq 2 then begin
    if global_label[0] ne 1 then stop
    global_label = global_label_file
    n_global = n_global_file
  endif
  ; Initialise the vector of global attribute values
  global_value = strarr( n_global )
  global_type = intarr( n_global )
  ; Iterate through requested global attributes
  for i_global = 0, n_global - 1 do begin
    ; If this global attribute exists in the file
    if max( global_label_file eq global_label[i_global] ) eq 1 then begin
      ncdf_attget, id_file, global_label[i_global], temp, global=1
      global_value[i_global] = strtrim( string( temp ), 2 )
      global_type[i_global] = var_type( temp )
      if global_type[i_global] eq 1 then global_type[i_global] = 7
    endif
  endfor
endif

;***********************************************************************
; Load Data from the File

; Iterate through variables in the file until we have found a variable
for i_var = 0, n_var - 1 do begin
  if not( keyword_set( data ) ) then begin
    ; Explore this variable only if requested
    id = where( var_label eq var_label_file[i_var], n_id )
    if n_id ne 0 then begin
      var_label = var_label[id[0]]
      ; Load data only if requested
      if no_data_opt eq 0 then begin
        ncdf_varget, id_file, i_var, data
        if dimension( data ) gt 0 then data = reform( data )
        ;; Convert integers to floating point
        ;if max( var_type( data ) eq [ 1, 7 ] ) eq 0 then begin
        ;  if max( var_type( data ) eq [ 4, 5 ] ) eq 0 then data = float( data )
        ;endif
        ;; Convert bytes to strings
        ;if var_type( data ) eq 1 then begin
        ;  ; This assumes that we have a vector of bytes
        ;  if dimension( data ) ne 2 then begin
        ;    if dimension( data ) eq 1 then begin
        ;      data = reform( data, n_elements( data ), 1 )
        ;    endif else begin
        ;      stop
        ;    endelse
        ;  endif
        ;  n_temp = n_elements( data[0,*] )
        ;  temp_data = strarr( n_temp )
        ;  for i_data = 0, n_temp - 1 do begin
        ;    temp_data[i_data] = string( data[*,i_data] )
        ;  endfor
        ;  data = temporary( temp_data )
        ;endif
      endif
      ; Count attribute requests
      n_attribute = n_elements( attribute_label )
      ; Get variable attribute information
      n_var_att = ( ncdf_varinq( id_file, var_label_file[i_var] ) ).natts
      ; Proceed only if there is attribute information
      if n_var_att gt 0 then begin
        ; Read attribute information
        var_att = strarr( n_var_att )
        for i_att = 0, n_var_att - 1 do begin
          var_att[i_att] = ncdf_attname( id_file, var_label_file[i_var], i_att )
        endfor
        ; Copy all variable attributes if requested
        if var_type( attribute_label ) eq 2 then begin
          if attribute_label[0] ne 1 then stop
          attribute_label = var_att
          n_attribute = n_var_att
        endif
        ; Initialise list of variable attribute values
        if n_attribute ne 0 then begin
          attribute_value = strarr( n_attribute )
          attribute_type = strarr( n_attribute )
        endif
        ; Check for missing value flag
        if keyword_set( missing_label ) and keyword_set( data ) then begin
          ; Proceed only if we have a replacement value
          check_missing = 0
          if n_elements( missing_replace ) ne 0 then begin
            check_missing = 1
          endif else begin
            temp = var_type( data )
            id = where( fix( missing_replace_list[*,0] ) eq temp, n_id )
            if n_id eq 1 then begin
              if temp eq 7 then begin
                missing_replace = missing_replace_list[id[0],1]
              endif else begin
                missing_replace = float( missing_replace_list[id[0],1] )
              endelse
              check_missing = 1
            endif
          endelse
          ; Proceed only if missing flag attribute and a replacement value exist
          if ( max( missing_label eq var_att ) eq 1 ) $
              and ( check_missing eq 1 ) then begin
            ; Determine missing value flag
            ncdf_attget, id_file, var_label_file[i_var], missing_label, $
                missing_value
            ; Convert missing value flag to standard NaN
            id = where( data eq missing_value, n_id )
            if n_id gt 0 then data[id] = missing_replace
          ; Otherwise warn of failure
          endif else begin
            if !except ne 0 then begin
              print, 'Warning: given missing value name or replacement does ' $
                  + 'not exist'
            endif
          endelse
        endif
        ; Adjust if instructed (scaling and shifting)
        if keyword_set( data ) then begin
          id = where( var_att eq 'scale_factor', n_id )
          if n_id ne 0 then begin
            ncdf_attget, id_file, var_label_file[i_var], 'scale_factor', temp
            data = data * temp
          endif else begin
            id = where( var_att eq 'scaling_factor', n_id )
            if n_id ne 0 then begin
              ncdf_attget, id_file, var_label_file[i_var], 'scaling_factor', $
                  temp
              data = data * temp
            endif
          endelse
          id = where( var_att eq 'add_offset', n_id )
          if n_id ne 0 then begin
            ncdf_attget, id_file, var_label_file[i_var], 'add_offset', temp
            data = data + temp
          endif
        endif
        ; Determine units
        id = where( var_att eq 'units', n_id )
        if n_id ne 0 then begin
          ncdf_attget, id_file, var_label_file[i_var], 'units', units
          units = string( units )
        endif
        ; Find and record requested attributes
        for i_attribute = 0, n_attribute - 1 do begin
          if max( var_att eq attribute_label[i_attribute] ) eq 1 then begin
            ncdf_attget, id_file, var_label_file[i_var], $
                attribute_label[i_attribute], temp
          endif else begin
            temp = ''
          endelse
          attribute_value[i_attribute] = strtrim( string( temp ), 2 )
          attribute_type[i_attribute] = var_type( temp )
          ; Assume that we should not have attributes of type byte
          ; (change to string)
          if attribute_type[i_attribute] eq 1 then begin
            attribute_type[i_attribute] = 7
          endif
        endfor
      endif
    endif
  endif
endfor

;***********************************************************************
; Close the File

; Close the file
ncdf_close, id_file

; Return zero if data not found
if n_elements( data ) eq 0 then data = !values.f_nan

;***********************************************************************
; The End

return, data
END
