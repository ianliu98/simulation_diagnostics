;+
; NAME:
;    netcdf_write_metadata.pro
;
; PURPOSE:
;    This procedure returns default NetCDF attribute values for a requested 
;    list of variables.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    netcdf_write_metadata, var_label, attribute_label=attribute_label, attribute_value=attribute_value
;
; INPUTS:
;    VAR_LABEL:  The required string scalar or vector listing the variables for 
;        which to return attribute values.  Of length N_VAR.
;    ATTRIBUTE_LABEL, ATTRIBUTE_VALUE
;
; KEYWORD PARAMETERS:
;    ATTRIBUTE_LABEL:  A required string scalar or array listing the attributes 
;        for the variables in VAR_LABEL.  Of size N_ATTRIBUTE*N_VAR, where 
;        N_ATTRIBUTE is the number of attributes.
;    ATTRIBUTE_VALUE:  An optional string scalar or array containing the values 
;        of the attributes listed in ATTRIBUTE_LABEL.  This procedure will only 
;        fill in empty values (if defaults exist), it will not overwrite 
;        existing values.  This keyword parameter is also a required output, 
;        returning the full list of attribute values.  Of size 
;        N_ATTRIBUTE*N_VAR.
;    ATTRIBUTE_TYPE:  An optional integer scalar or array listing the data type 
;        for the attribute values in ATTRIBUTE_VALUE.  For instance a value of 
;        "7" defines the value as a string.  The default is to assume all 
;        values are strings.
;    FILE_NETCDF_READ_GEO_VARINFO:  An optional string containing the name of 
;        the markup format file containing default NetCDF attribute values.  
;        The default is "netcdf_read_geo_varinfo.xml".
;
; OUTPUTS:
;    ATTRIBUTE_VALUE
;
; USES:
;    markup_read.pro
;    netcdf_read_geo_varinfo.xml
;
; PROCEDURE:
;    This procedure searches for default metadata values in a markup-formatted 
;    library file.
;
; EXAMPLE:
;    netcdf_write_metadata, 'pr', attribute_label=['units','standard_name'], attribute_value=attribute_value
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-10-24 (Added to IDL 
;        routine library)
;    Modified:  DAS, 2017-12-24 (Added ATTRIBUTE_TYPE keyword)
;-

;***********************************************************************

PRO NETCDF_WRITE_METADATA, $
    VAR_LABEL, $
    ATTRIBUTE_LABEL=attribute_label, ATTRIBUTE_TYPE=attribute_type, $
      ATTRIBUTE_VALUE=attribute_value, $
    FILE_NETCDF_READ_GEO_VARINFO=file_netcdf_read_geo_varinfo

;***********************************************************************
; Constants

; Define the file containing standard metadata values
n_file_netcdf_read_geo_varinfo = n_elements( file_netcdf_read_geo_varinfo )
if n_file_netcdf_read_geo_varinfo gt 1 then stop
if n_file_netcdf_read_geo_varinfo eq 0 then begin
  file_netcdf_read_geo_varinfo = file_which( 'netcdf_read_geo_varinfo.xml' )
endif

; Confirm there is a label for this variable
n_var = n_elements( var_label )
if n_elements( var_label ) eq 0 then stop

; Confirm there a consistent attribute label array
if n_elements( attribute_label ) eq 0 then stop
if n_elements( attribute_label[0,*] ) ne n_var then stop
n_attribute = n_elements( attribute_label[*,0] )

; Confirm there is a consistent attribute value array
if n_elements( attribute_value ) gt 0 then begin
  if n_elements( attribute_value[*,0] ) ne n_attribute then stop
  if n_elements( attribute_value[0,*] ) ne n_var then stop
; Or initialise an empty one if none was input
endif else begin
  attribute_value = strarr( n_attribute, n_var )
endelse

; Confirm there is a consistent attribute type array
if n_elements( attribute_type ) gt 0 then begin
  if n_elements( attribute_type[*,0] ) ne n_attribute then stop
  if n_elements( attribute_type[0,*] ) ne n_var then stop
; Or initialise an empty one if none was input
endif else begin
  attribute_type = -1 + intarr( n_attribute, n_var )
  id = where( attribute_value ne '', n_id )
  if n_id gt 0 then attribute_type[id] = 7
endelse

;***********************************************************************
; Load metadata

; Iterate through variables
for i_var = 0, n_var - 1 do begin
  ; Identify these unfilled attribute entries
  id_attribute = where( attribute_value[*,i_var] eq '', n_id_attribute )
  ; Fill empty attribute values with defaults
  if n_id_attribute gt 0 then begin
    markup_read, file_netcdf_read_geo_varinfo, comment_char=';', $
        select_headers=attribute_label[id_attribute,i_var], $
        select_values='label='+var_label[i_var], settings=temp_metadata
    attribute_value[id_attribute,i_var] = temp_metadata
    attribute_type[id_attribute,i_var] = 7
  endif
endfor

;***********************************************************************
; The end

return
END
