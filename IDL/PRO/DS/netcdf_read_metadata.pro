;+
; NAME:
;    netcdf_read_metadata.pro
;
; PURPOSE:
;    This function reads and returns the list of variables stored in the 
;    requested NetCDF file.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    result = netcdf_read( file_name )
;
; INPUTS:
;    FILE_NAME:  A required scalar string containing the name of the NetCDF 
;        file.
;
; KEYWORD PARAMETERS:
;    -
;
; OUTPUTS:
;    RESULT:  Returns a list of all of the variables in the NetCDF file.
;
; USES:
;    -
;
; PROCEDURE:
;    This function reads the metadata from a NetCDF file and extracts the list 
;    of variables.
;
; EXAMPLE:
;    result = netcdf_read_metadata( 'file.nc' )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-11-30
;-

;***********************************************************************

FUNCTION NETCDF_READ_METADATA, $
    FILE_NAME

;***********************************************************************
; Constants and variables

; Confirm a legal file name has been requested
if n_elements( file_name ) ne 1 then stop

;***********************************************************************
; Read the metadata from the file

; Open the file
id_file = ncdf_open( file_name[0] )

; Find out what is in the file
info = ncdf_inquire( id_file )
; The number of variables
n_var = info.nvars

; Find the variable labels
var_label = strarr( n_var )
for i_var = 0, n_var - 1 do begin
  var_info = ncdf_varinq( id_file, i_var )
  var_label[i_var] = var_info.name
endfor

; Close the file
ncdf_close, id_file

;***********************************************************************
; The end

return, var_label
END
