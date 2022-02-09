;+
; NAME:
;    netcdf_write
;
; PURPOSE:
;    This procedure writes data from one or more variables and their associated 
;    dimensions and metadata to a NetCDF file.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    netcdf_write, file_name
;
; INPUTS:
;    FILE_NAME:  An optional scalar string containing the name of the output 
;        NetCDF file.  The default is 'netcdf_write_data.nc'.
;    DATA_ARRAY, DATA_ATTRIBUTE_LABEL, DATA_ATTRIBUTE_VALUE, DATA_DIM, 
;      DATA_LABEL, DIM{1,2,3,4,5}_ATTRIBUTE_LABEL, 
;      DIM{1,2,3,4,5}_ATTRIBUTE_VALUE, DIM{1,2,3,4,5}_DIM_LABEL, 
;      DIM{1,2,3,4,5}_LABEL, DIM{1,2,3,4,5}_VECTOR, DRIVER_NAME, 
;      FILE_NETCDF_READ_GEO_VARINFO, GLOBAL_ATTRIBUTE_LABEL, 
;      GLOBAL_ATTRIBUTE_VALUE, USERNAME=username
;
; KEYWORD PARAMETERS:
;    DATA_ARRAY:  An optional numerical array containing the data for the main 
;        variables to be written in the file.  If input, then it should have 
;        the same dimensions as set in the DIM*_VECTOR lengths, in the order of 
;        the "*" numbers, plus a final dimension containing the number of 
;        variables.  For instance if DIM1_VECTOR is of length 3 and DIM2_VECTOR 
;        is of length 4, DIM3_VECTOR, DIM4_VECTOR, and DIM5_VECTOR are not 
;        input, and there are N_DATA_VAR=2 variables, then DATA_ARRAY should be 
;        of size 3*4*2.
;    DATA_ATTRIBUTE_LABEL:  An optional string array containing the list of 
;        labels of the attributes for the variables in DATA_ARRAY to be 
;        included in the file.  Of size N_DATA_ATTRIBUTE*N_DATA_VAR.  If 
;        DATA_ATTRIBUTE_VALUE is input, then this is required.
;    DATA_ATTRIBUTE_TYPE:  An optional integer array listing the data types for 
;        the N_DATA_ATTRIBUTE*N_DATA_VAR data variable attributes.  For 
;        instance a value of "7" would mean that the variable should be written 
;        as a string.  The default is to assume all values are strings.
;    DATA_ATTRIBUTE_VALUE:  An optional string array containing the list of 
;        values of the attributes for the variables in DATA_ARRAY to be 
;        included in the file.  Of size N_DATA_ATTRIBUTE*N_DATA_VAR.  The 
;        labels of the attributes are given in DATA_ATTRIBUTE_LABEL.  If a 
;        value is blank, then the procedure will adopt a default if available.  
;    DATA_DIM:  An optional vector of size N_DATA_DIM containing the indices 
;        of the dimensions used in DATA_ARRAY.  The default is to assume that 
;        all defined dimensions are used.
;    DATA_LABEL:  A string vector containing the labels for the variables with 
;        data in DATA_ARRAY.  Required if DATA_ARRAY is input.
;    DATA_TYPE:  An optional integer array listing the data types for the 
;        N_DATA_VAR variables.  Of length N_DATA_VAR.  For instance, a value of 
;        "7" would mean that the variable should be saved as a string.  The 
;        default is to use the data type of the DATA_ARRAY input data.
;    DIM{1,2,3,4,5}_ATTRIBUTE_LABEL:  An optional string array containing the 
;        list of labels of the attributes for the variables associated with 
;        dimension {1,2,3,4,5}.  Of size N_DIM{1,2,3,4,5}_ATTRIBUTE
;        *DIM{1,2,3,4,5}_VAR_NUM.  If the respective 
;        DIM{1,2,3,4,5}_ATTRIBUTE_VALUE is input, then this is required.
;    DIM{1,2,3,4,5}_ATTRIBUTE_TYPE:  An optional integer array listing the data 
;        types of the attributes for the variables associated with the 
;        dimension {1,2,3,4,5}.  Of size N_DIM{1,2,3,4,5}_ATTRIBUTE
;        *DIM{1,2,3,4,5}_VAR_NUM.  For instance a value of "7" would mean that 
;        the variable should be written as a string.  The default is to assume 
;        all values are strings.
;    DIM{1,2,3,4,5}_ATTRIBUTE_VALUE:  An optional string array containing the 
;        list of values of the attributes for the variables associated with 
;        dimension {1,2,3,4,5}.  Of size N_DIM{1,2,3,4,5}_ATTRIBUTE
;        *DIM{1,2,3,4,5}_VAR_NUM.  The labels of the attributes are given in 
;        DIM{1,2,3,4,5}_ATTRIBUTE_LABEL.  If a value is blank, then the 
;        procedure will adopt a default if available.
;    DIM{1,2,3,4,5}_DIM_LABEL:  An optional string containing the label for the 
;        dimension {1,2,3,4,5}.  If not input then the first entry in 
;        DIM{1,2,3,4,5}_LABEL is used.
;    DIM{1,2,3,4,5}_LABEL:  A string vector of length DIM{1,2,3,4,5}_VAR_NUM 
;        containing the labels for the variables associated with dimension 
;        {1,2,3,4,5}.  If DIM{1,2,3,4,5}_DIM_LABEL is not input, then this 
;        label for the first variable is also used as the label for the 
;        dimension itself.  Required when the respective DIM{1,2,3,4,5}_VECTOR 
;        is input.
;    DIM{1,2,3,4,5}_TYPE:  An optional integer array listing the data types for 
;        the N_DIM{1,2,3,4,5}_VAR variables associated with dimension 
;        {1,2,3,4,5}.  Of length N_DIM{1,2,3,4,5}_VAR.  For instance, a value 
;        of "7" would mean that the variable should be saved as a string.  The 
;        default is to use the data type of the DIM{1,2,3,4,5}_VECTOR input 
;        data.
;    DIM{1,2,3,4,5}_VECTOR:  An optional array of size 
;        DIM{1,2,3,4,5}_VECTOR_LEN*DIM{1,2,3,4,5}_VAR_NUM containing the data 
;        for the dimension variables associated with dimension {1,2,3,4,5}.  It 
;        can be of type string or a numerical type.
;    DRIVER_NAME:  An optional scalar string containing the name of the program 
;        calling this procedure.  The value will be included in the default 
;        output files history global attribute.  This is ignored if a 'history' 
;        global attribute is defined in GLOBAL_ATTRIBUTE_LABEL.
;    FILE_NETCDF_READ_GEO_VARINFO:  An optional scalar string containing the 
;        name of the file, including directory, containing default metadata.  
;        The default is defined in netcdf_write_metadata.pro.
;    GLOBAL_ATTRIBUTE_LABEL:  An optional string array containing the list of 
;        labels of the global attributes to be included in the file.  Of size 
;        N_GLOBAL_ATTRIBUTE.  If GLOBAL_ATTRIBUTE_VALUE is input, then this is 
;        required.
;    GLOBAL_ATTRIBUTE_TYPE:  An optional integer array listing the data 
;        types of the global attributes.  Of size N_GLOBAL_ATTRIBUTE  For 
;        instance a value of "7" would mean that the variable should be written 
;        as a string.  The default is to assume all values are strings.
;    GLOBAL_ATTRIBUTE_VALUE:  An optional string array containing the list of 
;        values of the global attributes to be included in the file.  Of size 
;        N_GLOBAL_ATTRIBUTE.  The labels of the attributes are given in 
;        GLOBAL_ATTRIBUTE_LABEL.  Values that are to be written as strings 
;        should be enclosed in double quotation marks, i.e. '"value"'.
;    USERNAME:  An optional scalar string containing the name of the user 
;        calling this procedure.  The value will be included in the default 
;        output files history global attribute.  The default is defined by the 
;        $USER environment variable in the operating system.  This is ignored 
;        if a 'history' global attribute is defined in GLOBAL_ATTRIBUTE_LABEL.
;
; OUTPUTS:
;    ---
;
; USES:
;    dimension.pro
;    netcdf_write_metadata.pro
;    var_type.pro
;
; PROCEDURE:
;    This procedure opens a NetCDF file for output, initialises dimensions and 
;    variables in the file, adds metadata, and adds the dimension variables and 
;    data variables.
;
; EXAMPLE:
;    This will write a single longitude vector to the file 'lon.nc'.
;      netcdf_write, 'lon.nc', dim1_vector=findgen(360), dim1_label='lon'
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2016-09-09, as 
;        climate_write.pro.
;    Modified:  DAS, 2017-10-19 (Branched from climate_write.pro;  
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2017-10-27 (Added DATA_TYPE, DIM{1,2,3,4,5}_TYPE, and 
;        USERNAME keyword parameters)
;    Modified:  DAS, 2017-12-23 (Added DATA_ATTRIBUTE_TYPE, 
;        DIM{1,2,3,4,5}_ATTRIBUTE_TYPE, and GLOBAL_ATTRIBUTE_TYPE keyword 
;        inputs).
;    Modified:  DAS, 2018-07-06 (Added DIM{1,2,3,4,5}_DIM_LABEL, DATA_DIM 
;        keyword inputs;  fixed bug with multiple variables per dimension)
;    Modified:  DAS, 2018-08-16 (Fixed issue with writing data variable)
;-

;***********************************************************************

PRO NETCDF_WRITE, $
    FILE_NAME, $
    DATA_ARRAY=data_array, DATA_LABEL=data_label, $
      DATA_ATTRIBUTE_LABEL=data_attribute_label, $
      DATA_ATTRIBUTE_TYPE=data_attribute_type, $
      DATA_ATTRIBUTE_VALUE=data_attribute_value, DATA_TYPE=data_type, $
      DATA_DIM=data_dim, $
    DIM1_VECTOR=dim1_vector, DIM1_LABEL=dim1_label, $
      DIM1_ATTRIBUTE_LABEL=dim1_attribute_label, $
      DIM1_ATTRIBUTE_TYPE=dim1_attribute_type, $
      DIM1_ATTRIBUTE_VALUE=dim1_attribute_value, DIM1_TYPE=dim1_type, $
      DIM1_DIM_LABEL=dim1_dim_label, $
    DIM2_VECTOR=dim2_vector, DIM2_LABEL=dim2_label, $
      DIM2_ATTRIBUTE_LABEL=dim2_attribute_label, $
      DIM2_ATTRIBUTE_TYPE=dim2_attribute_type, $
      DIM2_ATTRIBUTE_VALUE=dim2_attribute_value, DIM2_TYPE=dim2_type, $
      DIM2_DIM_LABEL=dim2_dim_label, $
    DIM3_VECTOR=dim3_vector, DIM3_LABEL=dim3_label, $
      DIM3_ATTRIBUTE_LABEL=dim3_attribute_label, $
      DIM3_ATTRIBUTE_TYPE=dim3_attribute_type, $
      DIM3_ATTRIBUTE_VALUE=dim3_attribute_value, DIM3_TYPE=dim3_type, $
      DIM3_DIM_LABEL=dim3_dim_label, $
    DIM4_VECTOR=dim4_vector, DIM4_LABEL=dim4_label, $
      DIM4_ATTRIBUTE_LABEL=dim4_attribute_label, $
      DIM4_ATTRIBUTE_TYPE=dim4_attribute_type, $
      DIM4_ATTRIBUTE_VALUE=dim4_attribute_value, DIM4_TYPE=dim4_type, $
      DIM4_DIM_LABEL=dim4_dim_label, $
    DIM5_VECTOR=dim5_vector, DIM5_LABEL=dim5_label, $
      DIM5_ATTRIBUTE_LABEL=dim5_attribute_label, $
      DIM5_ATTRIBUTE_TYPE=dim5_attribute_type, $
      DIM5_ATTRIBUTE_VALUE=dim5_attribute_value, DIM5_TYPE=dim5_type, $
      DIM5_DIM_LABEL=dim5_dim_label, $
    GLOBAL_ATTRIBUTE_LABEL=global_attribute_label, $
      GLOBAL_ATTRIBUTE_TYPE=global_attribute_type, $
      GLOBAL_ATTRIBUTE_VALUE=global_attribute_value, $
    DRIVER_NAME=driver_name, $
    FILE_NETCDF_READ_GEO_VARINFO=file_netcdf_read_geo_varinfo, $
    USERNAME=username

;***********************************************************************
; Constants

; Default filename
if not( keyword_set( file_name ) ) then file_name = 'netcdf_write_data.nc'

;; The default missing value
;if n_elements( missing_value ) eq 0 then missing_value = !values.f_nan

; Determine the number and length of dimension variables for each dimension
if n_elements( dim1_vector ) gt 0 then begin
  dim1_vector_len = n_elements( dim1_vector[*,0] )
  dim1_var_num = n_elements( dim1_vector[0,*] )
  if n_elements( dim1_label ) ne dim1_var_num then stop
endif else begin
  dim1_vector_len = 0
  dim1_var_num = 0
endelse
if n_elements( dim2_vector ) gt 0 then begin
  dim2_vector_len = n_elements( dim2_vector[*,0] )
  dim2_var_num = n_elements( dim2_vector[0,*] )
  if n_elements( dim2_label ) ne dim2_var_num then stop
endif else begin
  dim2_vector_len = 0
  dim2_var_num = 0
endelse
if n_elements( dim3_vector ) gt 0 then begin
  dim3_vector_len = n_elements( dim3_vector[*,0] )
  dim3_var_num = n_elements( dim3_vector[0,*] )
  if n_elements( dim3_label ) ne dim3_var_num then stop
endif else begin
  dim3_vector_len = 0
  dim3_var_num = 0
endelse
if n_elements( dim4_vector ) gt 0 then begin
  dim4_vector_len = n_elements( dim4_vector[*,0] )
  dim4_var_num = n_elements( dim4_vector[0,*] )
  if n_elements( dim4_label ) ne dim4_var_num then stop
endif else begin
  dim4_vector_len = 0
  dim4_var_num = 0
endelse
if n_elements( dim5_vector ) gt 0 then begin
  dim5_vector_len = n_elements( dim5_vector[*,0] )
  dim5_var_num = n_elements( dim5_vector[0,*] )
  if n_elements( dim5_label ) ne dim5_var_num then stop
endif else begin
  dim5_vector_len = 0
  dim5_var_num = 0
endelse

; Determine dimension sizes
dim_all_vector_len = [ dim1_vector_len, dim2_vector_len, dim3_vector_len, $
    dim4_vector_len, dim5_vector_len ]
dim_all_var_num = [ dim1_var_num, dim2_var_num, dim3_var_num, dim4_var_num, $
    dim5_var_num ]
n_dim_all = n_elements( dim_all_vector_len )
id_dim_use = where( dim_all_vector_len gt 0, n_dim_use )
if n_dim_use eq 0 then stop
id_dim_var = n_dim_use
; Copy dimension labels
dim_all_label = strarr( n_dim_all, max( dim_all_var_num ) )
if dim1_var_num gt 0 then dim_all_label[0,0:dim1_var_num-1] = dim1_label
if dim2_var_num gt 0 then dim_all_label[1,0:dim2_var_num-1] = dim2_label
if dim3_var_num gt 0 then dim_all_label[2,0:dim3_var_num-1] = dim3_label
if dim4_var_num gt 0 then dim_all_label[3,0:dim4_var_num-1] = dim4_label
if dim5_var_num gt 0 then dim_all_label[4,0:dim5_var_num-1] = dim5_label

; Adopt default dimension labels
if not( keyword_set( dim1_dim_label ) ) and keyword_set( dim1_label ) then begin
  dim1_dim_label = dim1_label[0]
endif
if not( keyword_set( dim2_dim_label ) ) and keyword_set( dim2_label ) then begin
  dim2_dim_label = dim2_label[0]
endif
if not( keyword_set( dim3_dim_label ) ) and keyword_set( dim3_label ) then begin
  dim3_dim_label = dim3_label[0]
endif
if not( keyword_set( dim4_dim_label ) ) and keyword_set( dim4_label ) then begin
  dim4_dim_label = dim4_label[0]
endif
if not( keyword_set( dim5_dim_label ) ) and keyword_set( dim5_label ) then begin
  dim5_dim_label = dim5_label[0]
endif
dim_all_dim_label = strarr( n_dim_all )
if keyword_set( dim1_dim_label ) then dim_all_dim_label[0] = dim1_dim_label
if keyword_set( dim2_dim_label ) then dim_all_dim_label[1] = dim2_dim_label
if keyword_set( dim3_dim_label ) then dim_all_dim_label[2] = dim3_dim_label
if keyword_set( dim4_dim_label ) then dim_all_dim_label[3] = dim4_dim_label
if keyword_set( dim5_dim_label ) then dim_all_dim_label[4] = dim5_dim_label

; Determine if we have no data array
if n_elements( data_array ) eq 0 then begin
  n_data_var = 0
; Otherwise determine size of data array
endif else begin
  data_array_size = size( data_array )
  n_data_array_size = dimension( data_array )
  data_array_size = data_array_size[1:n_data_array_size]
  if n_data_array_size ne id_dim_var + 1 then begin
    data_array_size = [ data_array_size, $
        1 + intarr( id_dim_var + 1 - n_data_array_size ) ]
  endif
  ; The number of variables
  n_data_var = data_array_size[id_dim_var]
  ; The indices of dimensions used
  if n_elements( data_dim ) gt 0 then begin
    id_dim_use_data = data_dim
  endif else begin
    id_dim_use_data = id_dim_use
  endelse
  n_dim_use_data = n_elements( id_dim_use_data )
  ; Check data array satisfies dimensions
  if n_data_array_size gt id_dim_var + 1 then stop
  for i_dim = 0, n_dim_use_data - 1 do begin
    if data_array_size[i_dim] ne dim_all_vector_len[id_dim_use_data[i_dim]] $
        then stop
  endfor
endelse

; Check data attribute variables
n_data_attribute = n_elements( data_attribute_label )
if n_data_attribute gt 0 then begin
  n_data_attribute = n_elements( data_attribute_label[*,0] )
  if n_elements( data_attribute_label[0,*] ) ne n_data_var then stop
  if n_elements( data_attribute_value ) gt 0 then begin
    if n_elements( data_attribute_value[0,*] ) ne n_data_var then stop
    if n_elements( data_attribute_value[*,0] ) ne n_data_attribute then stop
  endif
  if n_elements( data_attribute_type ) gt 0 then begin
    if n_elements( data_attribute_type[0,*] ) ne n_data_var then stop
    if n_elements( data_attribute_type[*,0] ) ne n_data_attribute then stop
  endif else begin
    ; Define the default attribute type (string)
    data_attribute_type = 7 + intarr( n_data_attribute, n_data_var )
  endelse
endif

; Check global attribute variables for consistency
n_global_attribute = n_elements( global_attribute_label )
if n_elements( global_attribute_value ) ne n_global_attribute then stop
; Adopt dummy defaults if necessary
if n_global_attribute eq 0 then begin
  global_attribute_label = ''
  global_attribute_value = ''
  n_global_attribute = 1
endif
; Ensure consistent global attribute type vector
if n_elements( global_attribute_type ) gt 0 then begin
  if n_elements( global_attribute_type ) ne n_global_attribute then stop
; Or define the default
endif else begin
  global_attribute_type = 7 + intarr( n_global_attribute )
endelse

; Define default history global attribute addition
if max( global_attribute_label eq 'history' ) eq 0 then begin
  global_attribute_label = [ global_attribute_label, 'history' ]
  if not( keyword_set( username ) ) then begin
    spawn, 'echo $USER', username
  endif
  temp_value = 'Created by ' + username + ' on ' + systime() + ' using '
  if keyword_set( driver_name ) then begin
    if n_elements( driver_name ) ne 1 then stop
    temp_value = temp_value + driver_name + ' and '
  endif
  temp_value = temp_value + 'netcdf_write.pro.'
  global_attribute_value = [ global_attribute_value, temp_value ]
  n_global_attribute = n_global_attribute + 1
  global_attribute_type = [ global_attribute_type, 7 ]
endif

; Adopt default data attributes if necessary and available
if n_data_attribute gt 0 then begin
  netcdf_write_metadata, data_label, attribute_label=data_attribute_label, $
      attribute_value=data_attribute_value, $
      attribute_type=data_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif
; Adopt default dimension variable attributes if necessary and available
if n_elements( dim1_attribute_label ) gt 0 then begin
  netcdf_write_metadata, dim1_label, attribute_label=dim1_attribute_label, $
      attribute_value=dim1_attribute_value, $
      attribute_type=dim1_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif
if n_elements( dim2_attribute_label ) gt 0 then begin
  netcdf_write_metadata, dim2_label, attribute_label=dim2_attribute_label, $
      attribute_value=dim2_attribute_value, $
      attribute_type=dim2_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif
if n_elements( dim3_attribute_label ) gt 0 then begin
  netcdf_write_metadata, dim3_label, attribute_label=dim3_attribute_label, $
      attribute_value=dim3_attribute_value, $
      attribute_type=dim3_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif
if n_elements( dim4_attribute_label ) gt 0 then begin
  netcdf_write_metadata, dim4_label, attribute_label=dim4_attribute_label, $
      attribute_value=dim4_attribute_value, $
      attribute_type=dim4_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif
if n_elements( dim5_attribute_label ) gt 0 then begin
  netcdf_write_metadata, dim5_label, attribute_label=dim5_attribute_label, $
      attribute_value=dim5_attribute_value, $
      attribute_type=dim5_attribute_type, $
      file_netcdf_read_geo_varinfo=file_netcdf_read_geo_varinfo
endif

; Determine data variable types
if n_data_var gt 0 then begin
  if n_elements( data_type ) ne n_data_var then begin
    if n_elements( data_type ) ne 0 then stop
    data_type = -1 + intarr( n_data_var )
  endif
  id = where( data_type eq -1, n_id )
  if n_id gt 0 then data_type[id] = var_type( data_array )
endif
; Determine dimension variable types
if dim1_var_num gt 0 then begin
  if n_elements( dim1_type ) ne dim1_var_num then begin
    if n_elements( dim1_type ) ne 0 then stop
    dim1_type = -1 + intarr( dim1_var_num )
  endif
  id = where( dim1_type eq -1, n_id )
  if n_id gt 0 then dim1_type[id] = var_type( dim1_vector )
endif
; Determine dimension variable types
if dim2_var_num gt 0 then begin
  if n_elements( dim2_type ) ne dim2_var_num then begin
    if n_elements( dim2_type ) ne 0 then stop
    dim2_type = -1 + intarr( dim2_var_num )
  endif
  id = where( dim2_type eq -1, n_id )
  if n_id gt 0 then dim2_type[id] = var_type( dim2_vector )
endif
; Determine dimension variable types
if dim3_var_num gt 0 then begin
  if n_elements( dim3_type ) ne dim3_var_num then begin
    if n_elements( dim3_type ) ne 0 then stop
    dim3_type = -1 + intarr( dim3_var_num )
  endif
  id = where( dim3_type eq -1, n_id )
  if n_id gt 0 then dim3_type[id] = var_type( dim3_vector )
endif
; Determine dimension variable types
if dim4_var_num gt 0 then begin
  if n_elements( dim4_type ) ne dim4_var_num then begin
    if n_elements( dim4_type ) ne 0 then stop
    dim4_type = -1 + intarr( dim4_var_num )
  endif
  id = where( dim4_type eq -1, n_id )
  if n_id gt 0 then dim4_type[id] = var_type( dim4_vector )
endif
; Determine dimension variable types
if dim5_var_num gt 0 then begin
  if n_elements( dim5_type ) ne dim5_var_num then begin
    if n_elements( dim5_type ) ne 0 then stop
    dim5_type = -1 + intarr( dim5_var_num )
  endif
  id = where( dim5_type eq -1, n_id )
  if n_id gt 0 then dim5_type[id] = var_type( dim5_vector )
endif

;***********************************************************************
; Create NetCDF File

; Open file for writing
id_file = ncdf_create( file_name, clobber=1 )

; Put the file into define mode
ncdf_control, id_file, fill=1

; Initialise array of identifiers for the dimensions and string lengths for the 
; dimension within the output file
forfile_dim_use_index = lonarr( n_dim_use )
forfile_dim_use_strlen = -1 $
    + lonarr( n_dim_use, max( dim_all_var_num[id_dim_use] ) )
forfile_dim_use_type = -1 $
    + intarr( n_dim_use, max( dim_all_var_num[id_dim_use] ) )
; Define the dimensions to be included in the file
for i_dim = 0, n_dim_use - 1 do begin
  ; Copy dimension data
  if id_dim_use[i_dim] eq 0 then begin
    temp_dim_vector = dim1_vector
    temp_dim_type = dim1_type
  endif else if id_dim_use[i_dim] eq 1 then begin
    temp_dim_vector = dim2_vector
    temp_dim_type = dim2_type
  endif else if id_dim_use[i_dim] eq 2 then begin
    temp_dim_vector = dim3_vector
    temp_dim_type = dim3_type
  endif else if id_dim_use[i_dim] eq 3 then begin
    temp_dim_vector = dim4_vector
    temp_dim_type = dim4_type
  endif else if id_dim_use[i_dim] eq 4 then begin
    temp_dim_vector = dim5_vector
    temp_dim_type = dim5_type
  endif else begin
    stop
  endelse
  ; Initialise the dimension in the file
  forfile_dim_use_index[i_dim] = ncdf_dimdef( id_file, $
      dim_all_dim_label[id_dim_use[i_dim]], $
      dim_all_vector_len[id_dim_use[i_dim],0] )
  ; Copy the variable types
  forfile_dim_use_type[i_dim,0:dim_all_var_num[id_dim_use[i_dim]]-1] $
      = temp_dim_type
  ; Iterate through variables for this dimension
  for i_var = 0, dim_all_var_num[id_dim_use[i_dim]] - 1 do begin
    ; If the dimension contains string data then determine the maximum string 
    ; length
    if temp_dim_type[i_var] eq 7 then begin
      temp = max( strlen( strtrim( string( temp_dim_vector[*,i_var] ), 2 ) ) )
      forfile_dim_use_strlen[i_dim,i_var] = ncdf_dimdef( id_file, $
          'string_length_' + dim_all_label[id_dim_use[i_dim],i_var], temp )
    endif
  endfor
endfor

; Define the variables to be put into the file.
; Initialise vector of dimension variable indices within the output file
forfile_dim_use_var_index = lonarr( n_dim_use, $
    max( dim_all_var_num[id_dim_use] ) )
; Iterate through dimensions
for i_dim = 0, n_dim_use - 1 do begin
  ; Iterate through variables associated with this dimension
  for i_var = 0, dim_all_var_num[id_dim_use[i_dim]] - 1 do begin
    ; Determine the variable type
    temp_type_short_opt = 0
    temp_type_long_opt = 0
    temp_type_float_opt = 0
    temp_type_double_opt = 0
    temp_type_string_opt = 0
    if forfile_dim_use_type[i_dim,i_var] eq 2 then begin
      temp_type_short_opt = 1
    endif else if forfile_dim_use_type[i_dim,i_var] eq 3 then begin
      temp_type_long_opt = 1
    endif else if forfile_dim_use_type[i_dim,i_var] eq 4 then begin
      temp_type_float_opt = 1
    endif else if forfile_dim_use_type[i_dim,i_var] eq 5 then begin
      temp_type_double_opt = 1
    endif else if forfile_dim_use_type[i_dim,i_var] eq 7 then begin
      temp_type_string_opt = 1
    endif else begin
      stop
    endelse
    ; If this variable is of a numerical type
    if forfile_dim_use_type[i_dim,i_var] ne 7 then begin
      ; Define the index in the file and link to the dimension
      forfile_dim_use_var_index[i_dim,i_var] = ncdf_vardef( id_file, $
          dim_all_label[id_dim_use[i_dim],i_var], $
          [forfile_dim_use_index[i_dim]], short=temp_type_short_opt, $
          long=temp_type_long_opt, float=temp_type_float_opt, $
          double=temp_type_double_opt )
    ; If this variable is of a string type
    endif else begin
      ; Define the index in the file and link to the dimension
      temp = [ forfile_dim_use_strlen[i_dim,i_var], $
          forfile_dim_use_index[i_dim] ]
      forfile_dim_use_var_index[i_dim,i_var] = ncdf_vardef( id_file, $
          dim_all_label[id_dim_use[i_dim],i_var], temp, char=1 )
    endelse
  endfor
endfor
; Define the index of the data arrays
if n_data_var gt 0 then begin
  forfile_data_var_index = lonarr( n_data_var )
  for i_var = 0, n_data_var - 1 do begin
    ; Determine the variable type
    temp_type_short_opt = 0
    temp_type_long_opt = 0
    temp_type_float_opt = 0
    temp_type_double_opt = 0
    temp_type_string_opt = 0
    if data_type[i_var] eq 2 then begin
      temp_type_short_opt = 1
    endif else if data_type[i_var] eq 3 then begin
      temp_type_long_opt = 1
    endif else if data_type[i_var] eq 4 then begin
      temp_type_float_opt = 1
    endif else if data_type[i_var] eq 5 then begin
      temp_type_double_opt = 1
    endif else if data_type[i_var] eq 7 then begin
      temp_type_string_opt = 1
    endif else begin
      stop
    endelse
    forfile_data_var_index[i_var] = ncdf_vardef( id_file, data_label[i_var], $
        [forfile_dim_use_index[id_dim_use_data]], short=temp_type_short_opt, $
        long=temp_type_long_opt, float=temp_type_float_opt, $
        double=temp_type_double_opt, char=temp_type_char_opt )
  endfor
endif

; Add dimension attributes to file.
; Iterate through dimensions
for i_dim = 0, n_dim_use - 1 do begin
  ; Determine the available attributes for this dimension
  temp_attribute_label = 0
  if id_dim_use[i_dim] eq 0 then begin
    if n_elements( dim1_attribute_label ) gt 0 then begin
      temp_attribute_label = dim1_attribute_label
      temp_attribute_value = dim1_attribute_value
      temp_attribute_type = dim1_attribute_type
    endif
  endif else if id_dim_use[i_dim] eq 1 then begin 
    if n_elements( dim2_attribute_label ) gt 0 then begin
      temp_attribute_label = dim2_attribute_label
      temp_attribute_value = dim2_attribute_value
      temp_attribute_type = dim2_attribute_type
    endif
  endif else if id_dim_use[i_dim] eq 2 then begin 
    if n_elements( dim3_attribute_label ) gt 0 then begin
      temp_attribute_label = dim3_attribute_label
      temp_attribute_value = dim3_attribute_value
      temp_attribute_type = dim3_attribute_type
    endif
  endif else if id_dim_use[i_dim] eq 3 then begin 
    if n_elements( dim4_attribute_label ) gt 0 then begin
      temp_attribute_label = dim4_attribute_label
      temp_attribute_value = dim4_attribute_value
      temp_attribute_type = dim4_attribute_type
    endif
  endif else if id_dim_use[i_dim] eq 4 then begin 
    if n_elements( dim5_attribute_label ) gt 0 then begin
      temp_attribute_label = dim5_attribute_label
      temp_attribute_value = dim5_attribute_value
      temp_attribute_type = dim5_attribute_type
    endif
  endif else begin
    stop
  endelse
  ; Iterate through variables associated with this dimension
  for i_var = 0, dim_all_var_num[id_dim_use[i_dim]] - 1 do begin
    ; Copy attributes for this variable
    if keyword_set( temp_attribute_label ) then begin
      id_attribute = where( temp_attribute_label[*,i_var] $
          + temp_attribute_value[*,i_var] ne '', n_id_attribute )
    endif else begin
      n_id_attribute = 0
    endelse
    ; Iterate through attributes
    for i_attribute = 0, n_id_attribute - 1 do begin
      ; Determine the variable type for the attribute value
      temp_type_short_opt = 0
      temp_type_long_opt = 0
      temp_type_float_opt = 0
      temp_type_double_opt = 0
      temp_type_string_opt = 0
      temp_type = temp_attribute_type[id_attribute[i_attribute],i_var]
      temp_value = temp_attribute_value[id_attribute[i_attribute],i_var]
      if temp_type eq 2 then begin
        temp_type_short_opt = 1
        temp_value = fix( temp_value )
      endif else if temp_type eq 3 then begin
        temp_type_long_opt = 1
        temp_value = long( temp_value )
      endif else if temp_type eq 4 then begin
        temp_type_float_opt = 1
        temp_value = float( temp_value )
      endif else if temp_type eq 5 then begin
        temp_type_double_opt = 1
        temp_value = double( temp_value )
      endif else if temp_type eq 7 then begin
        temp_type_string_opt = 1
      endif else begin
        stop
      endelse
      ; Add attributes to file
      ncdf_attput, id_file, forfile_dim_use_var_index[i_dim,i_var], $
          temp_attribute_label[id_attribute[i_attribute],i_var], $
          short=temp_type_short_opt, long=temp_type_long_opt, $
          float=temp_type_float_opt, double=temp_type_double_opt, $
          char=temp_type_string_opt, temp_value
    endfor
  endfor
endfor

; Define data variable attributes.
; Iterate through variables
for i_var = 0, n_data_var - 1 do begin
  ; Iterate through attributes
  for i_attribute = 0, n_data_attribute - 1 do begin
    ; Extract attribute label and value if it exists
    temp_label = data_attribute_label[i_attribute,i_var]
    temp_value = data_attribute_value[i_attribute,i_var]
    temp_type = data_attribute_type[i_attribute,i_var]
    if temp_label + temp_value ne '' then begin
      ; Determine the variable type for the attribute value
      temp_type_short_opt = 0
      temp_type_long_opt = 0
      temp_type_float_opt = 0
      temp_type_double_opt = 0
      temp_type_string_opt = 0
      if temp_type eq 2 then begin
        temp_type_short_opt = 1
        temp_value = fix( temp_value )
      endif else if temp_type eq 3 then begin
        temp_type_long_opt = 1
        temp_value = long( temp_value )
      endif else if temp_type eq 4 then begin
        temp_type_float_opt = 1
        temp_value = float( temp_value )
      endif else if temp_type eq 5 then begin
        temp_type_double_opt = 1
        temp_value = double( temp_value )
      endif else if temp_type eq 7 then begin
        temp_type_string_opt = 1
      endif else begin
        stop
      endelse
      ; Add attribute value to file
      ncdf_attput, id_file, forfile_data_var_index[i_var], temp_label, $
          short=temp_type_short_opt, long=temp_type_long_opt, $
          float=temp_type_float_opt, double=temp_type_double_opt, $
          char=temp_type_string_opt, temp_value
    endif
  endfor
endfor

; Define global attributes.
; Iterate through global attributes
for i_attribute = 0, n_global_attribute - 1 do begin
  ; Add global attribute value if it exists
  if global_attribute_label[i_attribute] + global_attribute_value[i_attribute] $
      ne '' then begin
    ; Determine the variable type for the attribute value
    temp_type_short_opt = 0
    temp_type_long_opt = 0
    temp_type_float_opt = 0
    temp_type_double_opt = 0
    temp_type_string_opt = 0
    temp_value = global_attribute_value[i_attribute]
    temp_type = global_attribute_type[i_attribute]
    if temp_type eq 2 then begin
      temp_type_short_opt = 1
      temp_value = fix( temp_value )
    endif else if temp_type eq 3 then begin
      temp_type_long_opt = 1
      temp_value = long( temp_value )
    endif else if temp_type eq 4 then begin
      temp_type_float_opt = 1
      temp_value = float( temp_value )
    endif else if temp_type eq 5 then begin
      temp_type_double_opt = 1
      temp_value = double( temp_value )
    endif else if temp_type eq 7 then begin
      temp_type_string_opt = 1
    endif else begin
      stop
    endelse
    ; Add attribute value to file
    ncdf_attput, id_file, global=1, global_attribute_label[i_attribute], $
        short=temp_type_short_opt, long=temp_type_long_opt, $
        float=temp_type_float_opt, double=temp_type_double_opt, $
        char=temp_type_string_opt, temp_value
  endif
endfor

; Put the file into data mode
ncdf_control, id_file, endef=1

; Write dimension variables to file.
; Iterate through dimensions and dimension variables
for i_dim = 0, n_dim_use - 1 do begin
  for i_var = 0, dim_all_var_num[id_dim_use[i_dim]] - 1 do begin
    ; Get data vector
    if id_dim_use[i_dim] eq 0 then begin
      temp_dim_vector = dim1_vector[*,i_var]
      temp_dim_type = dim1_type[i_var]
    endif else if id_dim_use[i_dim] eq 1 then begin
      temp_dim_vector = dim2_vector[*,i_var]
      temp_dim_type = dim2_type[i_var]
    endif else if id_dim_use[i_dim] eq 2 then begin
      temp_dim_vector = dim3_vector[*,i_var]
      temp_dim_type = dim3_type[i_var]
    endif else if id_dim_use[i_dim] eq 3 then begin
      temp_dim_vector = dim4_vector[*,i_var]
      temp_dim_type = dim4_type[i_var]
    endif else if id_dim_use[i_dim] eq 4 then begin
      temp_dim_vector = dim5_vector[*,i_var]
      temp_dim_type = dim5_type[i_var]
    endif else begin
      stop
    endelse
    ; Ensure data type
    if temp_dim_type eq 2 then begin
      temp_dim_vector = fix( temp_dim_vector )
    endif else if temp_dim_type eq 3 then begin
      temp_dim_vector = long( temp_dim_vector )
    endif else if temp_dim_type eq 4 then begin
      temp_dim_vector = float( temp_dim_vector )
    endif else if temp_dim_type eq 5 then begin
      temp_dim_vector = double( temp_dim_vector )
    endif else if temp_dim_type eq 7 then begin
      temp_dim_vector = strtrim( string( temp_dim_vector ), 2 )
    endif else begin
      stop
    endelse
    ; Write the data to file
    ncdf_varput, id_file, forfile_dim_use_var_index[i_dim,i_var], $
        temp_dim_vector
  endfor
endfor

; And variable data to file.
if n_data_var gt 0 then begin
  ; Modfity data array to convenient format
  temp_dim_all_vector_len = dim_all_vector_len[id_dim_use_data]
  if n_elements( temp_dim_all_vector_len ) lt n_dim_all then begin
    temp_dim_all_vector_len = [ temp_dim_all_vector_len, $
        1+ intarr( n_dim_all - n_elements( temp_dim_all_vector_len ) ) ]
  endif
  data_array_use = reform( data_array, [ temp_dim_all_vector_len, n_data_var ] )
  ; Iterate through variables
  for i_var = 0, n_data_var - 1 do begin
    ; Ensure data type
    if data_type[i_var] eq 2 then begin
      temp_data_array = fix( data_array_use[*,*,*,*,*,i_var] )
    endif else if data_type[i_var] eq 3 then begin
      temp_data_array = long( data_array_use[*,*,*,*,*,i_var] )
    endif else if data_type[i_var] eq 4 then begin
      temp_data_array = float( data_array_use[*,*,*,*,*,i_var] )
    endif else if data_type[i_var] eq 5 then begin
      temp_data_array = double( data_array_use[*,*,*,*,*,i_var] )
    endif else if data_type[i_var] eq 7 then begin
      temp_data_array = strtrim( string( data_array_use[*,*,*,*,*,i_var] ), 2 )
    endif else begin
      stop
    endelse
    ; Convert back to desired format
    temp_data_array = reform( temp_data_array, $
        dim_all_vector_len[id_dim_use_data] )
    ; Write data array to file
    ncdf_varput, id_file, forfile_data_var_index[i_var], temp_data_array
    temp_data_array = 0
  endfor
  data_array_use = 0
endif

; Close the file
ncdf_close, id_file

;***********************************************************************
; The End

return
END
