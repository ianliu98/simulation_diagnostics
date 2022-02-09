;+
; NAME:
;    markup_read
;
; PURPOSE:
;    This procedure reads entries from an input file following a markup 
;    language format.
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    markup_read, file_name, settings=settings
;
; INPUTS:
;    FILE_NAME:  A required scalar string containing the name of the markup 
;        format file, including the directory.
;    COMMENT_CHAR, SELECT_HEADERS, SELECT_VALUES
;
; KEYWORD PARAMETERS:
;    COMMENT_CHAR:  An optional string containing the first character used in 
;        any comment lines.  The default is '#'.
;    HEADERS:  Returns the list of headers (categories/tags) included in the 
;        markup file.  Of length N_HEADERS.  If SELECT_HEADERS is input then 
;        only the headers listed in SELECT_HEADERS will be returned.
;    SELECT_HEADERS:  An optional string vector containing the list of headers 
;        (categories/tags) for which to return entries in HEADERS and SETTINGS.
;    SELECT_VALUES:  An optional string array containing a list of category/tag 
;        values for which to return the entries.  Elements must be of the form 
;        'header=value'.  For instance, in order to only return entries which 
;        include the '<DOMAIN>global</DOMAIN>' setting, enter 'DOMAIN=global'.
;    SETTINGS:  Returns the list of values for all the (possibly selected) 
;        categories for all of the (possibly selected) entries in the markup 
;        file.  Of size N_HEARERS*N_SETTINGS.
;
; OUTPUTS:
;    HEADERS, SETTINGS
;
; USES:
;    isin.pro
;    sread.pro
;    str.pro
;
; PROCEDURE:
;    This procedure reads values from a nested markup file of the format:
;      <CATEGORY_1>
;        <CATEGORY_2>
;          Value_2
;        </CATEGORY_2>
;        <CATEGORY_3>
;          Value_3a, Value_3b
;        </CATEGORY_3>
;      </CATEGORY_1>
;
; EXAMPLE:
;    See netcdf_read_geo.pro's use of netcdf_read_map_var.txt.
;
; MODIFICATION HISTORY:
;    Written:  Daithi A. Stone (dastone@runbox.com), 2016-10-07  (As 
;        wraf4_settings_read.pro.
;    Modified:  DAS, 2016-11-15 (Switched from '&' string separator to '%'.)
;    Modified:  DAS, 2016-12-06 (Modified SELECT_VALUES implementation to 
;        require satisfaction of all criteria.)
;    Modified:  DAS, 2017-01-25 (Permitted meaningless spaces to be included 
;        around values in the requested and existing settings)
;    Modified:  DAS, 2017-09-25 (Branched from wraf4_settings_read.pro;  
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2017-10-22 (Changed method for selecting for headers and 
;        values such that the input order is retained;  Satisfaction of all 
;        criteria no longer required)
;    Modified:  DAS, 2017-11-08 (Removed stop when requested setting not found)
;    Modified:  DAS, 2017-12-11 (Changed the method for selecting according to 
;        requests in order to worked with nested levels)
;    Modified:  DAS, 2018-08-22 (Removed automatic stop when requested header 
;        not found)
;-

;***********************************************************************

PRO MARKUP_READ, $
    FILE_NAME, $
    COMMENT_CHAR=comment_char, $
    SELECT_HEADERS=select_headers, $
    SELECT_VALUES=select_values, $
    HEADERS=headers, SETTINGS=settings

;***********************************************************************
; Constants

; Abort if no file requested or does not exist
if not( keyword_set( file_name ) ) then begin
  stop, 'ERROR markup_read.pro:  No file name provided.'
endif
temp = file_search( file_name, count=n_temp )
if n_temp ne 1 then begin
  stop, 'ERROR markup_read.pro:  File not found (' + file_name[0] + ').'
endif

; The default comment flag character
if n_elements( comment_char ) eq 0 then comment_char = '#'

;***********************************************************************
; Load the source list

; Read the settings file
list = sread( file_name, no_columns=1, comment_char=comment_char, nocompress=1 )
list = strtrim( list, 2 )
id = where( list ne '', n_list )
if n_list eq 0 then stop
list = list[id]

; Clear any leaked header or settings arrays from the procedure call
headers = ''
settings = ''
; Initialise the count of the types and instances of environments
n_headers = 0
n_settings = 0

; Initialise list of current category tree
current_cat = ''
n_current_cat = 0

; Initialise skip lines flag (no skipping for now).
; Possible values are '' (regular no skipping), or '<HEADER>' or '<COMMENT>', 
; which means all lines are skipped until '</HEADER>' or '</COMMENT>' is 
; encountered respectively.
flag_skip = ''

; Iterate through the list
for i_list = 0l, n_list - 1l do begin
  ; Copy this list entry
  temp_list = list[i_list]
  ; If we are currently waiting to end a sequence of line skipping
  if flag_skip eq '<HEADER>' then begin
    if temp_list eq '</HEADER>' then flag_skip = ''
  endif else if flag_skip eq '<COMMENT>' then begin
    if temp_list eq '</COMMENT>' then flag_skip = ''
  ; If this is the beginning of a sequence of line skipping
  endif else if max( temp_list eq [ '<HEADER>', '<COMMENT>' ] ) eq 1 then begin
    flag_skip = temp_list
  ; If this is the end of a category
  endif else if strmid( temp_list, 0, 2 ) eq '</' then begin
    ; Ensure we are in a category
    if n_current_cat eq 0 then begin
      temp = 'ERROR markup_read.pro:  End of category encountered but we are ' $
          + 'not in a category (line ' + str( i_list ) + ' of file ' $
          + file_name + ').'
      stop, temp
    endif
    ; Ensure this end of category is for the current lowest category
    if '</' + current_cat[n_current_cat-1] + '>' ne temp_list then begin
      temp = 'ERROR markup_read.pro:  End of category encountered but does ' $
          + ' not match current category (line ' + str( i_list ) + ' of file ' $
          + file_name + ').'
      stop, temp
    endif
    ; Move up a category
    n_current_cat = n_current_cat - 1
    current_cat[n_current_cat] = ''
  ; If this is the beginning of a category
  endif else if strmid( temp_list, 0, 1 ) eq '<' then begin
    ; Extract the category label
    temp_list = strmid( temp_list, 1, strlen( temp_list ) - 2 )
    ; Add the category
    if n_elements( current_cat ) ge n_current_cat + 1 then begin
      current_cat[n_current_cat] = temp_list
    endif else if n_current_cat eq 0 then begin
      current_cat = temp_list
    endif else begin
      current_cat = [ current_cat, temp_list ]
    endelse
    n_current_cat = n_current_cat + 1
    ; Copy category label to header, if not already there
    if n_headers eq 0 then begin
      headers = temp_list
      n_headers = n_headers + 1
    endif else if max( temp_list eq headers ) eq 0 then begin
      headers = [ headers, temp_list ]
      n_headers = n_headers + 1
      if n_settings ge 1 then begin
        settings = [ settings, strarr( 1, n_settings ) ]
      endif else begin
        settings = strarr( n_headers )
      endelse
    endif
    ; Determine if a new line is needed for the settings.
    ; If this is the first line
    if n_settings eq 0 then begin
      settings = strarr( n_headers, n_settings + 1 )
      n_settings = n_settings + 1
    ; If this category was preceded by another instance of the the same category
    endif else if list[i_list-1] eq '</' + temp_list + '>' then begin
      settings = [ [ settings ], [ strarr( n_headers, 1 ) ] ]
      if n_current_cat ge 2 then begin
        id = where( headers eq current_cat[n_current_cat-2], n_id )
        if n_id ne 1 then stop
        settings[0:id[0],n_settings] = settings[0:id[0],n_settings-1]
      endif
      n_settings = n_settings + 1
    endif
  ; If we are within a category and this is a data string
  endif else if n_current_cat gt 0 then begin
    ; If this is the first setting then initialise the array
    ;if n_settings eq 0 then settings = strarr( n_headers, 1 )
    id = where( headers eq current_cat[n_current_cat-1], n_id )
    if n_id ne 1 then stop
    settings[id[0],n_settings-1] = temp_list
    ;if strlen( settings[id[0],n_settings-1] ) eq 0 then begin
    ;  settings[id[0],n_settings-1] = temp_list
    ;endif else begin
    ;  settings[id[0],n_settings-1] = settings[id[0],n_settings-1] + '&' $
    ;      + temp_list
    ;endelse
  ; Otherwise ignore this line
  endif
endfor

; Remove empty columns
temp_max = intarr( n_headers )
for i_header = 0, n_headers - 1 do begin
  temp_max[i_header] = max( strlen( settings[i_header,*] ) )
endfor
id = where( temp_max gt 0, n_headers )
if n_headers eq 0 then stop
headers = headers[id]
settings = settings[id,*]

;***********************************************************************
; Select requested information

; Restrict settings array to instances satisfying the requested criteria
n_select_values = n_elements( select_values )
if n_select_values gt 0 then begin
  ; Initialise array flagging entries we are keeping
  flag_settings = intarr( n_settings )
  ; Interate through values selections
  for i_value = 0, n_select_values - 1 do begin
    ; Parse command, of form category=value1,value2,value3...
    temp = strsplit( select_values[i_value], '=', extract=1, count=n_temp )
    if n_temp ne 2 then stop
    temp_header = temp[0]
    temp_values = strsplit( temp[1], ',', extract=1, count=n_temp_values )
    temp_values = strtrim( temp_values, 2 )
    ; Locate column for this category
    id_temp_header = where( headers eq temp_header, n_id )
    if n_id gt 1 then stop
    if n_id eq 1 then begin
      ; Find and copy values for requested values
      if i_value eq 0 then begin
        id_flag = indgen( n_settings )
      endif else begin
        id_flag = where( flag_settings eq 1, n_id_flag )
      endelse
      if n_id eq 0 then stop
      for i_temp_values = 0, n_temp_values - 1 do begin
        id = where( strpos( ',' + settings[id_temp_header[0],id_flag] + ',', $
            ',' + temp_values[i_temp_values] + ',' ) ge 0, n_id )
        if n_id ge 1 then flag_settings[id_flag[id]] = 1
      endfor
    endif
  endfor
  ; Copy new values array for output
  id = where( flag_settings eq 1, n_settings )
  if n_settings gt 0 then begin
    settings = settings[*,id]
    flag_settings = 0
  endif
endif

; Restrict header and settings arrays to requested categories only
n_select_headers = n_elements( select_headers )
if n_select_headers gt 0 then begin
  ; Initialise new output settings array
  settings_new = strarr( n_select_headers, n_settings )
  ; Interate through header selections
  for i_header = 0, n_select_headers - 1 do begin
    ; Find and copy values for requested header
    id = where( headers eq select_headers[i_header], n_id )
    if n_id gt 1 then stop
    if n_id eq 1 then settings_new[i_header,*] = settings[id[0],*]
  endfor
  ; Copy new header and values arrays for output
  headers = select_headers
  settings = settings_new
endif

;***********************************************************************
; The End

return
END
