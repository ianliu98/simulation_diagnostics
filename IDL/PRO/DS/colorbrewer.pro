;+
; NAME:
;    COLORBREWER
;
; PURPOSE:
;    This procedure implements the colour palettes suggested at 
;    http://colorbrewer2.org.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    colorbrewer, table_label, n_color, dir_table=dir_table
;
; INPUTS:
;    TABLE_LABEL:  The required label of the colour table used by 
;        colorbrewer2.org.  For instance 'RdBu' for the Red-Blue diverging 
;        colour table.  Of type string.
;    N_COLOR:  The required number of colours to include in the table.  If 
;        colorbrewer2.org defines a table with that number then that table is 
;        adopted;  otherwise a table is linearly interpolated from the existing 
;        defined table with the most colours.  Of type integer, and must be 
;        greater than 2.
;    ADD_RGB, DIR_TABLE, REVERSE, START_INDEX
;
; OUTPUT:
;    ---
;
; KEYWORD PARAMETERS:
;    ADD_RGB:  An optional integer array of size [3,N_ADD] providing the RGB 
;        values for N_ADD additional colour to append to the end of the new 
;        colour table.  Values in the first dimension are for the red, green, 
;        and blue colours respectively.
;    DIR_TABLE:  An optional string listing the directory holding the
;        colorbrewer2.org GIMP/Inkscape *.gpl colour tables.  These files are 
;        available as individual files at 
;        https://github.com/axismaps/colorbrewer/tree/master/export/gpl 
;        or via the full colorbrewer2.org downloand at 
;        https://github.com/axismaps/colorbrewer/archive/master.zip as the 
;        export/gpl/ subdirectory.  The default is "colorbrewer2/gpl/" within 
;        the local directory.
;    REVERSE:  If set then the order of colours is reversed, such that for 
;        instance a red-blue table becomes a blue-red one.
;    START_INDEX:  The optional colour index value at which to add the new 
;        colours to the active colour table.
;
; USES:
;    https://github.com/axismaps/colorbrewer/tree/master/export/gpl
;    sread.pro
;
; PROCEDURE:
;    This procedure reads the colour tables from the colorbrewer2.org 
;    GIMP/Inkscape export files provided at colorbrewer2.org..
;
; LICENSE:
;    The colorbrewer2.org license is described at 
;    http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (dastone@runbox.com), 2017-06-09
;    Modified:  DAS, 2018-08-15 (Add the DIR_TABLE keyword input)
;-

;***********************************************************************

PRO COLORBREWER, $
    TABLE_LABEL, N_COLOR, $
    ADD_RGB=add_rgb, $
    DIR_TABLE=dir_table, $
    REVERSE=reverse_opt, $
    START_INDEX=start_index

;***********************************************************************
; Settings and constants

; Ensure table and number of colours defined
if n_elements( table_label ) ne 1 then stop
if n_elements( n_color ) ne 1 then stop
if n_color lt 3 then stop

; The default start index
if n_elements( start_index ) eq 0 then start_index = 0

; Define the directory containing the colour tables
if n_elements( dir_table ) eq 0 then dir_table = 'colorbrewer2/gpl/'
dir_table_len = strlen( dir_table )

;***********************************************************************
; Build the colour table

; Find all versions of the requested table
file_table = file_search( dir_table + table_label + '_*.gpl', $
    count=n_file_table )
file_table_len = strlen( file_table )
file_table_n_color = strarr( n_file_table )
for i_file = 0, n_file_table - 1 do begin
  file_table[i_file] = strmid( file_table[i_file], dir_table_len, $
      file_table_len[i_file]-dir_table_len )
  temp = strsplit( file_table[i_file], '_.', extract=1, count=n_temp )
  if n_temp ne 3 then stop
  file_table_n_color[i_file] = temp[1]
endfor
file_table_n_color = fix( file_table_n_color )

; Find the file with the requested number of colours
id = where( file_table_n_color eq n_color, n_id )
if n_id eq 1 then begin
  file_table = file_table[id[0]]
  file_table_n_color = file_table_n_color[id[0]]
; Otherwise take the file with the largest number of colours
endif else begin
  id = where( file_table_n_color eq max( file_table_n_color ), n_id )
  file_table = file_table[id[0]]
  file_table_n_color = file_table_n_color[id[0]]
endelse

; Read colour indices from file
color_list = sread( dir_table + file_table, separator=' ', skip_lines=4 )
table_red = fix( reform( color_list[0,*] ) )
table_blue = fix( reform( color_list[2,*] ) )
table_green = fix( reform( color_list[1,*] ) )

; Interpolate to a larger number of colours of necessary
if n_color gt file_table_n_color then begin
  x_old = ( findgen( file_table_n_color ) + 0.5 ) / file_table_n_color
  x_new = ( findgen( n_color ) + 0.5 ) / n_color
  table_red = round( interpol( float( table_red ), x_old, x_new ) )
  id = where( table_red lt 0, n_id )
  if n_id gt 0 then table_red[id] = 0
  id = where( table_red gt 255, n_id )
  if n_id gt 0 then table_red[id] = 255
  table_green = round( interpol( float( table_green ), x_old, x_new ) )
  id = where( table_green lt 0, n_id )
  if n_id gt 0 then table_green[id] = 0
  id = where( table_green gt 255, n_id )
  if n_id gt 0 then table_green[id] = 255
  table_blue = round( interpol( float( table_blue ), x_old, x_new ) )
  id = where( table_blue lt 0, n_id )
  if n_id gt 0 then table_blue[id] = 0
  id = where( table_blue gt 255, n_id )
  if n_id gt 0 then table_blue[id] = 255
endif

; If we want to reverse the order of colours
if keyword_set( reverse_opt ) then begin
  table_red = reverse( table_red )
  table_green = reverse( table_green )
  table_blue = reverse( table_blue )
endif

; Add any requested additional colours
if keyword_set( add_rgb ) then begin
  if n_elements( add_rgb[*,0] ) ne 3 then stop
  table_red = [ table_red, reform( add_rgb[0,*] ) ]
  table_green = [ table_green, reform( add_rgb[1,*] ) ]
  table_blue = [ table_blue, reform( add_rgb[2,*] ) ]
endif

; Set the new colour table
tvlct, table_red, table_green, table_blue, start_index

;***********************************************************************
; The end

return
END

