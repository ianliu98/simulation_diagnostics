;+
; NAME:
;    c20c_dtos_v2_unnan
;
; PURPOSE:
;    This procedure fills in NaNs in a 2-dimensional spatial field, contained
;    in a 2-dimensional or 3-dimensional array, with the value of the spatially 
;    nearest non-NaN neighbours.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_unnan, data
;
; INPUTS:
;    DATA:  A required input floating point array of size 
;        N_X_DATA,N_Y_DATA,N_Z_DATA, where N_X_DATA is the size of the first 
;        spatial dimension (e.g. longitude), N_Y_DATA is the size of the second 
;        spatial dimension (e.g. latitude), and N_Z_DATA is the size of any 
;        non-spatial dimension (e.g. time).  Also returns the modified array.
;
; KEYWORD PARAMETERS:
;    WRAP_X:  If set then the X dimensional is considered to be circular, i.e. 
;        that the ends are in contact (e.g. longitude on a global polar grid).
;        The default is to be not set.
;    WRAP_Y:  If set then the Y dimensional is considered to be circular, i.e. 
;        that the ends are in contact.  The default is to be not set.
;
; OUTPUTS:
;    DATA
;
; USES:
;    dimension.pro
;
; PROCEDURE:
;    This procedure finds NaN values and samples ever larger squares around 
;    each NaN value until it finds non-NaN values, at which time it replaces 
;    the NaN value with the average of the non-NaN values.
;
; EXAMPLE:
;    ; Create an array of increasing numbers, removing the [2,2] value (12).
;      data = findgen( 5, 5 )
;      data[2,2] = !values.f_nan
     ; Substitute for the missing value
;      c20c_dtos_v2_unnan, data
;    ; The [2,2] value should be 12, in this case recovering the removed value.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2012-07-31, as 
;        c20c_unnan.pro.
;    Modified:  DAS, 2017-09-20 (Branched to c20c_dtos_v2_unnan.pro;  
;        standardised documentation)
;-

;***********************************************************************

PRO C20C_DTOS_V2_UNNAN, $
    DATA, $
    WRAP_X=wrap_x_opt, WRAP_Y=wrap_y_opt

;***********************************************************************
; Constants

; Confirm acceptable dimensions
if not( keyword_set( data ) ) then stop
temp = dimension( data )
if ( temp ne 2 ) and ( temp ne 3 ) then stop

; Confirm that data has some finite values
if max( finite( data ) ) eq 0 then stop

; The data array dimensions
n_x_data = n_elements( data[*,0,0] )
n_y_data = n_elements( data[0,*,0] )
n_z_data = n_elements( data[0,0,*] )

; Determine the wrapping options
wrap_x_opt = keyword_set( wrap_x_opt )
wrap_y_opt = keyword_set( wrap_y_opt )

;***********************************************************************
; Fill in the NaNs

; Iterate through the third dimension
for i_z_data = 0, n_z_data - 1 do begin
  ; Copy the data
  temp_data = data[*,*,i_z_data]
  ; Find values masked to NaN
  id_nan = where( finite( temp_data ) eq 0, n_id_nan )
  ; Determine coordinates in grid measure
  y_id_nan = id_nan / n_x_data
  x_id_nan = id_nan - y_id_nan * n_x_data  
  ; Iterate through missing values
  for i_nan = 0, n_id_nan - 1 do begin
    ; Initialise check and square size
    flag_check = 0
    square = 1
    ; Iterate through wider squares around value
    while flag_check eq 0 do begin
      ; Find points on square
      y_bottom = y_id_nan[i_nan] - square
      y_top = y_id_nan[i_nan] + square
      if y_id_nan[i_nan] - square lt 0 then begin
        if wrap_y_opt eq 1 then begin
          y_botton = y_bottom + n_y_data
        endif else begin
          y_bottom = 0
        endelse
      endif else if y_id_nan[i_nan] + square ge n_y_data then begin
        if wrap_y_opt eq 1 then begin
          y_top = y_top - n_y_data
        endif else begin
          y_top = n_y_data - 1
        endelse
      endif
      x_left = x_id_nan[i_nan] - square
      x_right = x_id_nan[i_nan] + square
      if x_id_nan[i_nan] - square lt 0 then begin
        if wrap_x_opt eq 1 then begin
          x_left = x_left + n_x_data
        endif else begin
          x_left = 0
        endelse
      endif else if x_id_nan[i_nan] + square ge n_x_data then begin
        if wrap_x_opt eq 1 then begin
          x_right = x_right - n_x_data
        endif else begin
          x_right = n_x_data - 1
        endelse
      endif
      if x_left lt x_right then begin
        square_val = [ temp_data[x_left:x_right,y_bottom], $
            temp_data[x_left:x_right,y_top] ]
      endif else begin
        square_val = [ temp_data[x_left:n_x_data-1,y_bottom], $
            temp_data[0:x_right,y_bottom], temp_data[x_left:n_x_data-1,y_top], $
            temp_data[0:x_right,y_top] ]
      endelse
      if abs( y_top - y_bottom ) ge 2 then begin
        if y_bottom lt y_top then begin
          square_val = [ square_val, $
              reform( temp_data[x_left,y_bottom+1:y_top-1] ), $
              reform( temp_data[x_right,y_bottom+1:y_top-1] ) ]
        endif else begin
          square_val = [ square_val, $
              reform( temp_data[x_left,y_bottom+1:n_y_data-1] ), $
              reform( temp_data[x_left,0:y_top-1] ), $
              reform( temp_data[x_right,y_bottom+1:n_y_data-1] ), $
              reform( temp_data[x_right,0:y_top-1] ) ]
        endelse
      endif
      ; Determine if any of these points have a value
      id_good = where( finite( square_val ) eq 1, n_good )
      ; If there are non-nan values
      if n_good ge 1 then begin
        ; Adopt the average value
        data[x_id_nan[i_nan],y_id_nan[i_nan],i_z_data] $
            = mean( square_val, nan=1 )
        ; Take us out of the loop
        flag_check = 1
      ; Otherwise move on to the next square
      endif else begin
        square = square + 1
      endelse
    endwhile
  endfor
endfor

;***********************************************************************
; The end

return
END
