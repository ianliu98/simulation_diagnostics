;+
; NAME:
;    MODE
;
; PURPOSE:
;    This function calculates the mode of an input array.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    RESULT = mode( DATA )
;
; INPUTS:
;    DATA:  A required input array of any numerical type, of any dimensino.
;
; KEYWORD PARAMETERS:
;    EXCLUDE_NAN:  If set, then the function will exclude NaN values from the 
;        analysis.  This cannot be set at the same time as INCLUDE_NAN.  The 
;        more efficient default is to be ignorant of NaNs, but with the 
;        possibility of strange behaviour when NaNs are actually present.
;    INCLUDE_NAN:  If set, then the function will consider NaN as a possible 
;        value to output.  This cannot be set at the same time as EXCLUDE_NAN. 
;        The more efficient default is to be ignorant of NaNs, but with the 
;        possibility of strange behaviour when NaNs are actually present.
;
; USES:
;
; PROCEDURE:
;    This function sort the values in the array, identifies where neighbours 
;    differ, and then the largest distance between pairs of differing 
;    neighbours.
;
; EXAMPLE:
;    ; The following should return 4.
;    result = mode( [ 2, 4, 5, 4, 3 ] )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (dstone@lbl.gov), 2016-10-20.
;-

;***********************************************************************

FUNCTION MODE, $
    DATA, $
    EXCLUDE_NAN=exclude_nan_opt, INCLUDE_NAN=include_nan_opt

;***********************************************************************
; Constants and options

; Determine input array size
n_data = n_elements( data )
if n_data eq 0 then begin
  print, 'ERROR mode.pro:  No data input.'
  stop
endif

; Option to consider NaN as a possible result
include_nan_opt = keyword_set( include_nan_opt )
; Option to exclude NaNs from analayis
exclude_nan_opt = keyword_set( exclude_nan_opt )
; Ensure both options are not set
if exclude_nan_opt + include_nan_opt eq 2 then begin
  print, 'ERROR mode.pro:  Only one of EXCLUDE_NAN and INCLUDE_NAN options can be set at once.'
  stop
endif

;***********************************************************************
; Calculate the mode

; Copy the array
temp_data = data

; If we are including NaNs in the count
if include_nan_opt eq 1 then begin
  ; Convert NaNs to a numerical value
  id = where( finite( temp_data ) eq 0, n_id )
  if n_id gt 0 then begin
    marker_nan = min( temp_data ) - 1
    temp_data[id] = marker_nan
  endif
endif

; If we are excluding NaNs
if exclude_nan_opt eq 1 then begin
  ; Remove NaNs
  id = where( finite( temp_data ) eq 1, n_id )
  if n_id eq 0 then begin
    print, 'ERROR mode.pro:  No non-NaN values in array.'
    stop
  endif else if n_id ne n_data then begin
    temp_data = temp_data[id]
  endif
endif

; Sort the values
id_sort = sort( temp_data )
temp_data = temp_data[id_sort]
; Identify non-identical neighbours
id = where( temp_data ne shift( temp_data, -1 ), n_id )
; If all values are identical
if n_id eq 0 then begin
  ; Take that value
  result = temp_data[0]
; If there are different values
endif else begin
  ; Identify the most frequent value
  ; (In case of ties, the smaller value is taken)
  temp = max( id - [ -1, id ], id_max )
  result = temp_data[id[id_max]]
endelse

; Set result to NaN if that was the most frequent value, if EXCLUDE_NAN is not 
; set
if n_elements( marker_nan ) eq 1 then begin
  if result eq marker_nan then result = !values.f_nan
endif

;***********************************************************************
; The End

return, result
END
