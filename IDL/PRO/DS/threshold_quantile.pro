;+
; NAME:
;    THRESHOLD_QUANTILE
;
; PURPOSE:
;    This function estimates the quantiles of the given thresholds in the input 
;    data array.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    result = threshold_quantile( data, threshold )
;
; INPUTS:
;    DATA:  A numerical array of data, of size N_DATA.
;    THRESHOLD:  A numerical vector listing thresholds to find in DATA and for 
;        which to estimate their quantiles.  Of length N_THRESHOLD
;
; KEYWORD PARAMETERS:
;    PRESORTED:  If set then the function assumes that the values in DATA have 
;        already been sorted in ascending order, thus running more efficiently. 
;        The default is to assume that they have not been sorted and thus to 
;        spend some time sorting them.
;
; OUTPUTS:
;    RESULT:  A numerical array containing the N_THRESHOLD quantile values 
;        corresponding to the N_THRESHOLD thresholds in DATA as specified in 
;        THRESHOLD.
;
; USES:
;    -
;
; PROCEDURE:
;    This function essentially inverts the method used in 
;    quantile_threshold.pro.
;
; EXAMPLE:
;    data = randomn( 2, 100 )
;    threshold = [ -1., 1. ]
;    result = threshold_quantile( data, threshold )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dstone@lbl.gov), 2012-06-01
;    Modified:  DAS, 2012-07-05 (Added PRESORTED keyword)
;-

FUNCTION THRESHOLD_QUANTILE, $
    DATA, THRESHOLD, $
    PRESORTED=presorted_opt

;***********************************************************************
; Constants and checks

; Count the thresholds
n_threshold = n_elements( threshold )

;***********************************************************************
; Estimate the quantiles

; Initialise the output vector of thresholds
quantile = fltarr( n_threshold )

; Identify valid values in DATA
id_data = where( finite( data ) eq 1, n_id_data )
if n_id_data eq 0 then stop
; Sort good values in data
data_sort = data[id_data]
if not( keyword_set( presorted_opt ) ) then begin
  id = sort( data_sort )
  data_sort = data_sort[id]
endif

; Iterate through thresholds
for i_threshold = 0, n_threshold - 1 do begin
  ; Find the nearest data value lower than the threshold
  id_low = max( where( data_sort le threshold[i_threshold] ) )
  ; If there is no value lower
  if id_low eq -1 then begin
    ; Just take 0
    quantile[i_threshold] = 1. / n_id_data / 4.
  ; If there is no higher value
  endif else if id_low eq n_id_data - 1 then begin
    ; Just take 1
    quantile[i_threshold] = 1. - 1. / n_id_data / 4.
  ; If we are not at the ends of the samples
  endif else begin
    ; Interpolate quantile
    quantile[i_threshold] = ( float( id_low ) + 0.5 $
        + ( threshold[i_threshold] - data_sort[id_low] ) $
        / ( data_sort[id_low+1] - data_sort[id_low] ) ) $
        / n_id_data
  endelse
endfor

;***********************************************************************
; The end

return, quantile
END
