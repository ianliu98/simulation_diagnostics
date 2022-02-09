;+
; NAME:
;    QUANTILE_THRESHOLD
;
; PURPOSE:
;    This function estimates the thresholds in the input data array 
;    corresponding to the specified quantiles.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    result = quantile_threshold( data, quantile )
;
; INPUTS:
;    DATA:  A numerical array of data, of size N_DATA.
;    QUANTILE:  A floating point vector listing quantiles to find in DATA.  
;        Values must be within [0,1].  Of length N_QUANTILE
;
; KEYWORD PARAMETERS:
;    PRESORTED:  If set then the function assumes that the values in DATA have 
;        already been sorted in ascending order, thus running more efficiently. 
;        The default is to assume that they have not been sorted and thus to 
;        spend some time sorting them.
;
; OUTPUTS:
;    RESULT:  A numerical array containing the N_QUANTILE threshold values 
;        corresponding to the N_QUANTILE quantiles in DATA as specified in 
;        QUANTILE.
;
; USES:
;    -
;
; PROCEDURE:
;    This function estimates the thresholds corresponding to the quantiles 
;    using method #5 as documented in R.
;
; EXAMPLE:
;    data = randomn( 2, 100 )
;    quantile = [ 0.05, 0.95 ]
;    result = quantile_threshold( data, quantile )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dstone@lbl.gov), 2012-05-31
;    Modified:  DAS, 2012-07-05 (Added PRESORTED keyword)
;-

FUNCTION QUANTILE_THRESHOLD, $
    DATA, QUANTILE, $
    PRESORTED=presorted_opt

;***********************************************************************
; Constants and checks

; Ensure valid quantile values
if ( min( quantile, max=temp ) lt 0. ) or ( temp gt 1. ) then stop

; Count the quantiles
n_quantile = n_elements( quantile )

;***********************************************************************
; Estimate the thresholds

; Initialise the output vector of thresholds
threshold = fltarr( n_quantile )

; Identify valid values in DATA
id_data = where( finite( data ) eq 1, n_id_data )
if n_id_data eq 0 then stop
; Sort good values in data
data_sort = data[id_data]
if not( keyword_set( presorted_opt ) ) then begin
  id = sort( data_sort )
  data_sort = data_sort[id]
endif

; Define the index of the lower values closest to the thresholds
index = quantile * n_id_data - 0.5
; Deal with unresolved tail quantiles
id = where( index lt 0., n_id )
if n_id gt 0 then threshold[id] = data_sort[0]
id = where( index gt n_id_data-1., n_id )
if n_id gt 0 then threshold[id] = data_sort[n_id_data-1]
id = where( ( index ge 0. ) and ( index le n_id_data - 1. ), n_id )
index = index[id]
index_floor = floor( index )
threshold[id] = data_sort[index_floor] $
    + ( index - index_floor ) $
    * ( data_sort[index_floor+1] - data_sort[index_floor] )

;***********************************************************************
; The end

return, threshold
END
