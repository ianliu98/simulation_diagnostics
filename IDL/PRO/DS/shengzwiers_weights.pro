;+
; NAME:
;    shengzwiers_weights.pro
;
; PURPOSE:
;    This function returns the interpolation weighting matrix used in the Sheng 
;    and Zwiers (1998) algorithm for variance adjustment for interpolation.
;
; CATEGORY:
;    Time Series Analysis
;
; CALLING SEQUENCE:
;    result = shengzwiers_weights( month_len )
;
; INPUTS:
;    MONTH_LEN:  An optional integer vector containing the number of days in  
;        each month during a period of interest.  If input then the number of 
;        values is taken as the size of the matrix.  If not given then the it 
;        is assumed that all months are the same length and the number of 
;        months is provided in the N_MONTH.  If not input then N_MONTH must be 
;        specified
;
; KEYWORD PARAMETERS:
;    DOUBLE:  If set then a double precision floating array is returned.  The 
;        default is single precision.
;    N_MONTH:  If MONTH_LEN is not given then this is a scalar integer giving 
;        the size for the matrix.  If MONTH_LEN is input then this is ignored.
;
; OUTPUTS:
;    Result:  An floating point array containing the 2-dimensional weighting 
;        matrix required for the Sheng and Zwiers algorithm.
;
; USES:
;    ---
;
; PROCEDURE:
;    This function sets the values in the weighting matrix as specified in 
;    Sheng and Zwiers (1998).
;
; EXAMPLE:
;    See shengzwiers.pro
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2016-09-21, as 
;        sheng_zwiers_weights.pro.
;    Modified:  DAS, 2017-09-20 (Branched to shengzwiers_weights.pro;  
;        standardised documentation)
;-

;***********************************************************************

FUNCTION SHENGZWIERS_WEIGHTS, $
    MONTH_LEN, $
    N_MONTH=n_month, $
    DOUBLE=double_opt

;***********************************************************************
; Create matrix

; Get the dimensions of the matrix
if keyword_set( month_len ) then n_month = n_elements( month_len )

; Ensure floating point lengths
if keyword_set( month_len ) then begin
  if keyword_set( double_opt ) then begin
    month_len = double( month_len )
  endif else begin
    month_len = float( month_len )
  endelse
endif

; Create an vector of the indices of the diagonal elements of the output matrix
index = indgen( n_month ) * ( n_month + 1 )

; Initialise the output matrix
if keyword_set( double_opt ) then begin
  result = dblarr( n_month, n_month )
endif else begin
  result = fltarr( n_month, n_month )
endelse

; If all months are the same length (faster)
if not( keyword_set( month_len ) ) then begin
  ; Do the end values
  result[0,0] = 7. / 8.
  result[n_month-1,n_month-1] = 7. / 8.
  ; Do the other values
  result[index[1:n_month-2]] = 3. / 4.
  result[index[0:n_month-2]+1] = 1. / 8.
  result[index[1:n_month-1]-1] = 1. / 8.
; If months have different lengths (calculate explicitly)
endif else begin
  ; Do the end values
  result[0,0] = ( 8. * month_len[1] + 6. * month_len[0] ) $
      / ( 8. * ( month_len[1] + month_len[0] ) )
  id = n_month - 1
  result[id,id] = ( 8. * month_len[id-1] + 6. * month_len[id] ) $
      / ( 8. * ( month_len[id-1] + month_len[id] ) )
  ; Do the other self-weightings
  if n_month gt 2 then begin
    id = 1 + indgen( n_month - 2 )
    result[index[id]] = ( 4. * month_len[id-1] * month_len[id+1] $
        + 3. * month_len[id] * ( month_len[id-1] + month_len[id+1] ) $
        + 2. * ( month_len[id] ^ 2. ) ) $
        / ( 4. * ( month_len[id-1] + month_len[id] ) $
        * ( month_len[id] + month_len[id+1] ) )
  endif
  ; Do the weighting for the prior months
  id = indgen( n_month - 1 )
  result[index[id]+1] = month_len[id+1] $
      / ( 4. * ( month_len[id] + month_len[id+1] ) )
  ; Do the weight for the following months
  result[index[id+1]-1] = month_len[id] $
      / ( 4. * ( month_len[id] + month_len[id+1] ) )
endelse

;***********************************************************************
; The end

return, result
END
