;+
; NAME:
;    TREND
;
; PURPOSE:
;    This function calculates the least-squares linear trend in the input 
;    vector.
;
; CATEGORY:
;    Time Series Analysis
;
; CALLING SEQUENCE:
;    result = trend( [DATA_X,] DATA_Y )
;
; INPUTS:
;    DATA_Y:  The input vector of type integer or floating point and length 
;        N_DATA.
;
; OPTIONAL INPUTS:
;    DATA_X:  An optional ordinate vector of type integer or floating point.  
;        It must be of length N_DATA.  If not set, DATA_Y's elements are 
;        positioned at [0,1,2,...,N-1].
;
; KEYWORD PARAMETERS:
;    ALPHA:  The significance level for statistics output.  The default value 
;        is 0.05.
;    CONFINT:  Returns the half-width of the (1-ALPHA) confidence interval of 
;        RESULT as a scalar floating-point number.  The confindence interval is 
;        calculated as [RESULT-CONFINT,RESULT+CONFINT].
;    CTREND:  The critical trend for the calculation of the one-sided 
;        significance level of the output trend.  The default is 0.  The 
;        significance level is returned in SIGLEV.
;    FIT:  Returns a vector of floating-point numbers containing the linear 
;        least-squares fit values, with DATA_X as the ordinate.
;    NAN:  If set, the function ignores NaNs (Not-a-Number) as missing values 
;        if any appear in DATA_X or DATA_Y.  The default is to return NaN.
;    SIGLEV:  Returns the significance level of the difference between RESULT 
;        and CTREND as a scalar floating-point number.  If TWOSIDED is set, 
;        this returns the two-sided significance level, otherwise it returns 
;        the one-sided significance level.
;    TWOSIDED:  If set, all significance tests are two-sided.  The default is 
;        for one-sided tests.  Note this does not affect CONFINT which is 
;        inherently two-sided.
;
; OUTPUTS:
;    Result:  Returns the least-squares linear trend of DATA_Y(DATA_X) as a 
;        floating-point number.
;    CONFINT, FIT, SIGLEV
;
; USES:
;    -
;
; PROCEDURE:
;    This function uses IDL's POLY_FIT.pro to calculate the least-squares 
;    linear trend to Y(X).  It also calculates statistics of the calculated 
;    trend value.
;
; EXAMPLE:
;    Create a vector of 20 random Gaussian-distributed data elements.
;      data_y = randomn( seed, 20 )
;    Calculate the linear trend in data_y.
;      result = trend( data_y )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (stoned@atm.ox.ac.uk), 2000-07-12.
;    Modified:  DAS, 2001-04-03 (fixed bug with SIGLEV and CONFINT).
;    Modified:  DAS, 2005-08-05 (replaced SUM.PRO use with TOTAL;  removed 
;        CONSTANTS.PRO use; edited coding style;  reversed definition of NAN 
;        keyword for consistency with other routines)
;    Modified:  DAS, 2005-11-23 (fixed bug with NAN option).
;    Modified:  DAS, 2015-05-28 (added TWOSIDED keyword;  modified formating;  
;        updated documentation)
;    Modified:  DAS, 2015-06-03 (fixed bug in recent formating modification)
;    Modified:  DAS, 2015-09-15 (fixed bug from June reformating which meant 
;        that FIT was not being returned)
;-

;***********************************************************************

FUNCTION TREND, $
    DATA_X, DATA_Y, $
    ALPHA=alpha, CTREND=ctrend, $
    CONFINT=confint, FIT=fit_y, SIGLEV=siglev, $
    TWOSIDED=twosided_opt, $
    NAN=nan_opt

;***********************************************************************
; Variables, Constants, and Options

; Load constants
nan = !values.f_nan

; Input vector
if n_params() eq 1 then begin
  data_y_use = data_x
endif else begin
  data_y_use = data_y
endelse
n_data = n_elements( data_y_use )

; Ordinate vector
if n_params() eq 1 then begin
  data_x_use = findgen( n_data )
endif else begin
  data_x_use = data_x
endelse
if n_elements( data_x_use ) ne n_data then begin
  print, 'Error in trend.pro:  Inconsistent input vector lengths'
  stop
endif

; The significance level
if not( keyword_set( alpha ) ) then alpha = 0.05

; The critical trend value for significance level calculation
if not( keyword_set( ctrend ) ) then ctrend = 0.

; Not-a-Number option
if keyword_set( nan_opt ) then begin
  nan_opt = 1
endif else begin
  nan_opt = 0
endelse

; Statistics output
if arg_present( confint ) then begin
  confint_opt = 1
endif else begin
  confint_opt = 0
endelse
if arg_present( siglev ) then begin
  siglev_opt = 1
endif else begin
  siglev_opt = 0
endelse

;***********************************************************************
; Setup

; Remove NaNs if NAN_OPT is set
if nan_opt eq 1 then begin
  ; Search for good values
  id = where( ( finite( data_y_use ) eq 1 ) and ( finite( data_x_use ) eq 1 ), $
      n_id )
  ; If we have some good values then remove the bad ones
  if n_id gt 1 then begin
    data_x_use = data_x_use[id]
    data_y_use = data_y_use[id]
  ; If we have no good values then return
  endif else begin
    return, nan
  endelse
endif

;***********************************************************************
; Calculate Trend

; Subtract mean (gives better accuracy)
mean_y = mean( data_y_use )
data_y_use = data_y_use - mean_y
mean_x = mean( data_x_use )
data_x_use = data_x_use - mean_x

; Calculate trend
trend_val = ( poly_fit( data_x_use, data_y_use, 1, fit_y ) )[1]
fit_y = fit_y + mean_y

;***********************************************************************
; Statistics

; Calculate confidence interval
if ( confint_opt eq 1 ) or ( siglev_opt eq 1 ) then begin
  se = sqrt( ( total( ( data_y_use + mean_y - fit_y ) ^ 2 ) ) $
      / ( n_data - 2. ) )
  sxx = sqrt( total( data_x_use ^ 2 ) )
  confint = t_cvf( alpha / 2., n_data - 2 ) * se / sxx
endif

; Calculate significance level
if siglev_opt eq 1 then begin
  t_val = ( trend_val - ctrend ) * sxx / se
  siglev = 1 - t_pdf( t_val, n_data - 2 )
  ; Convert to a two-sided significance level if requested
  if keyword_set( twosided_opt ) then begin
    if siglev gt 0.5 then siglev = 1. - siglev
    siglev = siglev * 2.
  endif
endif

;***********************************************************************
; The End

return, trend_val
END
