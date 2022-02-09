;+
; NAME:
;    SIGN
;
; PURPOSE:
;    This function returns the sign of the input variable.
;
; CATEGORY:
;    Mathematics, Optimal Detection Package v3.1.1
;
; CALLING SEQUENCE:
;    Result = SIGN( DATA )
;
; INPUTS:
;    DATA:  A scalar or array of type integer or floating point.
;
; OUTPUTS:
;    Result:  Returns the sign of DATA.
;
; PROCEDURE:
;    This function determines whether DATA <, or > 0.
;
; EXAMPLE:
;    Determine the sign of -3.
;      result = sign( -3 )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone, 2000-06-09.
;    Modified:  DAS, 2000-07-06 (removed use of length.pro).
;    Modified:  DAS, 2000-07-10 (removed for loops).
;    Modified:  DAS, 2000-07-24 (decided 0 is positive).
;    Modified:  DAS, 2011-11-06 (Inclusion in Optimal Detection Package 
;        category;  modified format and variable names)
;-

;***********************************************************************

FUNCTION SIGN, DATA

;***********************************************************************
; Define variables

; Size of input
n_data = n_elements( data )

; Initialise output
result = 0 * fix( data )

;***********************************************************************
; Determine the sign

; Find all positives
id = where( data ge 0, n_id )
if n_id gt 0 then result[id] = 1

; Find all negatives
id = where( data lt 0, n_id )
if n_id gt 0 then result[id] = -1

;***********************************************************************
;The End

return, result
END
