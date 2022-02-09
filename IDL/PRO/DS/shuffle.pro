;+
; NAME:
;    SHUFFLE
;
; PURPOSE:
;    This function shuffles the values in an array.
;
; CATEGORY:
;    Array, Optimal Detection Package v3.1.1
;
; CALLING SEQUENCE:
;    Result = SHUFFLE( X )
;
; INPUTS:
;    X:  An array of any type.
;
; KEYWORD PARAMETERS:
;    INDEX:  Returns the index values of the shuffled array.
;    REPLACE:  If set, the function shuffles with replacement.  The default is 
;        without replacement.
;    SEED:  A seed for the random number generator to use.
;
; OUTPUTS:
;    Result:  Returns the shuffled version of array X.
;    INDEX:  See above.
;    SEED:  See above.  Also returns a seed for the next implementation of the 
;        random number generator.
;
; PROCEDURE:
;    This function used the RANDOMU IDL function to produce an array of random 
;    index values.
;
; EXAMPLE:
;    Define a vector.
;      x = [0,1,2,3,4]
;    Shuffle the vector with replacement.
;      result = shuffle( x, replace=1, seed=1 )
;    This should give result = [2,0,3,2,4].
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@atm.ox.ac.uk), 2001-07-18.
;    Modified:  DAS, 2002-11-21 (added seed initialisation).
;    Modified:  DAS, 2003-02-17 (allowed input of very long vectors)
;    Modified:  DAS, 2005-01-03 (added SEED keyword)
;    Modified:  DAS, 2011-03-16 (corrected bug when long vectors are input;  
;        modified formating)
;    Modified:  DAS, 2011-11-06 (Inclusion in Optimal Detection Package 
;        category)
;-

;***********************************************************************

FUNCTION SHUFFLE, $
    X, $
    INDEX=index, $
    REPLACE=replace_opt, $
    SEED=seed

;***********************************************************************
; Constants and Variables

; Number of elements
n_x = n_elements( x )
; Whether we need to worry about long integer indices
long_opt = n_x gt 32000
if long_opt eq 1 then begin
  zero = 0l
endif else begin
  zero = 0
endelse

; Index of shuffled values
if long_opt eq 1 then begin
  index = lonarr( n_x )
endif else begin
  index = intarr( n_x )
endelse

; Working vector of index values
if long_opt eq 1 then begin
  id_left = lindgen( n_x + 2 ) - 1
endif else begin
  id_left = indgen( n_x + 2 ) - 1
endelse
n_left = n_elements( id_left )

; Replacement option
replace_opt = keyword_set( replace_opt )

; Initialise the random number generator
if not( keyword_set( seed ) ) then seed = long( systime( seconds=1 ) )

;***********************************************************************
; Shuffle Values

for i = zero, n_x - 1 do begin
  ; Determine how many values left to shuffle
  if not( replace_opt ) then n_left = n_elements( id_left )
  ; Get a random index value
  ind = floor( randomu( seed, uniform=1 ) * ( n_left - 2 ) ) + 1
  index[i] = id_left[ind]
  ; If not replacing values, remove the value
  if not( replace_opt ) then begin
    id_left = [ id_left[0:ind-1], id_left[ind+1:n_left-1] ]
  endif
endfor

; Create shuffled array
y = x[index]

; Reform output
x_dim = ( size( x ) )[1:n_elements(size(x))-3]
y = reform( y , x_dim )
index = reform( index , x_dim )

;***********************************************************************
; The End

return, y
END
