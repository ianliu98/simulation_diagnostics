;+
; NAME:
;    algebra_to_matrix
;
; COPYRIGHT:
;    Copyright (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This function parses algebraic statements into a matrix for use in matrix 
;    algebra.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.1 and later are compatible with GDL.
;
; CALLING SEQUENCE:
;    Result = algebra_to_matrix( Instruct )
;
; INPUTS:
;    INSTRUCT:  A required vector of length N_INSTRUCT containing algebraic 
;        instructions.  Two formats are possible.  If of type string then the 
;        format is of the form '{sign}{index}' or '{factor}*{index}'.  {sign} 
;        can be '+', '-', or '';  {factor} can be any integer or floating point 
;        number, including negative values;  and {index} is the index number of 
;        the variable, according to where it is to appear in the matrix.  
;        Sequences of these command formats are allowed too.  So 
;        ['0.5*0','-1+0.5*0'] would return coefficients for half of the first 
;        variable and the difference of half the first variable and all the 
;        second variable.  If type integer, then the input must be positive 
;        with the digits corresponding to the indices of the variables to add.  
;        Thus, 13 in the integer format is equivalent to '0+2' or '+1*0+1*2' in 
;        the string format.
;
; KEYWORD PARAMETERS:
;    N_VAR:  An optional integer input giving the number of algebraic 
;        variables.  If not given then N_VAR is estimated from highest variable 
;        index number provided in INSTRUCT.
;
; OUTPUTS:
;    Result:  A matrix of size N_VAR*N_INSTRUCT where N_VAR is the number of 
;        algebraic variables.  Element [i_var,i_instruct] gives the factor by 
;        which variable i_var should be multiplied in the algebraic instruction 
;        i_instruct.
;
; USES:
;    ---
;
; PROCEDURE:
;    This function parses through string or integer instructions and interprets 
;    those instructions as multiplicative factors and indices for insertion in 
;    a algebraic matrix.
;
; EXAMPLE:
;    Result = algebra_to_matrix( ['0.5*0','-2+0.5*0'] )
;      Result should be:
;      0.50000      0.00000      0.00000
;      0.50000      0.00000     -1.00000
;    Result = algebra_to_matrix( [13,2] )
;      Result should be:
;      1.00000      0.00000      1.00000
;      0.00000      1.00000      0.00000
;    Also see gendetec.pro for a demonstration of use.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@csag.uct.ac.za), 2011-10-14, based on 
;        code in gendetec.pro by Myles Allen (v.3.1.0).
;    Modified:  DAS, 2012-02-15 (Added N_VAR keyword)
;-

;***********************************************************************

FUNCTION ALGEBRA_TO_MATRIX, $
    INSTRUCT, $
    N_VAR=n_var_0

;***********************************************************************
; Constants

; Number of rotation requests
n_instruct = n_elements( instruct )

; Copy the stated number of algebraic variables
if keyword_set( n_var_0 ) then n_var = n_var_0

; Determine the input format type
instruct_format = size( instruct )
temp = n_elements( instruct_format )
instruct_format = instruct_format[temp-2]

;***********************************************************************
; For string-format input

; If the instruct input is in more capable string format
if instruct_format eq 7 then begin
  ; Iterate through instructions
  for i_instruct = 0, n_instruct - 1 do begin
    ; Prepare the instruction for parsing
    temp = strmid( instruct[i_instruct], 0, 1 )
    if ( temp ne '+' ) and ( temp ne '-' ) then begin
      instruct[i_instruct] = '+' + instruct[i_instruct]
    endif
    ; Parse the instruction
    pos_instruct = strsplit( instruct[i_instruct] + '+1', '+-' )
    ; Iterate through parts of the instruction
    for i_pos = 0, n_elements( pos_instruct ) - 2 do begin
      ; Extract this part of the instruction
      temp_instruct = strmid( instruct[i_instruct], pos_instruct[i_pos] - 1, $
          pos_instruct[i_pos+1] - pos_instruct[i_pos] )
      ; Extract a multiplicative scaling if given
      if strpos( temp_instruct, '*' ) ge 0 then begin
        temp_instruct = strsplit( temp_instruct, '*', extract=1 )
      ; Otherwise take a factor of -1 or +1
      endif else begin
        temp_instruct = [ strmid( temp_instruct, 0, 1 ) + '1', $
            strmid( temp_instruct, 1, $
            pos_instruct[i_pos+1] - pos_instruct[i_pos] - 1 ) ]
      endelse
      ; Identify the variable
      id_var = fix( temp_instruct[1] )
      ; Initialise or expand the output matrix if necessary
      if n_elements( result ) eq 0 then begin
        if not( keyword_set( n_var_0 ) ) then n_var = id_var + 1
        result = fltarr( n_var, n_instruct )
      endif else if id_var ge n_var then begin
        temp = result
        result = fltarr( id_var + 1, n_instruct )
        result[0:n_var-1,*] = temp
        temp = 0
        if not( keyword_set( n_var_0 ) ) then n_var = id_var + 1
      endif
      ; Record this part of the transformation
      result[id_var,i_instruct] = float( temp_instruct[0] )
    endfor
  endfor
endif
  
;***********************************************************************
; For integer-format input (for back-compatability with gendetec.pro)

; If the instruct input is in the old integer format
if instruct_format eq 2 then begin
  ; Estimate the number of variables
  if not( keyword_set( n_var_0 ) ) then begin
    n_var = max( fix( alog10( instruct ) ) ) + 1
  endif
  ; Initialise output matrix
  result = fltarr( n_var, n_instruct )
  ; Iterate through instructions
  for i_instruct = 0, n_instruct - 1 do begin
    ; Unpack the integers in instruct[i_instruct]
    for i_pos = 0, fix( alog10( instruct[i_instruct] ) ) do begin
      id_var = fix( instruct[i_instruct] / 10 ^ i_pos ) $
          - 10 * fix( instruct[i_instruct] / 10 ^ ( i_pos + 1 ) ) - 1
      if id_var ge n_var then begin
        temp = result
        result = fltarr( id_var + 1, n_instruct )
        result[0:n_var-1,*] = temp
        temp = 0
        if not( keyword_set( n_var_0 ) ) then n_var = id_var + 1        
      endif
      result[id_var,i_instruct] = 1.
    endfor
  endfor
endif

;***********************************************************************
; The end

return, result
END
