;+
; NAME:
;    string_substitute
;
; PURPOSE:
;    This function substitutes a character sequence within an input string.
;
; CATEGORY:
;    Strings
;
; CALLING SEQUENCE:
;    result = string_substitute( str_in, str_remove, str_add )
;
; INPUTS:
;    STR_IN:  A required string or string array of strings within which to 
;        substitute character sequences.
;    STR_REMOVE:  A required string containing the character sequence to remove 
;        from STR_IN.
;    STR_ADD:  An optional string containing the character sequence with which 
;        to replace STR_REMOVE.
;
; KEYWORD PARAMETERS:
;    REGEX:  The strsplit REGEX keyword option when removing STR_REMOVE from 
;        STR_IN.  Performs the same effect whether ROBUST is set or not.
;    ROBUST:  If set then the function uses strpos.pro and strmid.pro instead 
;        of strsplit.pro in order to identify and remove cases of STR_REMOVE.  
;        This is less efficient, but strsplit.pro has issues if non-English 
;        characters are included in STR_REMOVE.
;
; OUTPUTS:
;    RESULT
;
; USES:
;    ---
;
; PROCEDURE:
;    This function uses strsplit.pro to remove the character sequence requested 
;    for removal, and strjoin.pro to add the replacement character sequence.
;
; EXAMPLE:
;    result = string_substitute( 'abcde', 'cd', 'fg' )
;    ; This should return 'abfge'.
;
; MODIFICATION HISTORY:
;    Written:  Daithi A. Stone (dastone@runbox.com), 2018-06-03
;-

;***********************************************************************

FUNCTION STRING_SUBSTITUTE, $
    STR_IN, STR_REMOVE, STR_ADD, $
    REGEX=regex_opt, $
    ROBUST=robust_opt

;***********************************************************************
; Constants

; The number of strings to work on
n_str_in = n_elements( str_in )
if n_str_in eq 0 then stop

; Copy input string to output
str_out = str_in

; Ensure a single character sequence has been included for removal
if n_elements( str_remove ) ne 1 then stop
if str_remove eq '' then stop
len_remove = strlen( str_remove )
; Ensure no more than a single character sequence has been included for addition
if n_elements( str_add ) gt 1 then stop
if n_elements( str_add ) eq 0 then str_add = ''

; Determine if we are to use the more robust method
robust_opt = keyword_set( robust_opt )
; The strsplit regex option (also effectively used with the strpos/strmid 
; method)
regex_opt = keyword_set( regex_opt )

; If STR_ADD equals STR_REMOVE then just return
if keyword_set( str_add ) eq 1 then begin
  if str_add eq str_remove then return, str_out
endif

;***********************************************************************
; Peform the substitution

; Iterate through input strings
for i_str_in = 0, n_str_in - 1 do begin
  ; Copy the string
  temp_str = str_out[i_str_in]
  ; The robust method using strpos and strmid
  if robust_opt eq 1 then begin
    ; Iterate until we have replaced all instances of STR_REMOVE
    pos = 0
    while pos ge 0 do begin
      ; Identify the next instance of STR_REMOVE
      pos = strpos( temp_str, str_remove )
      ; If an instance has been found
      if pos ge 0 then begin
        ; Determine the length of the string
        temp_len = strlen( temp_str )
        ; Identify the extraction points for the special case of the beginning 
        ; of STR_IN
        if pos eq 0 then begin
          if regex_opt eq 1 then begin
            pos_cut = [ 0, 0, len_remove, temp_len - len_remove ]
          endif else begin
            pos_cut = [ len_remove, temp_len - len_remove ]
          endelse
        ; Identify the extraction points for the special case of the end of 
        ; STR_IN
        endif else if pos eq temp_len - 1 then begin
          if regex_opt eq 1 then begin
            pos_cut = [ 0, temp_len - len_remove, temp_len, 0 ]
          endif else begin
            pos_cut = [ 0, temp_len - len_remove ]
          endelse
        ; Identify the extraction points for the general case
        endif else begin
          pos_cut = [ 0, pos, pos + len_remove, temp_len - pos - len_remove ]
        endelse
        ; If one segment is to be kept
        if n_elements( pos_cut ) eq 2 then begin
          ; Just remove the instance of STR_REMOVE
          temp_str = strmid( temp_str, pos_cut[0], pos_cut[1] )
        ; If two segments are to be kept
        endif else begin
          ; Replace the instance of STR_REMOVE with STR_ADD
          temp_str = strmid( temp_str, pos_cut[0], pos_cut[1] ) + str_add $
              + strmid( temp_str, pos_cut[2], pos_cut[3] )
        endelse
      endif
    endwhile
  ; The efficient method using strsplit
  endif else begin
    ; Extract any incidences of the string to be removed
    temp_str = strsplit( temp_str, str_remove, extract=1, count=n_temp_str, $
        regex=regex_opt )
    ; If there are incidences of that string
    if n_temp_str gt 0 then begin
      ; Replace the omitted with the replacement string
      if keyword_set( str_add ) then begin
        temp_str = strjoin( temp_str, str_add )
      endif else begin
        temp_str = strjoin( temp_str )
      endelse
    endif
  endelse
  ; Record modified string
  str_out[i_str_in] = temp_str
endfor

;***********************************************************************
; The End

return, str_out
END
