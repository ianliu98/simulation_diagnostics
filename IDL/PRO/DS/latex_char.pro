;+
; NAME:
;    latex_char.pro
;
; PURPOSE:
;    This function substitutes LaTeX character codes for non-English characters.
;
; CATEGORY:
;    Strings
;
; CALLING SEQUENCE:
;    result = latex_char( str_in )
;
; INPUTS:
;    STR_IN:  A required string or string array within which to replace all 
;        supported non-English characters with their LaTeX codes.
;
; KEYWORD PARAMETERS:
;    ---
;
; OUTPUTS:
;    RESULT
;
; USES:
;    string_substitute.pro
;
; PROCEDURE:
;    This function uses string_substitute.pro to perform the substition of a 
;    list of various characters to their respective LaTeX codes..
;
; EXAMPLE:
;    result = latex_char( 'Dáithí' )
;    ; This should return "D\'aith\'{\i}" (without the quotes).
;
; MODIFICATION HISTORY:
;    Written:  Daithi A. Stone (dastone@runbox.com), 2018-06-03
;-

;***********************************************************************

FUNCTION LATEX_CHAR, $
    STR_IN

;***********************************************************************
; Constants

; Copy input string
str_out = str_in

; Generate a conversion chart
convert_chart = [ $
    [ 'Á', "\'A" ], $
    [ 'á', "\'a" ], $
    [ 'À', '\`A' ], $
    [ 'à', '\`a' ], $
    [ 'Â', '\^A' ], $
    [ 'â', '\^a' ], $
    [ 'Ã', '\~A' ], $
    [ 'ã', '\~a' ], $
    [ 'Ä', '\"A' ], $
    [ 'ä', '\"a' ], $
    [ 'Ā', '\=A' ], $
    [ 'ā', '\=a' ], $
    [ 'å', '\r{a}' ], $
    [ 'Ć', "\'C" ], $
    [ 'ć', "\'c" ], $
    [ 'Č', '\v{C}' ], $
    [ 'č', '\v{c}' ], $
    [ 'Ç', '\c{C}' ], $
    [ 'ç', '\c{c}' ], $
    [ 'É', "\'E" ], $
    [ 'é', "\'e" ], $
    [ 'È', "\`E" ], $
    [ 'è', "\`e" ], $
    [ 'Ê', '\^E' ], $
    [ 'ê', '\^e' ], $
    [ 'Ë', '\"E' ], $
    [ 'ë', '\"e' ], $
    [ 'Ē', '\=E' ], $
    [ 'ē', '\=e' ], $
    [ 'Ḥ', '\d{H}' ], $
    [ 'Í', "\'I" ], $
    [ 'í', "\'{\i}" ], $
    [ 'Ì', '\`I' ], $
    [ 'ì', '\`{\i}' ], $
    [ 'Î', '\^I' ], $
    [ 'î', '\^{\i}' ], $
    [ 'Ï', '\"I' ], $
    [ 'ï', '\"{\i}' ], $
    [ 'Ł', '\L ' ], $
    [ 'ł', '\l ' ], $
    [ 'Ñ', '\~N' ], $
    [ 'ñ', '\~n' ], $
    [ 'Ó', "\'O" ], $
    [ 'ó', "\'o" ], $
    [ 'Ò', '\`O' ], $
    [ 'ò', '\`o' ], $
    [ 'Ô', '\^O' ], $
    [ 'ô', '\^o' ], $
    [ 'Ö', '\"O' ], $
    [ 'ö', '\"o' ], $
    [ 'Ō', '\=O' ], $
    [ 'ō', '\=o' ], $
    [ 'Õ', '\~O' ], $
    [ 'õ', '\~o' ], $
    [ 'Ø', '\O ' ], $
    [ 'ø', '\o ' ], $
    [ 'Ú', "\'U" ], $
    [ 'ú', "\'u" ], $
    [ 'Ù', '\`U' ], $
    [ 'ù', '\`u' ], $
    [ 'Û', '\^U' ], $
    [ 'û', '\^u' ], $
    [ 'Ü', '\"U' ], $
    [ 'ü', '\"u' ], $
    [ 'ū', '\=u' ], $
    [ 'Š', '\v{S}' ], $
    [ 'š', '\v{s}' ] ]
n_convert = n_elements( convert_chart[0,*] )

;***********************************************************************
; Perform conversion

; Iterate through supported characters
for i_convert = 0, n_convert - 1 do begin
  ; Substitute the LaTeX code for this character
  str_out = string_substitute( str_out, convert_chart[0,i_convert], $
      convert_chart[1,i_convert], regex=1, robust=1 )
endfor

;***********************************************************************
; End

return, str_out
END
