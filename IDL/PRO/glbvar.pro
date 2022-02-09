;+
; PURPOSE:
;    Return constant and global variables
;
; OUTPUT
;    Structure
;
;    M.Hikishima  Jun 12, 2012
;-


PRO glbvar, var


var = {  $
         head:  Bytarr(4)                  , $          ; definition of header to read the data
         wsz:   [1024, 768]                , $          ; display size to plot
         dir:   '../dat/'                  , $          ; data directory
         omgc:  '!9' +STRING("167B)+ '!X'  , $          ; omega character (see Font symbol)
         lomgc: '!9' +STRING("127B)+ '!X'  , $          ; capital omega character
         pec:   '!9' +STRING("136B)+ '!X'  , $          ; perpendicular
         sigma: '!9' +STRING("163B)+ '!X'  , $          ; sigma
         tau:   '!9' +STRING("164B)+ '!X'  , $          ; tau
         me:    9.10938291d-31             , $          ; electron rest mass
         cvv:   2.99792458d8               , $          ; speed of light
         qe:    1.60217656d-19             , $          ; charge of electron
         qm:    -1.d0                         $          ; ratio of charge to mass in KEMPO
      }


END
