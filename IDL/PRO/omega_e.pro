;+
;
; PURPOSE:
;   Calculation of an Electron gyrofrequency
;
;-

FUNCTION   omega_e,   b0

     ; Load constant quantities
     GLBVAR, var

     omega_e = ABS( var.qm * b0 )   ; q * B0 / m

     return,   omega_e

END
