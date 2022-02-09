;+
;
; PURPOSE:
;   Calculation of a wavenumver
;     Vp = c * delta * xi = omega / wavenumber 
;
;                              Nov 8, 2016 M.Hikishima @ISAS
;-

FUNCTION   wavenumber,   omega, $
                         wp1, $
                         wp2, $
                         cv, $
                         b0
                   
   ; load variables
   GLBVAR, var

   wp = wp(wp1, wp2)
   xi = Xi( omega, wp, b0 )
   delta = Delta( xi ) 

   wavenumber = omega / (cv * delta * xi)

   return, wavenumber 

END
