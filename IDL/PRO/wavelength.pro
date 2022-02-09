;+
;
; PURPOSE:
;   Calculation of a wave length
;      wavenumber * wave length = 2pi
;
;-

FUNCTION   wavelength,   omega, $
                         wp1, $
                         wp2, $
                         cv, $
                         b0  

  return, wavelength = 2.*!DPI / Wavenumber( omega, wp1, wp2, cv, b0 ) 

END
