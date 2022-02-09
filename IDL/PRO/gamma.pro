;+
; PURPOSE:
;       This gives a Lorentz factor
;
;                                                    2016/11/08 M.Hikishima @ ISAS
;-
FUNCTION   Gamma,   vpa_n, $ 	; parallel velocity (normalized by the speed of ligtht)
                    vpe_n 		; perpendicular velocity (normalized by the speed of ligtht)


   ; Check parameters
   IF ((vpa_n GT 1.0D) OR (vpe_n GT 1.0D)) THEN   $
        PRINT, 'Velocities should be normalized by the speed of light !!!'

   gamma =   1.0D / SQRT( 1.0D - (vpa_n^2 + vpe_n^2) )


   RETURN,   gamma

END
