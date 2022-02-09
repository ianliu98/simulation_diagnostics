;+
;
; PURPOSE:
;       This is a calculation of \xi^2 which is dimensionless parameter related to despersion relation of plasa wave ([Omura et al., 2007])
;
; INPUT:
;       omega:   wave frequency
;       b0: ambient magnetic field B0
;
;-

FUNCTION   Xi,   omega, $
                 omega_pe, $
                 b0

     ; Load constant quantities
     GLBVAR, var


     ; Local electron gyrofrequency
     omega_e = omega_e( b0 )

     ; xi^2
     xi2 = omega * (omega_e - omega) / omega_pe^2
      
     return,   SQRT( xi2 )

END
