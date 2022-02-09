;+
;
; PURPOSE:
;       This is calculation of trapping velocity
;
; INPUT:
;       Below
; RETURN:
;
;                                                        Nov 8, 2016 M.Hikishima @ISAS
;-

FUNCTION   TrapV,   omega, $
                    vpe, $ 
                    bw_pe, $
                    omega_pe_cold, $
                    omega_pe_hot, $
                    b0


   ; load constant quantities
   GLBVAR, var

   ; Check parameters

   omega_e  = ABS( var.qm ) * b0	; local electron gyrofrequency
   xi		= Xi( omega, omega_pe_cold, omega_pe_hot, b0 )
   delta	= Delta( xi )
   

   vr = delta * xi * ( 1.0D - omega_e / (gamma * omega) )

   RETURN,   vr
END
