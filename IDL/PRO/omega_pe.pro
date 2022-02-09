;+
;
; PURPOSE:
;       This gives a electron plasma frequency
;
; INPUT:
;       omega:   wave frequency
;       b0: ambient magnetic field B0
;
;                                                        Nov 8, 2016 M.Hikishima @ISAS
;-

FUNCTION   Omega_pe,   omega_pe_cold, $
                       omega_pe_hot

     RETURN,   SQRT( omega_pe_cold^2 + omega_pe_hot^2 )

END
