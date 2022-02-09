;+
;
; PURPOSE:
;       This is a calculation of \delta whici is dimensionless parameter related to dispersion relation
;       [Omura et al., 2007]
;
; INPUT:
;       \xi: Dimensionless parameter (call a function Xi)
;
;                                                        Nov 8, 2016 M.Hikishima @ISAS
;-

FUNCTION   Delta,   xi


   delta2 = 1.0D / (1.0D + xi^2)

   RETURN,   SQRT( delta2 )

END
