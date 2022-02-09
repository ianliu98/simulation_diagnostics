;+
;
; PURPOSE:
;     Calculation of resonance velocity
;
; INPUT:
;     omg: Frequency 
;     vpe; Perpendicular velocity (array or single)
;     wp_h: Plasma frequency (hot)
;     wp_c: Plasma frequency (cold)
;     b0: ambient magnetic field
;
; Equation:
;     Considering the relativistic factor
;     -   Vr = 1/k * (omg - omg_e/gamma)
;     -   gamma = 1/ sqrt( Vr^2 + Vpe^2 )
;     Slove the Vr from the equations Vr and gamma above
;
; RETURN:
;      
;     Resonance velocity (normalized by the speed of light)
;       The return values are poitive because Vp = omg/ +k
;       The return values are array if vpe is array
;
;                                              Jul 26, 2019 M.Hikishima @ISAS
;-

FUNCTION   ResV,   omg, $
                   vpe, $
                   wp_c, $
                   wp_h, $
                   b0

  ; load constant quantities
  GLBVAR, var

  wp = omega_pe( wp_c, wp_h )
  we = omega_e( b0 )
  xi = Xi( omg, wp, b0 )
  delta = Delta( xi )
  vp = Vp( delta, xi )

  den  = omg^2 + we^2 * vp^2
  num = omg^2 - Sqrt( omg^4 - den * (omg^2 - we^2 * (1.0 - vpe^2)) )
  vr = vp * num / den

  RETURN,   vr
END
