;+
; PURPOSE:
;     This gives the Phase velocity (normalized by the speed of light)
;
;                          2019/07/26 M.Hikishima @ISAS
;-

FUNCTION   Vp,     delta, $
                   xi

  return,   vp = delta * xi
end
