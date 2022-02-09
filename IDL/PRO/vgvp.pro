;+
; PURPOSE:
;       This gives a Lorentz factor
;
;                                                    2016/11/08 M.Hikishima @ ISAS
;-

c = 100.
omg_e0 = 1.
wp = 15.
omg = findgen(100)/100

xi2 = omg * (omg_e0 - omg) / wp^2
delta2 = 1./ (1. + xi2)
delta = sqrt( delta2 )
xi = sqrt( xi2 )

vp = c * delta * xi
vg = c * xi / delta * (xi2 + omg_e0/2./(omg_e0-omg) )^(-1)

vp = vp/ c
vg = vg/ c

maxv = max( [max(vp), max(vg)] )

cgplot, omg, vp, xtitle='omega', ytitle='Vg,Vp/c', yrange=[0, maxv]
cgplot, omg, vg, linestyle=1, /overplot

stop
end
