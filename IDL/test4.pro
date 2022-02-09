cv  = kmp.cv / 100.0
utpara = kmp.path2 / 100.0
utperp = kmp.peth2 / 100.0
betaa   = 0.3
uperph = sqrt(!PI/2.0) * ((1.0 - betaa^(1.5) / (1.0 - betaa))) * utperp
vperp  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * uperph
vpara  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * utpara
gammaa  = 1.0 / sqrt(1.0 - (vperp^2.0 + vpara^2.0)/cv^2.0)

we0 = 1.d
wpc = kmp.wp1
wph = kmp.wp2
wpe = sqrt(wpc^2 + wph^2)
cvv = kmp.cv
dww = 3e-6
www = findgen(100000L) * dww + dww
xi2 = www * (we0 - www) / wpe^2
xi2 = abs(xi2)
xii = sqrt(xi2)
ch2 = 1.d / (1.d + xi2)
chi = sqrt(ch2)
kww = www / (cvv * chi * xii)
kww = kww * cvv
kma = kww[-1]

kk = n_elements(where(ku le kma))
ww = dblarr(kk)
xi2_s = dblarr(kk)
xii_s = dblarr(kk)
ch2_s = dblarr(kk)
chi_s = dblarr(kk)
for ik=1,kk do begin
	ww[ik-1] = www[n_elements(where(kww le ku[ik]))-1]
	xi2_s[ik-1] = xi2[n_elements(where(kww le ku[ik]))-1]
	xii_s[ik-1] = xii[n_elements(where(kww le ku[ik]))-1]
	ch2_s[ik-1] = ch2[n_elements(where(kww le ku[ik]))-1]
	chi_s[ik-1] = chi[n_elements(where(kww le ku[ik]))-1]
endfor

xi2_t = dblarr(kk, ret)
xii_t = dblarr(kk, ret)
ch2_t = dblarr(kk, ret)
chi_t = dblarr(kk, ret)
www_t = dblarr(kk, ret)
kkk_t = dblarr(kk, ret)

for it=0,ret-1 do begin
	xi2_t[*,it] = xi2_s[*]
	xii_t[*,it] = xii_s[*]
	ch2_t[*,it] = ch2_s[*]
	chi_t[*,it] = chi_s[*]
	www_t[*,it] = ww[*]
	kkk_t[*,it] = ku[0:kk-1]
endfor

fff = field[0:kk-1, *]

wtr = chi_t * sqrt( abs(kkk_t * vperp * fff / gammaa) )
delta_w1 = 4.d * wtr
delta_w2 = ch2_t * (we0 / (gammaa * www_t) - 1.d)
delta_w3 = xi2_t + we0 / (2.d * (we0 - www_t))
delta_w = delta_w1 / (1.d + delta_w2 * delta_w3)

delta_k = delta_w / (cv * xii_t * chi_t)




end
