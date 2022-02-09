; calculate S

FUNCTION scal, kmp, fre

Q   = 0.05
wpe = kmp.wp1
wph = kmp.wp2
cv  = kmp.cv / 100.0
utpara = kmp.path2 / 100.0
utperp = kmp.peth2 / 100.0
beta   = 0.3
uperph  = sqrt(!PI/2.0) * ((1.0 - beta^(1.5) / (1.0 - beta))) * utperp
vperp  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * uperph
vpara  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * utpara
gamma  = 1.0 / sqrt(1.0 - (vperp^2.0 + vpara^2.0)/cv^2.0);
Omega_e = 1.00

xi2 = abs(fre) * (Omega_e - abs(fre)) / wpe^2.0
xi2 = abs(xi2)
delta2 = 1.0 / (1.0 + xi2)

w_s = fre / Omega_e
vp_s = sqrt(xi2) * sqrt(delta2)
vperp_s = vperp / cv
vr_s = (w_s^2.0 - sqrt(w_s^4.0 + (w_s^2.0 + vp_s^2.0) * (1.0 - w_s^2.0 - vperp_s^2.0))) / (w_s^2.0 / vp_s + vp_s)
vp = vp_s * cv
vr = vr_s * cv
vr2 = (1.0 - Omega_e / (gamma * fre)) * vp

vg = (cv * sqrt(xi2) / sqrt(delta2)) / (xi2 + Omega_e / (2.0 * (Omega_e - fre)))

s0 = sqrt(delta2) * vperp / (sqrt(xi2) * cv)
s1 = gamma * (1.0 - vr / vg)^2.0

s2_term1 = 1.0 / (2.0 * sqrt(xi2) * sqrt(delta2))
s2_term2 = gamma * fre * (vperp/cv)^2.0 / Omega_e
s2_term3 = (2.0 + (fre / Omega_e) * delta2 * (Omega_e - gamma * fre) / (Omega_e - fre)) * vr * vp / cv^2.0
s2 = s2_term1 * (s2_term2 - s2_term3)

RETURN, s0

END

