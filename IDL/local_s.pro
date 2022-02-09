; plot S of a specific position

dname = ''
read, 'data file name: ', dname
filey = dname + '_by.csv'
filez = dname + '_bz.csv'

field_y = read_csv(filey)
field_y = field_y.FIELD1
field_z = read_csv(filez)
field_z = field_z.FIELD1

read, 'how many points needed: ', npnt
read, 'b02 = ', b02

;;; parameter ;;;
tskip = 4
xskip = 4
ifdiag = 16
dt = 0.004
ddt = dt * ifdiag * tskip
Q = 0.05
wpe = 15
wph = 0.3
cv = 1.0
utpara = 0.26
utperp = 0.3
beta = 0.3
uperph  = sqrt(!PI/2.0) * ((1.0 - beta^(1.5) / (1.0 - beta))) * utperp
vperp  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * uperph
vpara  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * utpara
gamma  = 1.0 / sqrt(1.0 - (vperp^2.0 + vpara^2.0)/cv^2.0);
Omega_e = 1.00
a = (b02/1.0 - 1.0) / (32768.0/100.0/2.0)^2.0
field_size = size(field_y, /DIMENSIONS)
S = Fltarr(field_size[1], npnt)
b_save = S
fre_save = S
dwdt_save = S

For i=0,npnt-1 do begin
	read, 'Enter the position of the' + string(i+1) +' point [0~28672] :', posi 
	fre = fre_cal(field_y[posi/xskip,*], field_z[posi/xskip,*], ddt)
	
	;; calculation S ;;
	xi2 = abs(fre * (Omega_e - fre) / wpe^2.0)
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
	
	dwdt = fre
	Len = N_ELEMENTS(fre)
	for j=2,Len-3 do begin
		dwdt[j] = (fre[j-2] - 8 * fre[j-1] + 8 * fre[j+1] - fre[j+2]) / (12 * ddt)
	endfor
	
	k1 = s1 * dwdt
	k2 = s2
	dbdh = 2.0 * a * abs((posi - 28672.0/2.0) / 100.0)
	k2 = k2 * dbdh * 1.0
	bw = sqrt(field_y[posi/xskip,*]^2 + field_z[posi/xskip,*]^2)
	k3 = s0 * fre * bw
	
	S[*,i] = -1.0 * (k1 + k2) / k3
	
	b_save[*,i] = bw
	fre_save[*,i] = fre
	dwdt_save[*,i] = dwdt

Endfor

END
