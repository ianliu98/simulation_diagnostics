LOADCT, 39

GLBVAR, var

wpia = 3
dt = 0.004
diag = 2

fy = './fy_flt_3.sav'	; by after separating and band-passing
fz = './fz_flt_3.sav'

field_width = 257	; width of field array
finfo = FILE_INFO(fy)
fsize = finfo.size
field_length = fsize / (field_width * 4 + 8)	; length of field array
						; field data stored in float (4 bytes) + 2 headers (8 bytes)
field_y = fltarr(field_width, field_length)
field_z = fltarr(field_width, field_length)

openr, 1, fy  &  point_lun, 1, 0
openr, 2, fz  &  point_lun, 2, 0

tmp1 = fltarr(field_width)
tmp2 = fltarr(field_width)
head = var.head

for i=0,field_length-1 do begin

	readu, 1, head  &  readu, 1, tmp1  &  readu, 1, head 
	readu, 2, head  &  readu, 2, tmp2  &  readu, 2, head 

	field_y[*, i] = tmp1  &  field_z[*, i] = tmp2

endfor

close, 1  &  close, 2

wpia_position = [4096, 8192, 12288, 16384, 20480, 24576, 28672]
x_offset = wpia_position[wpia]

particle = './particles0.sav'
data = dblarr(5)	; t, vx, vy, vz, x
ddt = dt * diag

openr, 1, particle  &  point_lun, 1, 0

; 2d histogram parameters --------------------------
v_bin = 0.25   &  z_bin = !dpi / 50.0	; bin width
nbin1 = 100    &  nbin2 = 101	; should change when bin width changes
step1 = 0      &  step2 = 0  &  loop = 0L
v_para = fltarr(33554432)    &   p_zeta = fltarr(33554432)
pmax = 1  &  pmin = -1
average_number = 1	; should correlate with gyrofrequency
hist2d_average = dblarr(nbin1, nbin2, average_number)
hist2d_initial = dblarr(nbin1, nbin2)
hist2d_tmp = dblarr(nbin1, nbin2)

; vpara correction --------------------------------
wc = 1.0  &  wp1 = 15.0  &  wp2 = 0.3  &  ww = 0.05  &  cv = 100
wp   = sqrt(wp1^2.0 + wp2^2.0)
xi2  = abs(ww * (wc - ww) / wp^2.0)
chi2 = 1.0 / (1.0 + xi2)
vph  = cv * sqrt(xi2) * sqrt(chi2)
vr_term_a = cv^2.0 + vph^2.0 * wc^2.0 / ww^2.0
vr_term_b = -cv * vph
vr_term_c = vph^2.0 * (1.0 - wc^2.0 / ww^2.0)
vr_zero = (-vr_term_b - sqrt(vr_term_b^2.0 - vr_term_a * vr_term_c)) / vr_term_a * cv

while(NOT EOF(1)) do begin

	readu, 1, data
	step1 = round(data[0] / ddt)
	if ((step1 gt step2)  &&  (step2 ne 0)) then begin

		truncate_part = where(v_para eq 0)  &   truncate_posi = truncate_part[0]

		hist2d = HIST_2D(v_para[0:truncate_posi], p_zeta[0:truncate_posi],  $
			bin1=v_bin, bin2=z_bin, max1=0, max2=2*!dpi, min1=-25, min2=0)
		hist2d = transpose(hist2d)
	
		hist2d_average[*, *, (step2-1) mod average_number] = hist2d	

		if (step2 eq average_number) then begin
			hist2d_initial = 0.0
			for i=0,average_number-1 do begin
				hist2d_initial = hist2d_initial + hist2d_average[*,*,i]
			endfor
			hist2d_initial = hist2d_initial / average_number
			hist2d_initial[where(hist2d_initial eq 0)] = 1d-30
		endif

		if (step2 gt average_number) then begin
			hist2d_tmp = 0.0
			for i=0,average_number-1 do begin
				hist2d_tmp = hist2d_tmp + hist2d_average[*,*,i]
			endfor
			hist2d_tmp = hist2d_tmp / average_number
			vzeta = (hist2d_tmp / hist2d_initial) - 1.0
			;vzeta = alog10(hist2d_tmp / hist2d_initial)
			vzeta[where(vzeta gt pmax)] = 0.95 * pmax 

			CGIMAGE, vzeta, POSITION=[0.2, 0.1, 0.8, 0.9],                  $
				/axes, xrange=[0, 2*!pi], yrange=[-25, 0],              $
			        maxvalue=pmax, minvalue=pmin, /interpolate	
			CGCOLORBAR, /vertical, /right, POSITION=[0.91, 0.1, 0.92, 0.9], $
				range = [pmin, pmax]
		endif	

		v_para[*] = 0.0  &   p_zeta[*] = 0.0  &  loop = 0L
		
	endif
	step2 = step1

	vx = data[1]  &  vy = data[2]  &  vz = data[3]  &  xe = data[4]

	b_pos = xe - x_offset + (field_width-1) / 2.0
	sf1 = b_pos - floor(b_pos)
	sf2 = 1.0 - sf1

	by  = field_y[floor(b_pos), step2] * sf2 + field_y[floor(b_pos)+1, step2] * sf1
	bz  = field_z[floor(b_pos), step2] * sf2 + field_z[floor(b_pos)+1, step2] * sf1

	; compute zeta
	vdotb = vy * by + vz * bz
	vcrsb = vy * bz - vz * by

	absb  = max([sqrt(by^2.0 + bz^2.0), 1d-30])
	absv  = max([sqrt(vy^2.0 + vz^2.0), 1d-30])

	cos_vb = vdotb / (absv*absb)
	cos_vb = max([cos_vb, -1d0])
	cos_vb = min([cos_vb,  1d0])

	zeta = acos(cos_vb)
	if (vcrsb lt 0) then begin
		zeta = -zeta + 2 * !dpi
	endif

	; vpara correction
	gamm = 1.0 / sqrt(1.0 - (vx^2.0 + vy^2.0 + vz^2.0) / cv^2)
	vr = abs((1.0 - wc/(gamm * ww)) * vph)
	vdiff = vx + vr
	vpara = vr_zero + vdiff

	v_para[loop] = vpara  &  p_zeta[loop] = zeta
	;v_para[loop] = vx  &  p_zeta[loop] = zeta

	loop = loop + 1

endwhile

close, 1

END	; end program
