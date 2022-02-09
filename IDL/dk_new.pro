LOADCT, 33
GLBVAR, var

; Read files
fname = ''      ; define fname as string
READ, 'main name of file ? : ', fname
READ, 'first number of files to input ? : ', firnum
READ, 'end number of files to input ? : ', endnum
njob = endnum - firnum + 1
fnum = STRING(FIX(firnum) + INDGEN(njob))
jobname = fname + fnum
jobname = STRCOMPRESS(jobname, /REMOVE_ALL)     ; delete some spaces
PRINT, ''
PRINT, 'open files ... ', jobname


fnames = ['.ex', '.ey', '.ez', '.bx', '.by', '.bz']
fnum1 = 25
fnum2 = 26

; data skip (dt = ifdiag * dskp)
;PRINT, 'default : dskp = 1'
;READ, 'time skip (1 or even) ? : ', dskp
dskp = 1

; Read files
fchk = 0
jjt = 0L
jjjt = 0L
head = BYTARR(8)     ; --- header of file ---
FOR ijob = 0, njob-1 DO BEGIN
	prmfile = '../../dat/' + fname +'/'+ jobname(ijob)+ '.prm'          ; parameter file
	fFile1 = '../../dat/' +  fname +'/'+ jobname(ijob) + fnames(fnum1 - 21)     ; data file
	fFile2 = '../../dat/' +  fname +'/'+ jobname(ijob) + fnames(fnum2 - 21)

	READ_KEMPOPRM, prmfile, kmp
	nx = kmp.nx
	dr = kmp.dr
	dt = kmp.dt
	jobnum = kmp.jobnum
	nstep = kmp.nstep
	ifdiag = kmp.ifdiag
	cv = kmp.cv

	nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
	nt = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot
	IF (fchk EQ 0) THEN BEGIN
		xcen = nx / 2
		xx = xcen	; decide region of fft
		xmin = 0
		xmax = nx - 1
		ttime = nt * njob
		print,ttime
		ret = ttime/ dskp
		start_time = (nprev + kmp.ifdiag) * kmp.dt
		field1 = complexarr(xx, ret)
		field2 = field1
		fld = Fltarr(nx)
		fld2 = fld[*]


		t = FLTARR(ret)
		tmp = FLTARR(nx)
		time = FLTARR(1)
		fchk = 1
		dk = 2.*!PI/(2*xx)/(dr/cv)
		nn = 4
		testf = sin(nn*dk*findgen(xx))
		testfft = fft(testf,-1)
	ENDIF

	PRINT, 'opening ', fFile1     ; open the file
	PRINT, 'opening ', fFile2
	OPENR, 1, fFile1
	POINT_LUN, 1, 0               ; beginning of the file
	OPENR, 2, fFile2
	POINT_LUN, 2, 0

	head = var.head
	FOR jt = 0, nt-1 DO BEGIN
		READU, 1, head  &  READU, 1, time  &  READU, 1, head
		READU, 1, head  &  READU, 1, fld   &  READU, 1, head      
		READU, 2, head  &  READU, 2, time  &  READU, 2, head
		READU, 2, head  &  READU, 2, fld2  &  READU, 2, head
		IF ((jjt MOD dskp) EQ 0) THEN BEGIN
			t(jjjt) = time
			divf = DIVFIELD(fld[*],fld2[*])
			fy_f = divf[*,0]
			fz_f = divf[*,1]
			fy_b = divf[*,2]
			fz_b = divf[*,3]
			fft_field1 = fft(fy_f[xmin:xmax-1],-1)
			fft_field2 = fft(fz_f[xmin:xmax-1],-1)
			field1(*,jjjt) = 2. * fft_field1[0:xx-1]
			field2(*,jjjt) = 2. * fft_field2[0:xx-1]
			jjjt = jjjt + 1
			ENDIF
			jjt = jjt + 1
	ENDFOR
	CLOSE, 1  &  CLOSE, 2  &  CLOSE, 3
ENDFOR     ; end of job loop
 
field1 = abs(field1)
field2 = abs(field2)

field = 2. * SQRT( field1^2 + field2^2 )
print,max(field)

x = findgen(nx)

k = FINDGEN(2*xx)
dk = 2.*!PI/ (2*xx)/ (dr/cv)
k = k(*) * dk
ku = k[0:xx-1]

; --------------------
; --- calculate dk ---
; --------------------
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

; Plot
cgwindow, wxsize=1024, wysize=768
ppos=[.1, .1, .9, .9]

; range
w_typical = 0.05
w_bot = 0.04
w_top = 0.06
grid_typical = n_elements( where(ww lt w_typical) ) - 1
grid_bot = n_elements( where(ww lt w_bot) ) - 1
grid_top = n_elements( where(ww lt w_top) ) - 1
k11 = grid_bot 
k22 = grid_top 
t11 = 0
t22 = 5000
minval = 0
maxval = max(delta_k[k11:k22, t11:t22])

xyaxes = { xrange: [ku[k11], ku[k22]], yrange: [t[t11], t[t22]]/10000.0, $
           xtitle: 'k [c!E-1!N!4X!X!De0!N]', ytitle:'t x10!U4!N[!4X!X!De0!N!E-1!N]'}

cgimage, delta_k[k11:k22,t11:t22], ku[k11:k22], t[t11:t22],  $
	maxvalue=maxval, minvalue=minval,  $
	/axes, position=ppos, /interpolate,    $
	axkeywords=xyaxes,/window

cgcolorbar, /vertical, /right, range=[minval,maxval], $
	position=[0.91,0.1,0.92,0.9], title='delta_k',/addcmd
	
end
