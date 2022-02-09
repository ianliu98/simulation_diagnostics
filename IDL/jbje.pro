;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   JB / JE 
;      option: bandpass
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


LOADCT, 39

GLBVAR, var

bandpass = 0
log_plt  = 1

tskip = 8
xskip = 8
diag  = 256
dt    = 0.004

INPUT_FILE, jobname, prefname, firnum, endnum, njob

FOR ijob = 0, njob-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  in job loop  ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	datjob   = '../../dat/' + prefname + '/' + jobname[ijob]
	jbf = datjob + '.j2b'
	jef = datjob + '.j2e'
	byf = datjob + '.by'
	bzf = datjob + '.bz'
	jbbwf = datjob + '.jbbw2'
	prm_file = datjob + '.prm'
	READ_KEMPOPRM, prm_file, kmp 
	nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
	ntstep = (kmp.nstep - nprev) /diag / tskip

	if (ijob eq 0) then begin
	;+++++++++++++++++++++++++++++++++++++++
	;++++++++++ preset for arrays ++++++++++
	;+++++++++++++++++++++++++++++++++++++++
		nnx	       = kmp.nx/2 - kmp.nxl
		total_tstep    = ntstep * njob
		time           = Fltarr(total_tstep)
		start_time     = (nprev + diag) * kmp.dt
		container_je   = Fltarr((kmp.nx-2*kmp.nxl)/xskip,total_tstep)
		container_jb   = Fltarr((kmp.nx-2*kmp.nxl)/xskip,total_tstep)
		container_jbbw = Fltarr((kmp.nx-2*kmp.nxl)/xskip,total_tstep)
		container_by   = Fltarr((kmp.nx-2*kmp.nxl)/xskip,total_tstep)
		container_bz   = Fltarr((kmp.nx-2*kmp.nxl)/xskip,total_tstep)
		container_jbkt = Complexarr(nnx,total_tstep)
		container_jekt = Complexarr(nnx,total_tstep)
		container_jbbwkt = Complexarr(nnx,total_tstep)
		jb             = container_jb[*,*]
		je             = container_jb[*,*]
		jbbw           = container_jb[*,*]
		forfy          = container_jb[*,*]
		forfz          = container_jb[*,*]
		jbbw_kt	       = Complexarr(nnx,total_tstep)
		time_tmp       = Fltarr(1)
		jb_tmp         = Fltarr(kmp.nx)
		je_tmp         = Fltarr(kmp.nx)
		jbbw_tmp       = Fltarr(kmp.nx)
		by_tmp         = Fltarr(kmp.nx)
		bz_tmp         = Fltarr(kmp.nx)
                jx             = Lindgen((kmp.nx-2*kmp.nxr)/xskip, increment=xskip, start=kmp.nxl)
		jjt            = 0L
	endif
	
	print, 'Opening ', jbf     &   Openr, 1, jbf     &   Point_lun, 1, 0
	print, 'Opening ', jef     &   Openr, 2, jef     &   Point_lun, 2, 0
	print, 'Opening ', jbbwf   &   Openr, 3, jbbwf   &   Point_lun, 3, 0
	print, 'Opening ', byf     &   Openr, 4, byf     &   Point_lun, 4, 0
	print, 'Opening ', bzf     &   Openr, 5, bzf     &   Point_lun, 5, 0
	head  = var.head

	for jt=0L, ntstep-1 do begin
		readU, 1, head  &  readU, 1, time_tmp     &  readU, 1, head
		readU, 1, head  &  readU, 1, jb_tmp       &  readU, 1, head
		readU, 2, head  &  readU, 2, time_tmp     &  readU, 2, head
		readU, 2, head  &  readU, 2, je_tmp       &  readU, 2, head
		readU, 3, head  &  readU, 3, time_tmp     &  readU, 3, head
		readU, 3, head  &  readU, 3, jbbw_tmp     &  readU, 3, head
		readU, 4, head  &  readU, 4, time_tmp     &  readU, 4, head
		readU, 4, head  &  readU, 4, by_tmp       &  readU, 4, head
		readU, 5, head  &  readU, 5, time_tmp     &  readU, 5, head
		readU, 5, head  &  readU, 5, bz_tmp       &  readU, 5, head
		time[jjt]  = time_tmp
		container_jb[*, jjt]   = jb_tmp[ jx[*] ]
		container_je[*, jjt]   = je_tmp[ jx[*] ]
		container_jbbw[*, jjt] = jbbw_tmp[ jx[*] ]
		divf = DIVFIELD(by_tmp[*], bz_tmp[*])
		fftmp = double(divf[*,0])
		container_by[*, jjt]   = fftmp[ jx[*] ]
		fftmp = double(divf[*,1])
		container_bz[*, jjt]   = fftmp[ jx[*] ]
		fft_jb   = fft(jb_tmp[kmp.nxl:kmp.nx-kmp.nxr-1], -1)
		fft_je   = fft(je_tmp[kmp.nxl:kmp.nx-kmp.nxr-1], -1)
		fft_jbbw = fft(jbbw_tmp[kmp.nxl:kmp.nx-kmp.nxr-1], -1)
		container_jbkt[*, jjt]   = 2. * fft_jb[0:nnx-1]
		container_jekt[*, jjt]   = 2. * fft_je[0:nnx-1]
		container_jbbwkt[*, jjt] = 2. * fft_jbbw[0:nnx-1]
		jjt = jjt + 1
	endfor

	close, 1  &  close, 2  &  close, 3  & close, 4  &  close, 5 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  end job loop  ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; band-pass filter

dts = kmp.dt * diag * tskip   ; Sampling time [Omega_e^-1]
ws = 2.d * !dpi/ dts         ; Sampling freq. (angular) [Omega_e]

fft_pnt_f = jjt
dw_f = ws/ fft_pnt_f ; Frequency resolution (angular) [Omega_e^-1]
omg_ar_f  = Findgen(fft_pnt_f) * dw_f - dw_f/2.0

minf = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]
maxf = [0.02, 0.04, 0.06, 0.08, 0.10, 0.30]

jwmin = Max( Where(omg_ar_f lt minf[2]) )
jwmax = Min( Where(omg_ar_f ge maxf[2]) )

print, 'jfmin, jfmax: ', jwmin, jwmax

ssize = size(container_jb)
xlen = ssize[1]

print, 'FFT...'

for ix=0, xlen-1 do begin

	tmp_jb   = FFT( container_jb[ix,*], -1)
	tmp_je   = FFT( container_je[ix,*], -1)
	tmp_jbbw = FFT( container_jbbw[ix,*], -1)
	tmp_by   = FFT( container_by[ix,*], -1)
	tmp_bz   = FFT( container_bz[ix,*], -1)

	tmp_jb[0:jwmin] = 0    &  tmp_jb[jwmax:*] = 0
	tmp_je[0:jwmin] = 0    &  tmp_je[jwmax:*] = 0
	tmp_jbbw[0:jwmin] = 0  &  tmp_jbbw[jwmax:*] = 0
	tmp_by[0:jwmin] = 0    &  tmp_by[jwmax:*] = 0
	tmp_bz[0:jwmin] = 0    &  tmp_bz[jwmax:*] = 0

	jb[ix, *]    = 2.d * real_part( FFT(tmp_jb, 1) )
	je[ix, *]    = 2.d * real_part( FFT(tmp_je, 1) )
	jbbw[ix, *]  = 2.d * real_part( FFT(tmp_jbbw, 1) )
	forfy[ix, *] = 2.d * real_part( FFT(tmp_by, 1) )
	forfz[ix, *] = 2.d * real_part( FFT(tmp_bz, 1) )

endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

field = sqrt(forfy^2 + forfz^2)
jbbw_2 = jb / field

jbje = jb / je

jb_kt = 2. * abs(container_jbkt)
je_kt = 2. * abs(container_jekt)
jbbw_kt = 2. * abs(container_jbbwkt)

k = FINDGEN(2*nnx)
dk = 2.*!PI/ (2*nnx)/ (kmp.dr/kmp.cv)
k = k[*] * dk
ku = k[0:nnx-1]

;---------------------------------------
;------------    plot   ----------------
;---------------------------------------

axis_format = {xrange: [-(kmp.nx-kmp.nxl)/200, (kmp.nx-kmp.nxr)/200], $
	       yrange: [0, total_tstep*diag*tskip*dt]/10000.0, $
	       xtitle: 'h [c' +var.lomgc+ '!S!De0!R!U-1!N]', $
	       ytitle: 't     x10!U4!N[' +var.lomgc+ '!S!De0!R!U-1!N]'}

cgWindow, WXSize=720, WYSize=720

maxv = 1
minv = 0
cgImage, abs(jbje), position=[.1, .1, .9, .9], minvalue=minv, maxvalue=maxv, $
	/axes, axkeywords=axis_format, /interpolate, /addcmd

cgColorbar, /vertical, /right, yminor=1, range=[minv,maxv],  charsize=1.5, title='jb/je',$
	    position=[.91, .1, .92, .9], /addcmd

kk1 = 100  &  kk2 = 400
tt1 = 0    &  tt2 = 2000
minkt = 1e-7  &  maxkt = 4e-7
axis_format2 = {xrange: [ku[kk1], ku[kk2]], $
	        yrange: [time[tt1], time[tt2]]/10000.0, $
	        xtitle: 'k', $
	        ytitle: 't     x10!U4!N[' +var.lomgc+ '!S!De0!R!U-1!N]'}
cgWindow, WXSize=720, WYSize=720
cgImage, jb_kt[kk1:kk2, tt1:tt2], ku[kk1:kk2], time[tt1:tt2], minvalue=minkt, maxvalue=maxkt, $
	/axes, axkeywords=axis_format2, pos=[.1, .1, .9, .9],/interpolate, /addcmd
cgColorbar, /vertical, /right, yminor=1, range=[minkt,maxkt],  charsize=1.5, title='jb',$
	    position=[.91, .1, .92, .9], /addcmd
end
