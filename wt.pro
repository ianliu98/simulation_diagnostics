;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;              Frequency - time plot
;
;    separate into Forward and Backward waves,
;
;    The data obtained by the kempo code is used
;    Nov 27, 2009  M. Hikishima
;
;
; - Frequency resolution
;   Sampling period:    ts = dt * ifdiag [Omega_e^-1] 
;   Sampling frequency: ws = 2*pi* fs
;                          = 2*pi* (1/ts)
;                          = 2*pi/ (dt * ifdiag) [Omega_e] 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

LOADCT, 39

wave='hiss'

;           df      FFT len

; hiss:   10-30Hz  50-100msec


; w_max = 2 pi / dt * ifdiag


; he,hz
omg_max 	= 0.3
fft_pnt		= 2048L    ;modified later
tskip = 1
fc = 13000 ; [Hz]
tslide = 32



;    fft_pnt = fft_pnt/ dataSkip
;------------------------------------

; const
log2 = Alog10(2.0)

; Option
     ; use window function
     wind = 2
     ;read, 'Use window function ? rect[1], Hanning[2] : ', wind
     ; add noise
     unoise = 0
     ; add ECH and MSW
     addwave = 0
     if (addwave eq 1) then  print, 'If add ECH, select E field plot'
;------------------------------------
  
; Setting of output option
ps = 1
; READ, 'Make PS file ?, y[1], n[2] : ', ps
                                                                                                       
; Get global variables
  GLBVAR, var

; Read files
  INPUT_FILE, jobname, prefname, firnum, endnum, njob

; Set output option
eb = 2
;  Read, 'Ew[1] or Bw[2] ? : ', eb

;-----------------------------------------------------------------------------
; Reading of file
once = 0
jjt  = 0L
jst  = 0L
jt2  = 0L

FOR ijob = 0, njob-1 do begin

  datjob = '../../dat/' + prefname + '/' + jobname[ijob]

  ; Open parameter file
  prm_file = datjob + '.prm'
  READ_KEMPOPRM, prm_file, kmp          ; parameters of Kempo

  ; Open data file
  if (eb eq 1) then begin
    file1 = datjob + '.ey'
    file2 = datjob + '.ez'
    coe_e = 1.0/ kmp.cv
  endif
  if (eb eq 2) then begin
    file1 = datjob + '.by'
    file2 = datjob + '.bz'
  endif

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum 
  ntstep = (kmp.nstep - nprev) /kmp.ifdiag ; The total number of time step
;  ntstep = ntstep / 2.

  ; call only once
  if (once eq 0) then begin
    total_tstep = ntstep * njob/ tskip
    time        = Fltarr(total_tstep)
    start_time = (nprev + kmp.ifdiag) * kmp.dt
    
    ;fft_pnt = kmp.nstep/10240.

    ; setting of posiotion to plot
    ix_g = Lonarr(5)
    ix_g[0] = kmp.nx/2 - 100 * kmp.cv
    ix_g[1] = kmp.nx/2 - 50 * kmp.cv
    ix_g[2] = kmp.nx/2
    ix_g[3] = kmp.nx/2 + 50 * kmp.cv
    ix_g[4] = kmp.nx/2 + 100 * kmp.cv
    px = N_Elements(ix_g)

    dts = kmp.dt * kmp.ifdiag * tskip   ; Sampling time [Omega_e^-1]
    ws = 2.d0 * !dpi/ dts   ; Sampling freq. (angular) [Omega_e]
    dw = ws/ fft_pnt ; Frequency resolution (angular) [Omega_e^-1]
    omg_ar  = Findgen(fft_pnt) * dw - dw/2.0
    ; Less than Nyquist
    IF (ws/2 LE omg_max) THEN BEGIN
      omg_max = ws/ 2
    ENDIF
    jwmax  = Max( Where(omg_ar le omg_max) )
    omg_ar = omg_ar[0:jwmax]
    nfreq = N_Elements(omg_ar)

    nt  = Long((total_tstep - fft_pnt)/ tslide) + 1
    time_fft = Fltarr(nt)
    wt       = Fltarr(nfreq, nt, px)
    wt_f     = wt[*,*,*]
    wt_b     = wt_f[*,*,*]
    print, 'Array size: ', nfreq, nt, px

      ;--- Definition of FFT parameters ---
      ; below an interesting frequency
;     FFTts = 2048     ; size of time array
;      FFTts = 512     ; size of time array
;     FFTts = 256     ; size of time array
;      tsliede = Fix( (total_tstep - fft_pnt)/ FFTts )


    ; Frequency array
    if (total_tstep le fft_pnt) then begin
       print, 'FFT point < total time ?'  &  stop
    endif

    c_we0 = '['+var.lomgc+ '!S!De0!N]'
    c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'

    print, ''
    if (ws le 3*kmp.wc) then begin
       print, '!!!   sampling fequency is too small   !!! : ', ws
       stop
    endif
    print, '>>>'
    print, '   Sampling fequency:   ' +String( ws,'(f5.2)') +' [Omega_e0]'
    print, '   dw:                    '  +String( dw,'(f7.5)') +' [Omega_e0]'
    print, '   total time:       ' +String( total_tstep,'(I8.0)') +' [Omega_e0^-1]'
    print, '   FFT position (gird):', ix_g[*]

    fy   = Fltarr(px, total_tstep)
    fz   = fy[*,*]
    fy_f = fy[*,*] & fz_f = fy[*,*]
    fy_b = fy[*,*] & fz_b = fy[*,*]

    head  = var.head
    fytmp = Fltarr(kmp.nx)
    fztmp = fytmp[*]
    ttmp  = FLTARR(1)
    x     = (Findgen(kmp.nx) - kmp.nx/2.) * kmp.dr/ kmp.cv
    once = 1
  endif 

  ; Open data file
  Openr, 2, file1   &   point_lun, 2, 0  &  print, ' '  &  print, 'Opening ', file1
  Openr, 3, file2   &   point_lun, 3, 0  &  print, 'Opening ', file2  &  print,' '

  ; Reading of data
  FOR jt=0L, ntstep-1 do begin
    Readu, 2, head  &  Readu, 2, ttmp   &  Readu, 2, head
    Readu, 2, head  &  Readu, 2, fytmp  &  Readu, 2, head
    Readu, 3, head  &  Readu, 3, ttmp   &  Readu, 3, head
    Readu, 3, head  &  Readu, 3, fztmp  &  Readu, 3, head

    if (jst mod tskip) eq 0 then begin
    
      time[ jjt ] = ttmp
      fy[*, jjt]  = fytmp[ ix_g[*] ] 
      fz[*, jjt]  = fztmp[ ix_g[*] ]
  
      divf = DIVFIELD( fytmp[*], fztmp[*] )
      ; for forward
      fy_f[*, jjt]	= divf[ ix_g[*], 0 ]
      fz_f[*, jjt]	= divf[ ix_g[*], 1 ]
      ; for backward
      fy_b[*, jjt]	= divf[ ix_g[*], 2 ]
      fz_b[*, jjt]	= divf[ ix_g[*], 3 ]
  
      jjt++

    endif

    jst++
    if ((jst Mod 40000) eq 0) then  print, 'Reading ..., time : ',jst

    ENDFOR
    Close, 1  &  Close, 2  &  Close, 3
  ENDFOR     ; job

; ----------------------------------------------------------------------------
;                                 Plot
; ----------------------------------------------------------------------------
imgdir = './IMG/wt.eps'
SET_DISP, ps, imgdir

!P.multi   = [0, 5, 3]
topd = !D.table_size-2

; ----------------------------------------------------------------------------
; Setting of FFT

; make window function
if (wind ne 1) then  wind = Hanning(fft_pnt, alpha=0.5)     ; 0.5 is Hanning

;------------------------------------------------------------------------------
; FFT

fy   = Transpose(fy)
fz   = Transpose(fz)
fy_f = Transpose(fy_f)
fz_f = Transpose(fz_f)
fy_b = Transpose(fy_b)
fz_b = Transpose(fz_b)


a = systime(/seconds)

for jx = 0, px-1 do begin
  jt1 = 0L
  for jt = 0L, nt-1 do begin

    jt2 = jt1 + fft_pnt - 1
    time_fft[jt] = 0.5 * ( time[jt1] + time[jt2] )

    tmp1 = FFT( fy[jt1:jt2, jx]*wind, -1 )
    tmp2 = FFT( fz[jt1:jt2, jx]*wind, -1 )
    tmp1 = tmp1[0:jwmax]
    tmp2 = tmp2[0:jwmax]
    wt[*,jt,jx] = log2 + 0.5 * Alog10(  REAL_PART(tmp1)^2 + IMAGINARY(tmp1)^2 $
                                      + REAL_PART(tmp2)^2 + IMAGINARY(tmp2)^2 )

    tmp1 = FFT( fy_f[jt1:jt2, jx]*wind, -1 )
    tmp2 = FFT( fz_f[jt1:jt2, jx]*wind, -1 )
    tmp1 = tmp1[0:jwmax]
    tmp2 = tmp2[0:jwmax]
    wt_f[*,jt,jx] = log2 + 0.5 * Alog10(  REAL_PART(tmp1)^2 + IMAGINARY(tmp1)^2 $
                                        + REAL_PART(tmp2)^2 + IMAGINARY(tmp2)^2 )

    tmp1 = FFT( fy_b[jt1:jt2, jx]*wind, -1 )
    tmp2 = FFT( fz_b[jt1:jt2, jx]*wind, -1 )
    tmp1 = tmp1[0:jwmax]
    tmp2 = tmp2[0:jwmax]
    wt_b[*,jt,jx] = log2 + 0.5 * Alog10(  REAL_PART(tmp1)^2 + IMAGINARY(tmp1)^2 $
                                        + REAL_PART(tmp2)^2 + IMAGINARY(tmp2)^2 )
    jt1 = jt1 + tslide
  endfor
endfor     ; fft

print, 'dd', systime(/seconds) - a

; case of E field
if (eb eq 1) then begin
logcoe = Alog10(coe_e)
wt   = wt[*,*,*]   + logcoe
wt_f = wt_f[*,*,*] + logcoe
wt_b = wt_b[*,*,*] + logcoe
endif

; transpose
wt   = Transpose( wt[*,*,*], [1,0,2] )
wt_f = Transpose( wt_f[*,*,*], [1,0,2] )
wt_b = Transpose( wt_b[*,*,*], [1,0,2] )

; Remove the DC
wt = wt[*, 1:*, *]
wt_f = wt_f[*, 1:*, *]
wt_b = wt_b[*, 1:*, *]
omg_ar = omg_ar[1:*]

; excluding a DC component
plmax = MAX( [ Max(wt[*,1:*,*]), Max(wt_f[*,1:*,*]), Max(wt_b[*,1:*,*])] )
print, 'Data max : ', plmax

; plot range setting
plmax = -4
plmin = -7

print, 'Plot max, min : ', plmax, plmin
print, ''

; Convert the Power spectrum to V^2/Hz
me = 9.10938356e-31
qe = 1.602e-19
b0eq = 2.*!pi*fc * me/ qe ; wc = qB/m
df = dw/ 2./!pi*fc
print, df
; nT^2/ Hz
plmax_pow2 = (10.^plmax * b0eq * 1e9)^2 / df 
plmin_pow2 = (10.^plmin * b0eq * 1e9)^2 / df 

print, 'V^2/Hz(max, min): ', plmax_pow2, plmin_pow2


; ----------------------------------------------------------------------------
;                                 Plot
; ----------------------------------------------------------------------------
imgdir = './IMG/wt.eps'
ps = 0
SET_DISP, ps, imgdir

topd = !D.table_size-2
multimargin = 1
;if (!d.name eq 'PS') then begin
  !X.omargin = [18, 15]
  !Y.omargin = [10, 20]
;endif else begin
;  !p.font = 4
;  !X.omargin = [18, 2]
;  !Y.omargin = [10, 25]
;endelse

; scale factor
sf = 3
if (time_fft[nt/2] ge 1e4) then sf = 4
if (time_fft[nt/2] ge 1e5) then sf = 5

time_fft = time_fft/ 10.^sf
dt2 = (tslide * dts)/ 2.0
dt2 = dt2/10.^sf

; common
axes_a = { xrange: [time_fft[0]-dt2, time_fft[nt-1]+dt2], $
           yrange: [omg_ar[0], omg_ar[ n_elements(omg_ar)-1]], $
           xstyle: 1, ystyle: 1, xminor: 1, ticklen: -0.02, charsize: 3}

; Upper panels (not separation forward) -------------------------------
for jx=0, px-1 do begin
  axes=axes_a
  axes=Create_struct(axes, 'xtickformat', "(A1)" )
  if (jx ge 1) then begin 
    axes=Create_struct(axes, 'ytickformat', "(A1)" )
  endif
  
  cgImage, wt[*,*, jx], /axes,  $
    minvalue=plmin, maxvalue=plmax, axkeywords=axes, multimargin=multimargin,  $
    top=topd
  xp = (ix_g[jx] - kmp.nx/2)/ kmp.cv
  cgText, !X.window[0]+0.03, !Y.window[1]+0.02, charsize=1.5, color='blue', $
   'h =' +String( xp,'(I5)'), /normal
endfor
xwindow_t = !X.window  &  ywindow_t = !Y.window

; unit of h
cgText, !X.window[1]-0.02, !Y.window[1]+0.02, charsize=1.5,  $
	   + '[ c' +var.lomgc+ '!S!De0!R!U-1!N]', color='blue', /normal

; Middle panels (forward) -------------------------------------------
for jx=0, px-1 do begin
  axes=axes_a
  axes=Create_struct(axes, 'xtickformat', "(A1)" )
  if (jx ge 1) then begin
    axes=Create_struct(axes, 'ytickformat', "(A1)" )
  endif

  cgImage, wt_f[*,*, jx], /axes,  $
    minvalue=plmin, maxvalue=plmax, axkeywords=axes, multimargin=multimargin, top=topd
endfor
xwindow_m = !X.window  &  ywindow_m = !Y.window

; bottom panels (backward) -------------------------------------------
for jx=0, px-1 do begin
  axes=axes_a
  if (jx ge 1) then begin
    axes=Create_struct(axes, 'ytickformat', "(A1)")
  endif

  cgImage, wt_b[*,*, jx], /axes,  $
    minvalue=plmin, maxvalue=plmax, axkeywords=axes, multimargin=multimargin, top=topd
endfor
xwindow_b = !X.window  &  ywindow_b = !Y.window

if (eb eq 1) then  bartitle='log!D10!N E!Dw!N/cB!D0!N'
if (eb eq 2) then  bartitle='log!D10!N B!Dw!N/B!D0!N'

; colorbar
cgColorbar, /vertical, /right, $
  range=[plmin, plmax], $
  pos=[xwindow_b[1]+0.01, ywindow_b[0], xwindow_b[1]+0.02, ywindow_t[1]],  $
  ncolors=topd+1, charsize=3, title=bartitle

; Annotation
cgText, 0.45, 0.01, 't     x10!e'+string(sf,'(I1)')+'!N '+ c_we0m, /normal, charsize=2

IF (!D.name EQ 'PS') THEN BEGIN
  print, ' ' & print, 'The image is outputed to', imgdir & print, ' '
  Device, /close 
ENDIF

END
