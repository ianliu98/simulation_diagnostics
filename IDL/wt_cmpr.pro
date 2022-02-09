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
omg_max 	= 0.15
fft_pnt		= 2048L    ;modified later
tskip = 1
fc = 13000 ; [Hz]
tslide = 2L

text = ['Case 1', 'Case 2', 'Case 3']
text2 = ['Case 1', 'Case 4', 'Case 5']

ocont = 0
interp = 0

;    fft_pnt = fft_pnt/ dataSkip
;------------------------------------

; const
log2 = Alog10(2.0)

; Option
     ; use window function
;     wind = 2
     ;read, 'Use window function ? rect[1], Hanning[2] : ', wind
     ; add noise
     unoise = 0
     ; add ECH and MSW
     addwave = 0
     if (addwave eq 1) then  print, 'If add ECH, select E field plot'

     pltrg = 0
     read, 'fixed colorbar range? yes[0]  no[1] : ', pltrg

     cmpr = 1
     read, 'how many cases are gonna be compared with?: ', cmpr

     if (cmpr gt 1) then begin
       read, 'gradient[1] / density[2]: ', grde
     endif
;------------------------------------
  
; Setting of output option
ps = 1
; READ, 'Make PS file ?, y[1], n[2] : ', ps
                                                                                                       
; Get global variables
  GLBVAR, var

for cmpr_ind=0, cmpr-1 do begin

;if (cmpr_ind eq 0) then begin
;  tskip = 1
;endif else begin
;  tskip = 16
;endelse

wind = 2
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
  mltstp = kmp.mltstp

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
    test_use1 = fy[*, *] & test_use2 = fy[*, *]

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
      ;fytmp = smooth(fytmp, round(kmp.nx/mltstp))
      ;fztmp = smooth(fztmp, round(kmp.nx/mltstp))
      fy[*, jjt]  = fytmp[ ix_g[*] ] 
      fz[*, jjt]  = fztmp[ ix_g[*] ]
      test_use1[*, jjt] = fytmp[ ix_g[*] ]
      test_use2[*, jjt] = fztmp[ ix_g[*] ]
  
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
; Setting of FFT

; make window function
;fft_pnt = 256L
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
    fft_tmp1 = tmp1
    fft_tmp2 = tmp2
;    tmp1[0:round(jwmax/10)] = 0L & tmp1[length(tmp1)-round(jwmax/10):length(tmp1)-1] = 0L
;    tmp2[0:round(jwmax/10)] = 0L & tmp2[length(tmp2)-round(jwmax/10):length(tmp2)-1] = 0L
    tmp1 = tmp1[0:jwmax]
    tmp2 = tmp2[0:jwmax]
    wt[*,jt,jx] = log2 + 0.5 * Alog10(  REAL_PART(tmp1)^2 + IMAGINARY(tmp1)^2 $
                                      + REAL_PART(tmp2)^2 + IMAGINARY(tmp2)^2 )

    tmp1 = FFT( fy_f[jt1:jt2, jx]*wind, -1 )
    tmp2 = FFT( fz_f[jt1:jt2, jx]*wind, -1 )
;    tmp1[0:round(jwmax/10)] = 0L & tmp1[length(tmp1)-round(jwmax/10):length(tmp1)-1] = 0L
;    tmp2[0:round(jwmax/10)] = 0L & tmp2[length(tmp2)-round(jwmax/10):length(tmp2)-1] = 0L
    tmp1 = tmp1[0:jwmax]
    tmp2 = tmp2[0:jwmax]
    wt_f[*,jt,jx] = log2 + 0.5 * Alog10(  REAL_PART(tmp1)^2 + IMAGINARY(tmp1)^2 $
                                        + REAL_PART(tmp2)^2 + IMAGINARY(tmp2)^2 )

    tmp1 = FFT( fy_b[jt1:jt2, jx]*wind, -1 )
    tmp2 = FFT( fz_b[jt1:jt2, jx]*wind, -1 )
;    tmp1[0:round(jwmax/10)] = 0L & tmp1[length(tmp1)-round(jwmax/10):length(tmp1)-1] = 0L
;    tmp2[0:round(jwmax/10)] = 0L & tmp2[length(tmp2)-round(jwmax/10):length(tmp2)-1] = 0L
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
if (pltrg eq 0) then begin
  plmax = -3.5
  plmin = -6.5
endif else begin
  plmax = plmax
  plmin = plmax - 3
endelse

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

if (cmpr_ind eq 0) then begin
  xsize = 1440
  ysize = 960
  xmargin = [20, 20]
  xmargin2 = xmargin * 4
  ymargin = [10, 10]
  ymargin2 = ymargin * 6
  cgwindow, woxmargin = xmargin, woymargin = ymargin, wxsize = xsize, wysize = ysize
  mmargin = [0.75, 1, 0.75, 1]
endif
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
           xstyle: 1, ystyle: 1, xminor: 1, ticklen: -0.02, charsize: 1.5}


;(forward) -------------------------------------------
if (cmpr_ind eq 0) then begin
  for jx=0, px-1 do begin
    axes=axes_a
    axes=Create_struct(axes, 'xtickformat', "(A1)" )
    if (jx ge 1) then begin
      axes=Create_struct(axes, 'ytickformat', "(A1)" )
    endif
    if (jx eq 0) then begin
      cgImage, wt_f[*,*, jx], /axes, minvalue=plmin, maxvalue=plmax, interpolate=interp, $
	     axkeywords=axes, multimargin=mmargin, wmulti=[0,5,cmpr], /window
    endif else begin
      cgImage, wt_f[*,*, jx], /axes, minvalue=plmin, maxvalue=plmax, interpolate=interp, $
    	     axkeywords=axes, multimargin=mmargin, /addcmd
    endelse
    if (ocont eq 1) then begin
      ncontour = 1
      contourLevels =  cgConLevels( wt_f[*,*, jx], NLevels=ncontour, MinValue=plmax-1, MaxValue=plmax )
      cgContour, wt_f[*,*, jx], /OnImage, label=0, Levels=contourLevels, Thick=1, /addcmd
    endif

    xwindow = !X.window * xsize  &  ywindow = !Y.window * ysize
    if (xwindow[0] lt xsize/2) then begin
      xwindow_nw = xwindow[0]/(xsize/2.0) * (xsize/2.0-xmargin2[0]) + xmargin2[0]
    endif else begin
      xwindow_nw  = (xwindow[0]-xsize/2.0)/(xsize/2.0) * (xsize/2.0-xmargin2[1]) + xsize/2.0
      xwindow_nw2 = (xwindow[1]-xsize/2.0)/(xsize/2.0) * (xsize/2.0-xmargin2[1]) + xsize/2.0
    endelse
      ywindow_nw = (ywindow[1]-ysize/2.0)/(ysize/2.0) * (ysize/2.0-ymargin2[1]) + ysize/2.0

    xp = (ix_g[jx] - kmp.nx/2)/ kmp.cv
    cgText, xwindow_nw/xsize+0.01, ywindow_nw/ysize+0.01, charsize=1, color='blue', $
	 'h =' +String( xp,'(I5)'), /normal, /addcmd

  endfor

  xwindow_t = xwindow_nw2/xsize  &  ywindow_t = ywindow_nw/ysize
  ; unit of h
  cgText, xwindow_t-0.04, ywindow_t+0.01, charsize=1,  $
	  + '[ c' +var.lomgc+ '!S!De0!R!U-1!N]', color='blue', /normal, /addcmd

endif else begin

  for jx=0, px-1 do begin
    axes=axes_a
    if (cmpr_ind lt cmpr-1) then begin
      axes=Create_struct(axes, 'xtickformat', "(A1)" )
    endif
    if (jx ge 1) then begin
      axes=Create_struct(axes, 'ytickformat', "(A1)")
    endif
    cgImage, wt_f[*,*, jx], /axes, minvalue=plmin, maxvalue=plmax, interpolate=interp, $
	    axkeywords=axes, multimargin=1, /addcmd
    if (ocont eq 1) then begin
      ncontour = 1
      contourLevels =  cgConLevels( wt_f[*,*, jx], NLevels=ncontour, MinValue=plmax-1, MaxValue=plmax )
      cgContour, wt_f[*,*, jx], /OnImage, label=0, Levels=contourLevels, Thick=1, /addcmd
    endif
  endfor
  xwindow = !X.window * xsize  &  ywindow = !Y.window * ysize
  xwindow_nw2 = (xwindow[1]-xsize/2.0)/(xsize/2.0) * (xsize/2.0-xmargin2[1]) + xsize/2.0
  ywindow_nw  = ywindow[0]/(ysize/2.0) * (ysize/2.0-ymargin2[0]) + ymargin2[0]
  xwindow_b = xwindow_nw2/xsize  &  ywindow_b = ywindow_nw/ysize

endelse

endfor

if (eb eq 1) then  bartitle='log!D10!N E!Dw!N/cB!D0!N'
if (eb eq 2) then  bartitle='log!D10!N B!Dw!N/B!D0!N'

; colorbar
cgColorbar, /vertical, /right, $
  range=[plmin, plmax], $
  pos=[xwindow_b + 0.05, ywindow_b, xwindow_b + 0.06, ywindow_t],  $
  charsize=1.5, title=bartitle, /addcmd

; Annotation
cgText, 0.45, 0.04, 't     x10!e'+string(sf,'(I1)')+'!N '+ c_we0m, /normal, charsize=1.5, /addcmd
cgText, 0.04, 0.45, var.omgc +' ' + c_we0, orientation=90, /normal, charsize=1.5, /addcmd

;if (grde eq 1) then begin
;  cgText, 0.03, 0.75, charsize=0.75, text[0], /normal, /addcmd
;  cgText, 0.03, 0.50, charsize=0.75, text[1], /normal, /addcmd
;  cgText, 0.03, 0.35, charsize=0.75, text[2], /normal, /addcmd
;endif else begin
;  cgText, 0.03, 0.75, charsize=0.75, text2[0], /normal, /addcmd
;  cgText, 0.03, 0.50, charsize=0.75, text2[1], /normal, /addcmd
;  cgText, 0.03, 0.35, charsize=0.75, text2[2], /normal, /addcmd
;endelse

END
