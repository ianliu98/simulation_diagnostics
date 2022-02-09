;
;      Make a figure as Figure 5 of Omura2015
;

LOADCT, 39

tmp=['chorus', 'hiss']
;read, 'plot for chorus[0] or hiss[1] ?: ', pf
pf = 1
wave=tmp[pf]

;           df      FFT len
; chorus:
; hiss:   10-30Hz  50-100msec
;

; w_max = 2 pi / dt * ifdiag
if wave eq 'chorus' then begin 
  fft_pnt		= 1024L
  omg_max 	= 1.0
  fc = 10000 ; [Hz]
endif else begin
  fft_pnt		= 4096L
;  omg_max 	= i
  omg_max 	= 0.3
  fc = 13000 ; [Hz]
endelse


;    fft_pnt = fft_pnt/ dataSkip
;------------------------------------

; const
log2 = Alog10(2.0)

;------------------------------------
  
;dummy = ''
;read, 'OK ?', dummy

; Setting of output option
; ps = 2
; READ, 'Make PS file ?, y[1], n[2] : ', ps
; psname = './bw.eps'

; Get global variables
  GLBVAR, var

  !P.charsize = 3
; Read files
  INPUT_FILE, jobname, prefname, firnum, endnum, njob

; Set output option
eb = 2
;  Read, 'Ew[1] or Bw[2] ? : ', eb

;-----------------------------------------------------------------------------
; Reading of file
once = 0
jjt  = 0L
jt2  = 0L

FOR ijob = 0, njob-1 do begin

;   datjob = '../dat/' + jobname[ijob]
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

  ; call only once
  if (once eq 0) then begin
    nt = ntstep * njob

;nt = nt/32
fft_pnt = nt







    time        = Fltarr(nt)
    start_time = (nprev + kmp.ifdiag) * kmp.dt
     
    ; setting of posiotion to plot;
    ix_g = kmp.nx/ 2
    ;ix_g =  kmp.nx - 3000 

    dts = kmp.dt * kmp.ifdiag   ; Sampling time [Omega_e^-1]
    ws = 2.d0 * !dpi/ dts   ; Sampling freq. (angular) [Omega_e]
    dw = ws/ fft_pnt ; Frequency resolution (angular) [Omega_e^-1]
    omg_ar  = Findgen(fft_pnt) * dw - dw/2.0
    ; Less than Nyquist
;    IF (ws/2 LE omg_max) THEN BEGIN
;      omg_max = ws/ 2
;    ENDIF
;    jwmax  = Max( Where(omg_ar le omg_max) )
;    omg_ar = omg_ar[0:jwmax]
;    nfreq = N_Elements(omg_ar)

;    nt  = Long(total_tstep)
;    time_fft = Fltarr(nt)
;    wt       = Fltarr(nfreq)
;    wt_f     = wt[*]
;    wt_b     = wt_f[*]
;    print, 'Array size: ', nfreq

    ; Frequency array
;    if (total_tstep le fft_pnt) then begin
;       print, 'FFT point < total time ?'  &  stop
;    endif

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
    print, '   total time:       ' +String( nt,'(I8.0)') +' [Omega_e0^-1]'
    print, '   FFT position (gird):', ix_g

    fy   = Fltarr(nt)
    fz   = fy[*]
    fy_f = fy[*] & fz_f = fy[*]
    fy_b = fy[*] & fz_b = fy[*]

    head  = var.head
    fytmp = Fltarr(kmp.nx)
    fztmp = fytmp[*]
    ttmp  = FLTARR(1)
    x     = (Findgen(kmp.nx) - kmp.nx/2.) * kmp.dr/ kmp.cv
    once = 1
  endif 

  ; Open data file
  Openr, 2, file1   &   point_lun, 2, 0  &  print,' '  &  print, 'Opening ', file1
  Openr, 3, file2   &   point_lun, 3, 0  &  print, 'Opening ', file2  &  print,' '

  ; Reading of data
  FOR jt=0L, ntstep-1 do begin
    Readu, 2, head  &  Readu, 2, ttmp   &  Readu, 2, head
    Readu, 2, head  &  Readu, 2, fytmp  &  Readu, 2, head
    Readu, 3, head  &  Readu, 3, ttmp   &  Readu, 3, head
    Readu, 3, head  &  Readu, 3, fztmp  &  Readu, 3, head

    time[ jjt ] = ttmp

;    fy[jjt]  = fytmp[ ix_g ] 
;    fz[jjt]  = fztmp[ ix_g ]
;    fy_f[jjt]  = fy[ jjt ] 
;    fz_f[jjt]  = fz[ jjt ] 

    divf = DIVFIELD( fytmp[*], fztmp[*] )
    ; for forward
    fy_f[jjt]	= divf[ ix_g, 0 ]
    fz_f[jjt]	= divf[ ix_g, 1 ]
    ; for backward
    fy_b[jjt]	= divf[ ix_g, 2 ]
    fz_b[jjt]	= divf[ ix_g, 3 ]

;    if (jjt eq nt-1) then   break
    jjt++
    if ((jjt Mod 40000) eq 0) then  print, 'Reading ..., time : ',jjt

    ENDFOR
    Close, 1  &  Close, 2  &  Close, 3
  ENDFOR     ; job







print, ' '
print, ' '
print, ' '
print, ' '
print, ' '
print, 'divide is skipped'
print, ' '
print, ' '
print, ' '
print, ' '
print, ' '
print, ' '
print, ' '
print, ' '
print, ' '

; ----------------------------------------------------------------------------
;                                 Plot
; ----------------------------------------------------------------------------
;imgdir = './bomg.eps'
;SET_DISP, ps, psname

;!P.multi   = [0, 5, 3]
;topd = !D.table_size-2

; ----------------------------------------------------------------------------
; Setting of FFT


; make window function
;if (wind ne 1) then  wind = Hanning(fft_pnt, alpha=0.5)     ; 0.5 is Hanning

;------------------------------------------------------------------------------
; FFT

;fy   = Transpose(fy)
;fz   = Transpose(fz)
;fy_f = Transpose(fy_f)
;fz_f = Transpose(fz_f)
;fy_b = Transpose(fy_b)
;fz_b = Transpose(fz_b)


;for jx = 0, px-1 do begin
  jt1 = 0L
;  for jt = 0L, nt-1 do begin

;    jt2 = jt1 + fft_pnt - 1

;    time_fft[jt] = 0.5 * ( time[jt1] + time[jt2] )

    ; Energies are doubled, because energies at (fft_pnt/2 - fft_pnt),
    ;    corresponds to -w aren't considered (not plotted)
    ; -----
    ; The original calculation below is log10( sqrt( |2xfy*|^2 + |2xfz*|^2 ) )

; f = 0 - 0.3 fc (df = 0.01fc)
jf1 = 1L ; 0 is DC
;windlen = 128
windlen = nt/640







print, 'dw of pass filter: ', windlen * dw

; for concatenating
omg_plot = fltarr(1,2)

; FFT
fy_fft =  FFT( fy_f, -1 )
fz_fft =  FFT( fz_f, -1 )

; spectrum
;fdis = 2. * sqrt( abs(fy_fft)^2 + abs(fz_fft)^2 )

for jf = 0, 1000 do begin

    tmp1 = fy_fft
    tmp2 = fz_fft

    ; Bandpass at +,- freq.
    ; dc       n/2     n-1    (n=8)
    ; 0  1 2 3  4  5 6  7

    ; set 0 to DC and remain positive, nagative frequency components
    tmp1[0] = 0
    tmp2[0] = 0
    tmp1[0:jf1-1] = 0
    tmp2[0:jf1-1] = 0
    ; negative freq = 0
    tmp1[jf1+windlen: *] = 0
    tmp2[jf1+windlen: *] = 0
;    tmp1[jf1+windlen: nt-jf1-windlen] = 0
;    tmp2[jf1+windlen: nt-jf1-windlen] = 0
    if (jf1 ne 1) then begin
      tmp1[nt-jf1+1:*] = 0
      tmp2[nt-jf1+1:*] = 0
    end

    ; IFFT
    wave_fy = 2. *Real_part( FFT( tmp1, 1) )
    wave_fz = 2. *Real_part( FFT( tmp2, 1) )

    ; calculate instantaneous frequeny and its amplitude
    omg = fre_cal( wave_fy[*], wave_fz[*], dts, /amplitude )

    ; get data within the pass frequency
    ind = Where( (omg[*,0] ge omg_ar[jf1]) and (omg[*,0] lt omg_ar[jf1+windlen]))

;    ; reduce data points
    ind = congrid( ind, 3000 )


    ; array concatenation
    tmp = omg[ind,*]
    omg_plot = [omg_plot, tmp]

    print, "omg pass: ", jf1, omg_ar[jf1], omg_ar[jf1+windlen], omg_ar[jf1+windlen]-omg_ar[jf1]

    jf1 = jf1 + windlen
    if (omg_ar[jf1+windlen] gt omg_max) then     break
      
endfor


omg_plot = omg_plot[1:*,*]

;!p.multi = [0,1,2]
!p.charsize=2  
c_we0 = '['+var.lomgc+ '!S!De0!N]'
c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'

ymax = 1e-3
ymin = 1e-8

omg_plot1 = omg_plot[*,1]
j = where( omg_plot1 ge ymin)
omg_plot1 = omg_plot[j,1]
omg_plot0 = omg_plot[j,0]
omg_plot = [ [omg_plot0], [omg_plot1] ]
help,omg_plot0
leng = n_elements(omg_plot0)

;rebin
omg_plot0 = congrid( omg_plot[*,0], 65536)
omg_plot1 = congrid( omg_plot[*,1], 65536)

cgWindow, woxmargin=[5,5], woymargin=[5,5], wxsize=1024, wysize=960
cgscatter2D, omg_plot0, omg_plot1, $
  xrange=[0, omg_max], yrange=[ymin, ymax], $
   xst=1, yst=1, psym=3, fit=0,/ylog,xtitle=var.omgc+ ' '+ c_we0, $
   ytitle='B!Dw!N [B!D0!N]', /window

;width = nt/windlen
width = 640. / (leng/65536)
ii = 0
ww = FIX(65536/width)
topline = fltarr(ww,2)
for jj=0,ww-1 do begin

  seg0 = omg_plot0[ii:ii+width-1]
  seg1 = omg_plot1[ii:ii+width-1]
  ttmp = where(seg1 EQ max(seg1))
  topline[jj,*] = [seg0[ttmp[0]],seg1[ttmp[0]]]
  ii = ii+width

endfor

;cgplot,topline[*,0],topline[*,1],/yolog,/overplot,color='red'
 
;;;;;;;;;;;;;;;
;file = 'iamm.csv'
;data = Read_csv( file )
;freq = data.field1[*]
;bopt = data.field2[*]
;bth = data.field3[*]
; bopt
;cgplot, freq, bopt, color='blu7', xrange=[0, omg_max], $
;  yrange=[ymin, ymax], thick=3,  position =[0.2,0.2,0.9,0.55], $
;  xtitle= var.omgc +' ' + c_we0, ytitle='B!Dw!N [B!D0!N]', /ylog
; bth
;cgplot, freq, bth, color='red7', thick=3, /overplot


;IF (!D.name EQ 'PS') THEN BEGIN
;  print, ' ' & print, 'The image is outputed to', imgdir & print, ' '
;  Device, /close 
;ENDIF

END
