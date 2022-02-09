;----------------------------------------------------------------------------------
;
;                         Time - Space plot of field data
;
; OPTION:
;   - Separate into Forward and Backward waves,
;   - Distribution of wave frequency 
;       zero crossing method is used by a routine. 
;
; The data obtained by the kempo code is used
; Nov 27, 2009  M. Hikishima
;
;----------------------------------------------------------------------------------

; Shrink an array
tskip = 8
xskip = 8
LOADCT, 39

; Global variables
     GLBVAR, var

; Read files
     INPUT_FILE, jobname, prefname, firnum, endnum, njob

; Plot option
     tmp = ['ew', 'bw']
     i=1
     eb = tmp[i]

     tmp = ['notseparation', 'separation']
     i=1
     sep_mode = tmp[i]

     freqdis = 0

; Setting of output option

    read, 'time range of prospective wavepacket: t_top = ', ttop
    read, 'time range of prospective wavepacket: t_bot = ', tbot

; ---------------------------------------------------------------------------------

var.dir = '../../dat/'
; Reading of files
jjt = 0L
FOR ijob = 0, njob-1 do begin

  ;--- Read input parameters in KEMPO from .input_idl file ---
  prm_file = var.dir +prefname +'/' +jobName(ijob) +'.prm'
  READ_KEMPOPRM, prm_file, kmp

  fname = var.dir +prefname +'/' + jobname(ijob) + ['.by', '.bz']

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum 
  ntime = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot

  IF (ijob eq 0) then begin
    total_time = ntime * njob

    ;--- Reduced size of time and space array ----------------------
    ; time
    rnt = total_time/ tskip
    t = Fltarr(rnt)
    fsample =  2.d0*!dpi/ (kmp.ifdiag * kmp.dt)
    fsample_r = fsample/ tskip
    omg = 0.3
    b0 = 1

    rnx = kmp.nx / xskip

    lambda = Wavelength( omg, kmp.wp1, kmp.wp2, kmp.cv, b0)
    if xskip ge lambda/2 then begin
      print, 'xskip is too larger'
      stop
    endif

    print, '>>>'
    print, '   Wavelength: ', lambda
    print, '   Spatial skip: ', xskip 
    print, '   Then, ' 
    print, '   The output array t-x: '
    print, total_time, kmp.nx, ' ->', rnt, rnx 
    print,''  &  print, '+++ sampling frequency +++ : ', fsample, ' Omega_e0'
    print, fsample, ' ->', fsample_r, ' Omega_e0'
    print,''  &  print,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    if (fsample_r le 3.0*Abs(kmp.wc)) then begin
      print, '!!! sampling frequency is small        !!!!!'
      stop
    endif
    start_time = (nprev + kmp.ifdiag) * kmp.dt

    tmin_f = tbot / tskip
    tmax_f = ttop / tskip

    xmin_f = kmp.nxl/xskip
    xmax_f = (kmp.nx-kmp.nxr)/xskip

    tf = fltarr((tmax_f-tmin_f))
    forfy = Fltarr( xmax_f-xmin_f, n_elements(tf))
    forfz = forfy
    pasfy = forfy
    pasfz = forfy
    rcd1 = pasfy
    rcd2 = pasfy
    rcd3 = pasfy
    rcd5 = pasfy

    head = var.head
    ttmp = Fltarr(1)
    xtmp11 = Fltarr(kmp.nx)
    xtmp22 = xtmp11[*]

    jx = Lindgen(rnx, increment=xskip, start=0)
    x_ax = ( Findgen(rnx)*xskip - (kmp.nx-1.)/ 2. ) * kmp.dr/ kmp.cv

    x_ax_f = x_ax[xmin_f:xmax_f-1]

    jnt = 0L
    jkt = 0L
    jntl = 0L
    jntr = 0L
    tmin_all = tmin_f
    tmax_all = tmax_f
  endif 

  ; Open data file
  Openr, 2, fname[0]   &   Point_lun, 2, 0     ; beginning of the file
  Openr, 3, fname[1]   &   Point_lun, 3, 0
  print, 'Opening : ', fname

  for jt = 0L, ntime-1 do begin
    Readu, 2, head   &   Readu, 2, ttmp    &   Readu, 2, head
    Readu, 2, head   &   Readu, 2, xtmp11  &   Readu, 2, head
    Readu, 3, head   &   Readu, 3, ttmp    &   Readu, 3, head
    Readu, 3, head   &   Readu, 3, xtmp22  &   Readu, 3, head

    if ( (jjt MOD tskip) eq 0 ) then begin
      if ((jkt ge tmin_all) and (jkt lt tmax_all)) then begin

          xtmp1 = xtmp11[ jx[*] ]
          xtmp2 = xtmp22[ jx[*] ]

          divf = DIVFIELD(xtmp1[*], xtmp2[*])
          divf = Float( divf[*,*]) 

          tf[jntr] = ttmp
          tmp = divf[*,0]
          forfy[*, jntr] = tmp[xmin_f : xmax_f-1]
          tmp = divf[*,1]
          forfz[*, jntr] = tmp[xmin_f : xmax_f-1]

          jntr = ++jntr
      endif
      jkt = ++jkt
    endif     ; mod

    jjt = ++jjt
    IF ((jjt Mod 10000) EQ 0) THEN  PRINT, 'Reading.., time : ',jjt, '/', total_time
  ENDFOR
  CLOSE, 1   &   CLOSE, 2   &   CLOSE, 3
ENDFOR


plmax = -3.
plmin =  plmax - 3.0

print, 'Max, Min (log plot): ', plmax, plmin
print, '' & print, '--------------------------------------------------------'

c_we0 = '['+var.lomgc+ '!S!De0!N]'
c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'


tf_rd = tf
; release the memory
xtmp11 = 0
xtmp22 = 0
tf = tf[where(tf ne 0)]

passfilt = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FFT

dts = kmp.dt * kmp.ifdiag * tskip   ; Sampling time [Omega_e^-1]
ws = 2.d0 * !dpi/ dts   ; Sampling freq. (angular) [Omega_e]

; fft_pnt_b = jntl
 fft_pnt_f = jntr
dw_f = ws/ fft_pnt_f ; Frequency resolution (angular) [Omega_e^-1]
omg_ar_f  = Findgen(fft_pnt_f) * dw_f - dw_f/2.0

minf = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]
maxf = [0.02, 0.04, 0.06, 0.08, 0.10, 0.30]
rtx = fltarr( n_elements(minf) )
rty = fltarr( n_elements(minf) )

for jf=n_elements(minf)-1, 0, -1 do begin

  jwmin = Max( Where(omg_ar_f lt minf[jf]) )
  jwmax = Min( Where(omg_ar_f ge maxf[jf]) )

  print, 'jfmin, jfmax: ', jwmin, jwmax

  ; fft
  ; arrays are replaced with complex type
  xlen = xmax_f - xmin_f

  print, 'FFT...'

  for ix=0, xlen-1 do begin

    tmpy = FFT( forfy[ix,*], -1)
    tmpz = FFT( forfz[ix,*], -1)

    tmpy[0:jwmin] = 0  &  tmpy[jwmax:*] = 0
    tmpz[0:jwmin] = 0  &  tmpz[jwmax:*] = 0

    pasfy[ix, *] = 2*real_part( FFT( tmpy, 1) )
    pasfz[ix, *] = 2*real_part( FFT( tmpz, 1) )

  endfor

  ; ------------------------------------------------------------------------------
  ;                                  Plot
  ; ------------------------------------------------------------------------------

  pasfy = alog10( sqrt( temporary(pasfy)^2 + pasfz^2 ) )

  if (jf eq 2) then begin
    rcd2 = pasfy
  endif

  ; Image plot ---------------------------------------------------------------

  axisset={ticklen: -0.02, xst:1, yst:1}
  c_we0 = '['+var.lomgc+ '!S!De0!N]'

  if (jf eq 5) then begin
    cgwindow, woxmargin = [10,20], woymargin = [10,10], wxsize = 1440, wysize=1080
    !P.charsize = 3
    mmargin = 1
    topd = !D.table_size - 2
    cgImage, pasfy, /axes, AXKEYWORDS=axisset, $
      minvalue=plmin, maxvalue=plmax, top=topd, $
      xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0, $
      xtitle='h [c' +var.lomgc+ '!S!De0!R!U-1!N]', $
      title=var.omgc +'=' +string(minf[jf], '(f4.2)') +'-' +string(maxf[jf], '(f4.2)') +'  ' +c_we0, $
      ytitle = 't     x10!U4!N  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
      multimargin=[5,7,5,7],wmulti=[0,3,2],charsize=1.5, /window
  endif else begin
    cgImage, pasfy, /axes, AXKEYWORDS=axisset, $
      minvalue=plmin, maxvalue=plmax, top=topd, $
      xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0, $
      xtitle='h [c' +var.lomgc+ '!S!De0!R!U-1!N]', $
      title=var.omgc +'=' +string(minf[jf], '(f4.2)') +'-' +string(maxf[jf], '(f4.2)') +'  ' +c_we0, $
      ytitle = 't     x10!U4!N  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
      multimargin=[5,7,5,7],charsize=1.5,/addcmd
  endelse  

endfor ; f loop


bartitle='log!D10!N B!Dw!N/B!D0!N'

cgColorbar, /vertical, /right,  $
    yminor=1, range=[plmin, plmax], title=bartitle, ncolors=topd+1,  $
    position=[0.9, 0.3, 0.91, 0.7], charsize=1.5, /addcmd

end
