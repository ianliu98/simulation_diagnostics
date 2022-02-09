
; Shrink an array
tskip = 8
xskip = 8
LOADCT, 39

; Global variables
GLBVAR, var

read, 'how many cases to compare: ', cmpr
read, 'time range of prospective wavepacket: t_top = ', ttop
read, 'time range of prospective wavepacket: t_bot = ', tbot

savf = 0  ; 0 -> normal;  1 -> only save by and bz at specific position
sav  = 0  ; 0 -> plot;  1 -> plot and save;  2 -> save
; ---------------------------------------------------------------------------------

var.dir = '../../dat/'

FOR cmpr_ind = 0, cmpr-1 DO BEGIN

; Read files
INPUT_FILE, jobname, prefname, firnum, endnum, njob

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
    field = pasfy
    rcd1 = pasfy
    rcd2 = pasfy
    rcd3 = pasfy
    rcd5 = pasfy

    head = var.head
    ttmp = Fltarr(1)
    xtmp11 = Fltarr(kmp.nx)
    xtmp22 = xtmp11[*]

    jx = Lindgen(rnx, increment=xskip, start=0)
    jx2 = jx - (rnx / 2.0) * xskip
    aa = (kmp.b02/1.0 - 1.0) / (kmp.nx/kmp.cv/2.0)^2.0
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

      if ((jkt ge tmin_f) and (jkt lt tmax_f)) then begin

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
  
  CLOSE, 2   &   CLOSE, 3

ENDFOR

; release the memory
xtmp11 = 0
xtmp22 = 0
tmp    = 0

tf = tf[where(tf ne 0)] ; time

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; band-pass filter

dts = kmp.dt * kmp.ifdiag * tskip   ; Sampling time [Omega_e^-1]
ws = 2.d0 * !dpi/ dts               ; Sampling freq. (angular) [Omega_e]

fft_pnt_f = jntr
dw_f = ws/ fft_pnt_f ; Frequency resolution (angular) [Omega_e^-1]
omg_ar_f  = Findgen(fft_pnt_f) * dw_f - dw_f/2.0

minf = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]
maxf = [0.02, 0.04, 0.06, 0.08, 0.10, 0.30]

jwmin = Max( Where(omg_ar_f lt minf[2]) )
jwmax = Min( Where(omg_ar_f ge maxf[2]) )

print, 'jfmin, jfmax: ', jwmin, jwmax

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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


if (savf eq 1) then begin
  savfy_tmp = './csv/center_by' + string(round(cmpr_ind)) + '.csv'
  savfz_tmp = './csv/center_bz' + string(round(cmpr_ind)) + '.csv'
  savfy_tmp = strcompress(savfy_tmp, /REMOVE_ALL)
  savfz_tmp = strcompress(savfz_tmp, /REMOVE_ALL)
  dimensions = size(forfy, /DIMENSIONS)
  write_csv, savfy_tmp, pasfy[dimensions[0]/2,*]
  write_csv, savfz_tmp, pasfz[dimensions[0]/2,*]
endif else begin


; field
field = alog10( sqrt( pasfy^2 + pasfz^2 ))

; frequency
trans_forf  = Transpose(pasfy[*,*])
trans_forf1 = Transpose(pasfz[*,*])
forf = trans_forf
dt = tskip * kmp.ifdiag * kmp.dt
for jx = 0, xmax_f-xmin_f-1 do begin
  forf[*, jx] = fre_cal(trans_forf[*,jx], trans_forf1[*,jx], dt)
  if ((jx mod 512) eq 0) Then print, 'continuing...'
endfor
forf = Transpose(forf[*,*])

; release memory
trans_forf  = 0
trans_forf1 = 0
pasfz       = 0
tmpy        = 0
tmpz        = 0


; -------------------
; --- calculate S ---
; -------------------
Q   = 0.1
wpe = kmp.wp1
wph = kmp.wp2
cv  = kmp.cv / 100.0
utpara = kmp.path2 / 100.0
utperp = kmp.peth2 / 100.0
beta   = 0.3
uperph = sqrt(!PI/2.0) * ((1.0 - beta^(1.5) / (1.0 - beta))) * utperp
vperp  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * uperph
vpara  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * utpara
gamma  = 1.0 / sqrt(1.0 - (vperp^2.0 + vpara^2.0)/cv^2.0)
Omega_e = 1.00

xi2 = abs(forf * (Omega_e - forf) / wpe^2.0)
delta2 = 1.0 / (1.0 + xi2)

w_s = forf / Omega_e
vp_s = sqrt(xi2) * sqrt(delta2)
vperp_s = vperp / cv
vr_s = (w_s^2.0 - sqrt(w_s^4.0 + (w_s^2.0 + vp_s^2.0) * (1.0 - w_s^2.0 - vperp_s^2.0))) / (w_s^2.0 / vp_s + vp_s)
vp = vp_s * cv
vr = vr_s * cv
vr2 = (1.0 - Omega_e / (gamma * forf)) * vp

vg = (cv * sqrt(xi2) / sqrt(delta2)) / (xi2 + Omega_e / (2.0 * (Omega_e - forf)))

s0 = sqrt(delta2) * vperp / (sqrt(xi2) * cv)
s1 = gamma * (1.0 - vr / vg)^2.0

s2_term1 = 1.0 / (2.0 * sqrt(xi2) * sqrt(delta2))
s2_term2 = gamma * forf * (vperp/cv)^2.0 / Omega_e
s2_term3 = (2.0 + (forf / Omega_e) * delta2 * (Omega_e - gamma * forf) / (Omega_e - forf)) * vr * vp / cv^2.0
s2 = s2_term1 * (s2_term2 - s2_term3)

dwdt = forf
size_dim = size(forf, /DIMENSIONS)
size_x = size_dim[0]
size_t = size_dim[1]
for ii=0, size_x-1 do begin
  for jj=1, size_t-1 do begin
    dwdt[ii,jj] = (forf[ii,jj] - forf[ii,jj-1]) / dt
  endfor
endfor
dwdt[*,0] = dwdt[*,1]
dbdh = 2.0 * aa * abs(jx2)

; clear memory
vg = 0
xi2 = 0
delta2 = 0
s2_term1 = 0
s2_term2 = 0
s2_term3 = 0

k1 = s1 * dwdt
s1 = 0
k2 = s2
for ii=0, size_t-1 do begin
  k2[*,ii] = 1.0 * s2[*,ii] * dbdh[kmp.nxl/xskip:-kmp.nxr/xskip-1]
endfor
s2 = 0
k3 = s0 * forf * 10^field
s0 =0

S = -1.0 * (k1 + k2) / k3


if (sav ne 2) then begin
; ------------------------------------------------------------------------------
;                                  Plot
; ------------------------------------------------------------------------------

if (cmpr_ind eq 0) then begin

  axisset={ticklen : -0.02, $
	   xstyle  : 1,     $
	   ystyle  : 1,     $
	   charsize: 1.5,   $
	   xminor  : 1}

  pfmax = max(field)
  pfmin = pfmax - 2

  psmax = 1
  psmin = 0

  c_we0 = '['+var.lomgc+ '!S!De0!N]'
  c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'

  xsize = 1440
  ysize = 960
  xmargin = [20, 20]
  ymargin = [10, 10]
  mmargin = [3, 5, 3, 5]

  edge = 200

  cgwindow, woxmargin = xmargin, woymargin = ymargin, wxsize = xsize, wysize = ysize

endif

if (cmpr_ind eq 0) then begin
  ; field
  cgImage, field[*, edge:-edge], /axes, AXKEYWORDS=axisset, minvalue=pfmin, maxvalue=pfmax,   $
    xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0,                 $
    ytitle = 't     x10!U4!N  [' +var.lomgc+ '!S!De0!R!U-1!N]',                            $
    multimargin=mmargin, wmulti=[0,3,2], layout=[3,2,1], /window

  ; s
  cgImage, abs(S[*, edge:-edge]), /axes, AXKEYWORDS=axisset, minvalue=psmin, maxvalue=psmax,  $
    xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0,                 $
    xtitle = 'h [c' +var.lomgc+ '!S!De0!R!U-1!N]',                                         $
    ytitle = 't     x10!U4!N  [' +var.lomgc+ '!S!De0!R!U-1!N]',                            $
    multimargin=mmargin, wmulti=[0,3,2], layout=[3,2,4], /addcmd

endif else begin
  ; field
  cgImage, field[*, edge:-edge], /axes, AXKEYWORDS=axisset, minvalue=pfmin, maxvalue=pfmax,   $
    xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0,                 $
    multimargin=mmargin, wmulti=[0,3,2], layout=[3,2,cmpr_ind+1], /addcmd

  ; s 
  cgImage, abs(S[*, edge:-edge]), /axes, AXKEYWORDS=axisset, minvalue=psmin, maxvalue=psmax,  $
    xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0,                 $
    xtitle = 'h [c' +var.lomgc+ '!S!De0!R!U-1!N]',                                         $
    multimargin=mmargin, wmulti=[0,3,2], layout=[3,2,cmpr_ind+4], /addcmd

endelse  

; save
if (sav eq 1) then begin
  f1 = 'field' + string(round(cmpr_ind)) + '.csv'
  f2 = 'S'  + string(round(cmpr_ind)) + '.csv'
  f1 = strcompress(f1, /REMOVE_ALL)
  f2 = strcompress(f2, /REMOVE_ALL)
  write_csv, f1, field
  write_csv, f2, S
endif

endif

if (sav eq 2) then begin
  f1 = 'field' + string(round(cmpr_ind)) + '.csv'
  f2 = 'S'  + string(round(cmpr_ind)) + '.csv'
  f1 = strcompress(f1, /REMOVE_ALL)
  f2 = strcompress(f2, /REMOVE_ALL)
  write_csv, f1, field
  write_csv, f2, S
endif

endelse  ; savf

ENDFOR

if (savf eq 0) then begin

if (sav ne 2) then begin

fbartitle = 'log!D10!N B!Dw!N/B!D0!N'
sbartitle = '|S|' 

cgColorbar, /vertical, /right,  $
  yminor=1, range=[pfmin, pfmax], title=fbartitle,  $
  position=[0.91, 0.54, 0.92, 0.86], charsize=1.5, /addcmd

cgColorbar, /vertical, /right,  $
  yminor=1, range=[psmin, psmax], title=sbartitle,  $
  position=[0.91, 0.14, 0.92, 0.46], charsize=1.5, /addcmd

endif

endif

end
