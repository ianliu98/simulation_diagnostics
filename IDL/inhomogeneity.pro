;----------------------------------------------------------------------------------
;   calculate inhomogeneity factor S
;----------------------------------------------------------------------------------

; Shrink an array
tskip = 4 
xskip = 4
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

;read, 'time range of prospective wavepacket: t_top = ', ttop
;read, 'time range of prospective wavepacket: t_bot = ', tbot

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

    ttop = floor(total_time)
    tbot = 0
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
    Readu, 2, head   &   Readu, 2, ttmp   &   Readu, 2, head
    Readu, 2, head   &   Readu, 2, xtmp11  &   Readu, 2, head
    Readu, 3, head   &   Readu, 3, ttmp   &   Readu, 3, head
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


plmax = 0.06
plmin =  0.04

print, 'Max, Min (log plot): ', plmax, plmin
print, '' & print, '--------------------------------------------------------'

c_we0 = '['+var.lomgc+ '!S!De0!N]'
c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'


tf_rd = tf
; release the memory
xtmp11 = 0
xtmp22 = 0
xtmp1 = 0
xtmp2 = 0
dicf = 0

;tb = tb[where(tb ne 0)]
tf = tf[where(tf ne 0)]

passfilt = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FFT

dts = kmp.dt * kmp.ifdiag * tskip   ; Sampling time [Omega_e^-1]
ws = 2.d0 * !dpi/ dts   ; Sampling freq. (angular) [Omega_e]

fft_pnt_f = jntr
dw_f = ws/ fft_pnt_f ; Frequency resolution (angular) [Omega_e^-1]
omg_ar_f  = Findgen(fft_pnt_f) * dw_f - dw_f/2.0

minf = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]
maxf = [0.02, 0.04, 0.06, 0.08, 0.10, 0.30]
rtx = fltarr( n_elements(minf) )
rty = fltarr( n_elements(minf) )

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
  ; ifft
  pasfy[ix, *] = 2*real_part( FFT( tmpy, 1) )
  pasfz[ix, *] = 2*real_part( FFT( tmpz, 1) )

endfor
pasfy_fre = pasfy
pasfz_fre = pasfz

; field
field = sqrt( temporary(pasfy)^2 + pasfz^2 )

; frequency
trans_forf  = Transpose(pasfy_fre)
trans_forf1 = Transpose(pasfz_fre)
forf = trans_forf
dt = tskip * kmp.ifdiag * kmp.dt

; ---------------
; method 1  ----
; --------------
;truncate = 0.025
;truncate_pnt = Min(Where( omg_ar_f ge truncate))
;for jx=0, xmax_f-xmin_f-1 do begin
;  fre_tmp = Omg_cal( trans_forf[*,jx], trans_forf1[*,jx], dt )
;  ; fft to filter frequency
;  fre_tmp_fft = FFT( fre_tmp[1:*], -1)
;  fre_tmp_fft[truncate_pnt:*] = 0
;  forf[0, jx] = 0.0
;  forf[1:*, jx] = real_part( FFT( fre_tmp_fft, 1))
;endfor

; -------------
; method 2-----
; -------------
for jx = 0, xmax_f-xmin_f-1 do begin
	forf[*, jx] = fre_cal(trans_forf[*,jx], trans_forf1[*,jx], dt)
	if ((jx mod 512) eq 0) Then print, 'continuing...'
endfor
forf = Transpose(forf[*,*])
fre = forf

; release memory
forf = 0
trans_forf = 0
trans_forf1 = 0
pasfy_fre = 0
pasfz_fre = 0

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
uperph  = sqrt(!PI/2.0) * ((1.0 - beta^(1.5) / (1.0 - beta))) * utperp
vperp  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * uperph
vpara  = cv / sqrt(cv^2.0 + (utpara^2.0 + uperph^2.0)) * utpara
gamma  = 1.0 / sqrt(1.0 - (vperp^2.0 + vpara^2.0)/cv^2.0);
Omega_e = 1.00

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
size_dim = size(fre, /DIMENSIONS)
size_x = size_dim[0]
size_t = size_dim[1]
for ii=0, size_x-1 do begin
  for jj=1, size_t-1 do begin
    dwdt[ii,jj] = (fre[ii,jj] - fre[ii,jj-1]) / dt
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
help, k1
s1 = 0
k2 = s2
for ii=0, size_t-1 do begin
  k2[*,ii] = 1.0 * s2[*,ii] * dbdh[2048/xskip:-2048/xskip-1]
endfor
help, k2
s2 = 0
k3 = s0 * fre * field
help, k3
s0 = 0

S = -1.0 * (k1 + k2) / k3
help, S
;S = -1.0 * (s1 * dwdt + 1.0 * s2 * dbdh[2048/xskip:-2048/xskip-1]) / (s0 * fre * field)

  ; ------------------------------------------------------------------------------
  ;                                  Plot
  ; ------------------------------------------------------------------------------
mmax = 1
cgwindow, wxsize = 960, wysize=960
!P.charsize = 2
cgImage, abs(S), minvalue=0, maxvalue=mmax, /axes, $
  xrange = [min(x_ax_f), max(x_ax_f)], yrange=[min(tf),max(tf)]/10000.0, $
  xtitle='h [c' +var.lomgc+ '!S!De0!R!U-1!N]', $
  ytitle = 't     x10!U4!N  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
  position=[0.2,0.1,0.8,0.9], /window
  
cgColorbar, /vertical, /right, range=[0,mmax],title='|S|', position=[0.81,0.1,0.82,0.9], /addcmd



end
