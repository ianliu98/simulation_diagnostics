;----------------------------------------------------------------------------------
;   calculate inhomogeneity factor S
;----------------------------------------------------------------------------------

; Shrink an array
tskip = 1
xskip = 1

; Global variables
GLBVAR, var

; Read files
INPUT_FILE, jobname, prefname, firnum, endnum, njob

; fft position
wpi_posi = 3
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
    rnt = total_time/ tskip
    rnx = kmp.nx / xskip

    ; info
    fsample =  2.d0*!dpi/ (kmp.ifdiag * kmp.dt)
    fsample_r = fsample/ tskip
    omg = 0.3  &  b0 = 1
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


    t = Fltarr(rnt)

    xmin_f = kmp.nxl/xskip
    xmax_f = (kmp.nx-kmp.nxr)/xskip

    forfy = Fltarr( rnx, rnt )
    forfz = forfy

    head = var.head
    ttmp = Fltarr(1)
    xtmp11 = Fltarr(kmp.nx)
    xtmp22 = xtmp11[*]

    jx = Lindgen(rnx, increment=xskip, start=0)

    jnt = 0L
    jkt = 0L
    jntl = 0L
    jntr = 0L
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

        t[jntr] = ttmp

        divf = DIVFIELD(xtmp11, xtmp22)
        divf = Float( divf[*,*])

        tmp = divf[*,0]
        forfy[*, jntr] = tmp[ jx[*] ]
        tmp = divf[*,1]
        forfz[*, jntr] = tmp[ jx[*] ]

        jntr = ++jntr

    endif     ; mod

    jjt = ++jjt

    IF ((jjt Mod 10000) EQ 0) THEN  PRINT, 'Reading.., time : ',jjt, '/', total_time

  ENDFOR
  CLOSE, 1   &   CLOSE, 2   &   CLOSE, 3
ENDFOR

; release the memory
xtmp11 = 0
xtmp22 = 0
tmp = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FFT

dts = kmp.dt * kmp.ifdiag * tskip   ; Sampling time [Omega_e^-1]
ws = 2.d0 * !dpi/ dts   ; Sampling freq. (angular) [Omega_e]

fft_pnt_f = jntr
dw_f = ws/ fft_pnt_f ; Frequency resolution (angular) [Omega_e^-1]
omg_ar_f  = Findgen(fft_pnt_f) * dw_f - dw_f/2.0

minf = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]
maxf = [0.02, 0.04, 0.06, 0.08, 0.10, 0.30]

jwmin = Max( Where(omg_ar_f lt minf[2]) )
jwmax = Min( Where(omg_ar_f ge maxf[2]) )

print, 'jfmin, jfmax: ', jwmin, jwmax

wpi_wind = 256
wpi_poss = [4096, 8192, 12288, 16384, 20480, 24576, 28672]

left = wpi_poss[wpi_posi]/xskip - wpi_wind / 2
righ = wpi_poss[wpi_posi]/xskip + wpi_wind / 2  ; xskip should be 1

fy_flt = Fltarr(righ-left, rnt)
fz_flt = fy_flt

print, 'FFT...'

for ix=left, righ-1 do begin

  tmpy = FFT( forfy[ix,*], -1)
  tmpz = FFT( forfz[ix,*], -1)

  tmpy[0:jwmin] = 0  &  tmpy[jwmax:*] = 0
  tmpz[0:jwmin] = 0  &  tmpz[jwmax:*] = 0

  fy_flt[ix-left, *] = 2*real_part( FFT( tmpy, 1) )
  fz_flt[ix-left, *] = 2*real_part( FFT( tmpz, 1) )

endfor

data_save_y = './fy_flt_3.sav'
data_save_z = './fz_flt_3.sav'

openw, 1, data_save_y  &  point_lun, 1, 0
openw, 2, data_save_z  &  point_lun, 2, 0

for it=0, rnt-1 do begin

  writeu, 1, head  &  writeu, 1, fy_flt[*, it]  &  writeu, 1, head
  writeu, 2, head  &  writeu, 2, fz_flt[*, it]  &  writeu, 2, head

endfor

close, 1  &  close, 2

end
