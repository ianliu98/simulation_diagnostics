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


; ---------------------------------------------------------------------------------

var.dir = '../../dat/'
; Reading of files
jjt = 0L
FOR ijob = 0, njob-1 do begin

  ;--- Read input parameters in KEMPO from .input_idl file ---
  prm_file = var.dir +prefname +'/' +jobName(ijob) +'.prm'
  READ_KEMPOPRM, prm_file, kmp

  ; replace kmp
;  kmp = {      $
;  jobnum : 1,  $
;  nstep : 3735552, $
;  ifdiag : 2, $
;  dt : 0.004, $
;  wp1 : 15.0, $
;  wp2 : 0.3,  $
;  cv : 100.0, $
;  nx : 32768, $
;  wc : 1.0    $
;  }


  fname = var.dir +prefname +'/' + jobname(ijob) + ['.by', '.bz']

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
  ntime = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot

  IF (ijob eq 0) then begin

    total_time = ntime * njob

    ;--- Reduced size of time and space array ----------------------
    rnt = total_time / tskip
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


    t = fltarr(rnt)

    forfy = fltarr( rnx, rnt )
    forfz = forfy

    head = var.head
    ttmp = fltarr(1)
    xtmp11 = fltarr(rnx)
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

        forfy[*, jntr] = xtmp11[ jx[*] ]
        forfz[*, jntr] = xtmp22[ jx[*] ]

        jntr = ++jntr

    endif     ; mod

    jjt = ++jjt

    IF ((jjt Mod 10000) EQ 0) THEN  PRINT, 'Reading.., time : ',jjt, '/', total_time

  ENDFOR
  CLOSE, 1   &   CLOSE, 2   &   CLOSE, 3
ENDFOR

end
