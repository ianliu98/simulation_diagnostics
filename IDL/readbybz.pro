;----------------------------------------------------------------------------------
;   calculate inhomogeneity factor S
;----------------------------------------------------------------------------------

; Shrink an array
tskip = 1
xskip = 1
LOADCT, 39

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

  fname = var.dir +prefname +'/' + jobname(ijob) + ['.by', '.bz']

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
  ntime = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot

  IF (ijob eq 0) then begin

    total_time = ntime * njob
    rnt = total_time/ tskip

    rnx = kmp.nx / xskip
    xmin_f = kmp.nxl/xskip
    xmax_f = (kmp.nx-kmp.nxr)/xskip

    forfy = Fltarr( xmax_f-xmin_f, rnt)
    forfz = forfy
    tmpy = forfy
    tmpz = forfz

    t = Fltarr(rnt)
    head = var.head
    ttmp = Fltarr(1)
    xtmp11 = Fltarr(kmp.nx)
    xtmp22 = xtmp11[*]

    jx = Lindgen(xmax_f-xmin_f, increment=xskip, start=xmin_f)

    jnt = 0L
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

xlen = xmax_f - xmin_f

for ix=0, xlen-1 do begin

  tmpy = FFT( forfy[ix,*], -1)
  tmpz = FFT( forfz[ix,*], -1)

endfor

end
