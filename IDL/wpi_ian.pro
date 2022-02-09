; wpi
; 

LOADCT, 39

; Global variables
GLBVAR, var

; Read files
INPUT_FILE, jobname, prefname, firnum, endnum, njob

read, 'ncpu = ', ncpu
read, 'position[1~3] = ', pstn
read, 'forward[1] backward[0]', fbwd_tmp
if (fbwd_tmp eq 1) then begin
  fbwd = 'f'
endif else begin
  fbwd = 'b'
endelse
read, 'tmax[0~1]: ', tmax_tmp
read, 'tmin[0~1]: ', tmin_tmp

position = [8192, 16384, 24576]
dt = 0.004
var.dir = '../../dat/'

jjt = 0L
FOR ijob = 0, njob-1 do begin

  ;--- Read input parameters in KEMPO from .input_idl file ---
  prm_file = var.dir +prefname +'/' +jobName(ijob) +'.prm'
  READ_KEMPOPRM, prm_file, kmp

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
  ntime = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot

  IF (ijob eq 0) then begin
    
    tmax = kmp.nstep * tmax_tmp
    tmin = kmp.nstep * tmin_tmp
    head = var.head
    print, 'tmax = ', tmax
    print, 'tmin = ', tmin
    wpi_dat = fltarr(4)
    xwind = 64 / 4
    
  endif

  wpi_ts = './wpi/wpi_ts' + string(ijob)
  wpi_xs = './wpi/wpi_xs' + string(ijob)
  wpi_vxs = './wpi/wpi_vxs' + string(ijob)
  wpi_zetas = './wpi/wpi_zetas' + string(ijob)
  wpi_ts = strcompress(wpi_ts, /REMOVE_ALL)
  wpi_xs = strcompress(wpi_xs, /REMOVE_ALL)
  wpi_vxs = strcompress(wpi_vxs, /REMOVE_ALL)
  wpi_zetas = strcompress(wpi_zetas, /REMOVE_ALL)
  Openw, 6, wpi_ts     &  Point_lun, 6, 0
  Openw, 7, wpi_xs     &  Point_lun, 7, 0
  Openw, 8, wpi_vxs    &  Point_lun, 8, 0
  Openw, 9, wpi_zetas  &  Point_lun, 9, 0
  ; ------------ wpi ------------------
  for ii = 0L, ncpu-1 do begin
    print, 'cpu: ', ii
    wpi_name = var.dir +prefname +'/' +jobName(ijob) + '.vzeta'+ string(round(pstn)) + fbwd + string(ii)
    wpi_name = strcompress(wpi_name, /REMOVE_ALL)
    Openr, 4, wpi_name  &  Point_lun, 4, 0

    while (NOT EOF(4)) do begin
      Readu, 4, head  &  Readu, 4, wpi_dat  &  Readu, 4, head
      if ((wpi_dat[0] lt tmin*dt) or (wpi_dat[0] gt tmax*dt)) then begin
	continue
      endif
      posi_tmp = abs(position[pstn-1] - wpi_dat[1])
      if (posi_tmp le xwind/2) then begin
    	writeu, 6, wpi_dat[0]
	writeu, 7, wpi_dat[1]
	writeu, 8, wpi_dat[2]
	writeu, 9, wpi_dat[3]
      endif
    endwhile
    close, 4  
  endfor
 
  close, 6  &  close, 7  &  close, 8  &  close, 9
  
ENDFOR

END
