; wpi
; 

LOADCT, 39

; Global variables
GLBVAR, var

; Read files
read, 'jobnum: ', njob

; video
video_s = 1
IF video_s eq 1 then begin
  vsize = [960, 960]
  fps = 10
  video = idlffvideowrite('./wpi/wpi.mp4')
  stream = video.addvideostream(vsize[0], vsize[1], fps)
  set_plot, 'z', /copy
  device, set_resolution=vsize, set_pixel_depth=24, decomposed=0
ENDIF


jjt = 0L
FOR ijob = 0, njob-1 do begin

  IF (ijob eq 0) then begin
    head = var.head
  endif

  ;wpi_ts = './wpi/wpi_ts' + string(ijob)
  ;wpi_xs = './wpi/wpi_xs' + string(ijob)
  ;wpi_vxs = './wpi/wpi_vxs' + string(ijob)
  ;wpi_zetas = './wpi/wpi_zetas' + string(ijob)
  wpi_ts = './wpi/wpi_ts' + string(ijob)
  wpi_xs = './wpi/wpi_vxs' + string(ijob)
  wpi_vxs = './wpi/wpi_vys' + string(ijob)
  wpi_zetas = './wpi/wpi_vzs' + string(ijob)
  wpi_ts = strcompress(wpi_ts, /REMOVE_ALL)
  wpi_xs = strcompress(wpi_xs, /REMOVE_ALL)
  wpi_vxs = strcompress(wpi_vxs, /REMOVE_ALL)
  wpi_zetas = strcompress(wpi_zetas, /REMOVE_ALL)
  length = FILE_INFO(wpi_ts)
  length = length.size / 4L
  wpi_t = FLTARR(length)
  wpi_x = FLTARR(length)
  wpi_vx = FLTARR(length)
  wpi_zeta = FLTARR(length)
  Openr, 6, wpi_ts   &  Point_lun, 6, 0
  Openr, 7, wpi_xs  &  Point_lun, 7, 0
  Openr, 8, wpi_vxs  &  Point_lun, 8, 0
  Openr, 9, wpi_zetas  &  Point_lun, 9, 0
  ; ------------ wpi ------------------

    while (NOT EOF(6)) do begin
      Readu, 6, wpi_t   &  Readu, 7, wpi_x
      Readu, 8, wpi_vx  &  Readu, 9, wpi_zeta
    endwhile
 
  close, 6  &  close, 7  &  close, 8  &  close, 9

  ddt = 0.004 * 16
  wpi_tt = round(wpi_t / ddt)
  t_uniq = wpi_tt[uniq(wpi_tt, sort(wpi_tt))]
  wind = 1
  ttmp = 0L

  ;for iii = 0L, N_ELEMENTS(t_uniq)-2 do begin
  for iii = 0L, 5000L do begin
    wpi_vxx = []
    wpi_zetaa = []
    for jjj = iii, iii+wind-1 do begin
      tt = t_uniq[jjj]
      if (t_uniq[jjj+1] - tt) eq 1 then begin
        wpi_vxx = [wpi_vxx, wpi_vx[where(wpi_tt eq tt)]]
        wpi_zetaa = [wpi_zetaa, wpi_zeta[where(wpi_tt eq tt)]]
      endif
    endfor
    
    ; calculation
    wpi_vxx = wpi_vxx / 100
    
    ; plot
    ymin = -0.5
    ymax = 0
    cgscatter2D, wpi_zetaa[where(wpi_vxx ge ymin)], wpi_vxx[where(wpi_vxx ge ymin)], $
	    yrange=[ymin, ymax], xrange=[0, !PI*2], position=[0.1, 0.1, 0.9, 0.9], $
	    fit=0, PSym=3, ytitle='v_perp', xtitle='zeta'

    if (video_s eq 1) then begin
      timestamp = video.put(stream, tvrd(true=1))
      print, 'count: ', ttmp
      ttmp = ttmp + 1
    endif else begin    
      dummy = ''
      read, 'press Enter: ', dummy
    endelse
    
  endfor
  
ENDFOR

if (video_s eq 1) then begin
  device, /close
  video.cleanup
endif

END
