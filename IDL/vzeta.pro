; wpi
;

LOADCT, 39

GLBVAR, var

INPUT_FILE, jobname, prefname, firnum, endnum, njob

; video
video_s = 0
IF video_s eq 1 then begin
  vsize = [960, 960]
  fps = 10
  video = idlffvideowrite('./wpi/vzeta.mp4')
  stream = video.addvideostream(vsize[0], vsize[1], fps)
  set_plot, 'z', /copy
  device, set_resolution=vsize, set_pixel_depth=24, decomposed=0
ENDIF

var.dir = '../../dat/'
  
vdiv = 200
zdiv = 100
ivzeta = 8
ave = 1

FOR ijob = 0, njob-1 do begin
  
  prm_file = var.dir + prefname + '/' + jobname(ijob) + '.prm'
  READ_KEMPOPRM, prm_file, kmp
  
  file = var.dir + prefname + '/' + jobname(ijob) + '.vzeta'

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
  ntime = (kmp.nstep - nprev) / ivzeta
  ntime = ntime / 2
  
  IF (ijob eq 0) then begin
    total_time = ntime * njob
    head = var.head
    vzeta = DBLARR(zdiv+1, vdiv+1)
    test = vzeta
    vzeta_ave = DBLARR(zdiv+1, vdiv+1, ave)
    vzeta_ini = vzeta[*,*]
    potential = vzeta[*,*]
    zeta_tmp = DBLARR(zdiv+1)
    tmp = 0

    ttmp = vzeta[*,*]

  ENDIF
  
  Openr, 2, file  &  Point_lun, 2, 0 

  for jt=0, ntime-1 do begin
    
    for jv=0, vdiv do begin
      readu, 2, head  &  readu, 2, zeta_tmp  &  readu, 2, head
      vzeta[*, jv] = zeta_tmp[*]
    endfor

    ttmp = vzeta

    test = vzeta
    vzeta_ave[*, *, jt mod ave] = vzeta

    if jt eq ave-1 then begin
    ;if jt eq 12500 then begin
      vzeta_ini = 0.0
      for ii = 0, ave-1 do begin
	vzeta_ini = vzeta_ini + vzeta_ave[*, *, ii]
      endfor
      vzeta_ini = vzeta_ini / ave
      maxv = 0.1
      minv = -0.1
    endif

    if jt ge ave-1 then begin
      vzeta = 0.0
      for jj = 0, ave-1 do begin
	vzeta = vzeta + vzeta_ave[*, *, jj]
      endfor
      vzeta = vzeta / ave
      vzeta = (vzeta - vzeta_ini) / vzeta_ini
      ;vzeta = alog10(vzeta/vzeta_ini)
      vzeta[where(vzeta gt maxv+0.05)] = maxv
      ;vzeta[where(vzeta lt 0)] = 0
      ; plot
      cgImage, vzeta[0:zdiv, 0:vdiv/2], /axes, interpolate=1, minvalue=minv-0.05, maxvalue=maxv+0.05, $
	  font=1, charsize=1, position=[0.1, 0.1, 0.9, 0.9]
      cgColorbar, /vertical, /right, range=[minv-0.05, maxv+0.05], position=[0.91, 0.1, 0.92, 0.9]
    
      if (video_s eq 1) then begin
        timestamp = video.put(stream, tvrd(true=1))
      endif else begin
        dummy = ''
        read, 'press Enter: ', dummy
    endelse

    endif

  endfor
    close, 2

ENDFOR

if (video_s eq 1) then begin
  device, /close
  video.cleanup
endif

END
