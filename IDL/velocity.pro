; initial
LOADCT, 39

tskip = 1

GLBVAR, var
INPUT_FILE, jobname, prefname, firnum, endnum, njob
species = ['c', 'h']
component = ['.vx', '.vperp', '.vy', '.vz']
Read, 'cold[1], hot[2]: ', species_in
Read, 'vpara[1], vperp[2], v_y[3], v_z[4]: ', component_in

if (component_in eq 1) then begin
  v_range = 10001
  xtick = findgen(v_range)/50 - 100
endif else begin
  v_range = 5001
  xtick = findgen(v_range)/50
endelse

suffix = component[component_in-1] + species[species_in-1]


;video
vsize = [1920, 1080]
fps = 10
video=idlffvideowrite('./save_file/'+jobname[0]+'.mp4')
stream=video.addvideostream(vsize[0], vsize[1], fps)
set_plot, 'z', /copy
device, set_resolution=vsize, set_pixel_depth=24, decomposed=0
frameskip = 512



jjt = 0L
FOR ijob = 0, njob-1 do begin
  
  datjob = '../../dat/' + prefname + '/' + jobname[ijob]
  file = datjob + suffix
  prm_file = datjob +'.prm'
  READ_KEMPOPRM, prm_file, kmp
  kmp.ifdiag = 256
  
  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum 
  ntime = (kmp.nstep - nprev) /kmp.ifdiag 

  IF (ijob eq 0) then begin
    total_time = ntime * njob
	rnt = total_time/ tskip
	start_time = (nprev + kmp.ifdiag) * kmp.dt
	
	time = Dblarr(rnt)
	v = lonarr(rnt, v_range)
	time_tmp = Dblarr(1)
	v_tmp = lonarr(v_range)

	head = var.head
  ENDIF
  
  OPENR, 1, file  &   Point_lun, 1, 0
  print, 'Opening : ', file
  
  for jt=0L, ntime-1 do begin
    readU, 1, head  &  readU, 1, time_tmp, v_tmp &  readU, 1, head
    pos = [0.1, 0.2, 0.83, 0.85]
    if ((species_in eq 1) && (component_in eq 1)) then begin
      cgplot, xtick[4750:5250], v_tmp[4750:5250], /axes, ytitle='number', xtitle='velocity', title = 'time='+strtrim(time_tmp)
    endif else begin
      cgplot, xtick, v_tmp, /axes, ytitle='number', xtitle='velocity', title = 'time='+strtrim(time_tmp)
    endelse 
    if (jjt mod frameskip) eq 0 then begin
      timestamp = video.put(stream, tvrd(true=1))
    endif
    time[jjt] = time_tmp
    v[jjt,*] = v_tmp
    jjt = jjt + 1
  endfor
  
  close,1
  
ENDFOR

device, /close
video.cleanup

end
