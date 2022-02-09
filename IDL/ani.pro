;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;		Anisotropy & miss count
;			ian 2021.1
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

LOADCT, 39

GLBVAR, var
INPUT_FILE, jobname, prefname, firnum, endnum, njob

once  = 0
tskip = 1
diag  = 256

FOR ijob = 0, njob-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  in job loop  ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	datjob   = '../../dat/' + prefname + '/' + jobname[ijob]
	file     = datjob + '.ani'
	prm_file = datjob + '.prm'
	READ_KEMPOPRM, prm_file, kmp 
	nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
	ntstep = (kmp.nstep - nprev) /diag

	if (once eq 0) then begin
	;+++++++++++++++++++++++++++++++++++++++
	;++++++++++ preset for arrays ++++++++++
	;+++++++++++++++++++++++++++++++++++++++
		total_tstep = ntstep * njob/ tskip
		time        = Dblarr(total_tstep)
		anisotropy  = time
		miss_count  = lonarr(total_tstep)

		time_tmp    = Dblarr(1)
		anis_tmp    = time_tmp
		miss_tmp    = lonarr(1) 
		jjt         = 0L
		
		once = 1
	endif
	
	print, 'Opening ', file   &   Openr, 1, file   &   Point_lun, 1, 0
	head  = var.head

	for jt=0L, (ntstep/tskip)-1 do begin
		readU, 1, head  &  readU, 1, time_tmp, anis_tmp, miss_tmp   &  readU, 1, head
		time[jjt]        =  time_tmp
		anisotropy[jjt]  =  anis_tmp
		miss_count[jjt]  =  miss_tmp
		jjt = jjt + 1
	endfor

	close,1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  end job loop  ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ENDFOR

cgWindow, woxmargin = [5, 5], woymargin = [5, 5], WXSize=1080, WYSize=960
cgplot, time, anisotropy, yrange=[min(anisotropy),max(anisotropy)], /axes, $
	wmulti=[0,1,2], ytitle='anisotropy', xtitle='time', title = file,  $
	/window
cgplot, time, miss_count, /axes, ytitle='miss count', xtitle='time',       $
	title = file, /addcmd

end
