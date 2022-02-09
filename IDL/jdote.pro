;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    J dot E
;      option: bandpass
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


LOADCT, 39

GLBVAR, var
READ, 'Which component?, J*E[1], J*B[2]            : ', eorb
READ, 'Cold or Hot or All?, Cold[1], Hot[2], All[3]: ', corh
READ, 'how many cases in comparison: ', cases

components = [['.j1e', '.j1b'], $
			  ['.j2e', '.j2b'], $
			  ['.je' , '.jb' ]]

bandpass = 0
log_plt  = 0

tskip = 4
xskip = 4
diag  = 256
dt    = 0.004

FOR case_ind = 0, cases-1 do begin

INPUT_FILE, jobname, prefname, firnum, endnum, njob

FOR ijob = 0, njob-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  in job loop  ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	datjob   = '../../dat/' + prefname + '/' + jobname[ijob]
	file     = datjob + components[eorb-1,corh-1]
	file2    = datjob + ['.by', '.bz']
	prm_file = datjob + '.prm'
	READ_KEMPOPRM, prm_file, kmp 
	nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
	ntstep = (kmp.nstep - nprev) /diag / tskip

	if (ijob eq 0) then begin
	;+++++++++++++++++++++++++++++++++++++++
	;++++++++++ preset for arrays ++++++++++
	;+++++++++++++++++++++++++++++++++++++++
		total_tstep = ntstep * njob
		time        = Fltarr(total_tstep)
		start_time  = (nprev + diag) * kmp.dt
		jeb         = Fltarr(kmp.nx/xskip,total_tstep)
		je_fft      = jeb[*,*]
		by	    = jeb[*,*]
		bz	    = jeb[*,*]
		time_tmp    = Fltarr(1)
		jeb_tmp     = Fltarr(kmp.nx)
		by_tmp      = Fltarr(kmp.nx)
		bz_tmp      = Fltarr(kmp.nx)
                jx          = Lindgen(kmp.nx/xskip, increment=xskip, start=0)
		jjt         = 0L
	endif
	
	print, 'Opening ', file      &   Openr, 1, file      &   Point_lun, 1, 0
	print, 'Opening ', file2[0]  &   Openr, 2, file2[0]  &   Point_lun, 2, 0
	print, 'Opening ', file2[1]  &   Openr, 3, file2[1]  &   Point_lun, 3, 0
	head  = var.head

	for jt=0L, ntstep-1 do begin
		readU, 1, head  &  readU, 1, time_tmp   &  readU, 1, head
		readU, 1, head  &  readU, 1, jeb_tmp    &  readU, 1, head
		readU, 2, head  &  readU, 2, time_tmp   &  readU, 3, head
		readU, 2, head  &  readU, 2, by_tmp     &  readU, 3, head
		readU, 3, head  &  readU, 3, time_tmp   &  readU, 3, head
		readU, 3, head  &  readU, 3, bz_tmp     &  readU, 3, head
		time[jjt]  = time_tmp
		jeb[*,jjt] = jeb_tmp[ jx[*] ]
		by[*,jjt]  = by_tmp[ jx[*] ]
		bz[*,jjt]  = bz_tmp[ jx[*] ]
		jjt = jjt + 1
	endfor

	close,1
	close,2
	close,3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;  end job loop  ;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ENDFOR

;jeb = jeb[255:4096-256,*]

;if (eorb eq 2) then begin
;  jeb = jeb * (dt / 2)
;endif

fld = sqrt(by^2 + bz^2)
jeb = jeb / fld


if (bandpass eq 1) then begin

; bandpass filter
dts = dt * diag * tskip
ws  = 2.0 * !dpi / dts
fft_pnt_f = jjt
dw_f = ws / fft_pnt_f
omg_ar_f = Findgen(fft_pnt_f) * dw_f - dw_f / 2.0

w_l = 0.04
w_h = 0.06

jwmin = Max( where(omg_ar_f lt w_l) )
jwmax = Min( where(omg_ar_f ge w_h) )

for ix=0, N_ELEMENTS(jx)-1 do begin
  tmpj = FFT( jeb[ix,*], -1)
  tmpj[0:jwmin] = 0  &  tmpj[jwmax:*] = 0
  je_fft[ix, *] = 2*real_part( FFT( tmpj, 1) )
endfor

endif

;---------------------------------------
;------------    plot   ----------------
;---------------------------------------

if (case_ind eq 0) then begin

axis_format = {xrange: [-(kmp.nx-512)/200, (kmp.nx-512)/200], $
	       yrange: [0, total_tstep*diag*tskip*dt]/10000.0, $
	       font  : 1, $
	       charsize: 3, $
	       xtitle: 'h [c' +var.lomgc+ '!S!De0!R!U-1!N]', $
	       ytitle: 't     x10!U4!N[' +var.lomgc+ '!S!De0!R!U-1!N]'}

mmargin = [3, 6, 3, 6]

cgWindow, woxmargin=[10,15], woymargin=[5,5], WXSize=1080, WYSize=640

endif

if (bandpass eq 1) then begin

  if (log_plt eq 1) then begin

    ; log value
    je_abs     = abs(je_fft)
    je_abs_log = alog10(je_abs)

    if (case_ind eq 0) then begin
      pmax = max(je_abs_log)
      pmin = pmax - 1.5
      cgimage, je_abs_log, minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, wmulti=[0,cases,1],multimargin=mmargin, /window
    endif else begin
      axis_format.ytitle = ''
      cgimage, je_abs_log, minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, multimargin=mmargin, /addcmd
      xwindow = !X.window  &  ywindow = !Y.window
    endelse 

  endif else begin
    
    if (case_ind eq 0) then begin
      pmax = max(abs(je_fft)) / 2.0
      pmin = 0
      cgimage, abs(je_fft), minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, wmulti=[0,cases,1],multimargin=mmargin, /window
    endif else begin
      axis_format.ytitle = ''
      cgimage, abs(je_fft), minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, multimargin=mmargin, /addcmd
      xwindow = !X.window  &  ywindow = !Y.window
    endelse 

  endelse

endif else begin


  if (log_plt eq 1) then begin

    ; log value
    je_abs     = abs(jeb)
    je_abs_log = alog10(je_abs)

    if (case_ind eq 0) then begin
      pmax = max(je_abs_log) -0.5
      pmin = pmax - 3
      cgimage, je_abs_log, minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, wmulti=[0,cases,1],multimargin=mmargin, /window
    endif else begin
      axis_format.ytitle = ''
      cgimage, je_abs_log, minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, multimargin=mmargin, /addcmd
      xwindow = !X.window  &  ywindow = !Y.window
    endelse 

  endif else begin
    
    if (case_ind eq 0) then begin
      pmax = max(abs(jeb)) / 2.0
      pmin = 0 
      pmax = 0.1 
      pmin = 0
      cgimage, abs(jeb), minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, wmulti=[0,cases,1],multimargin=mmargin, /window
    endif else begin
      axis_format.ytitle = ''
      cgimage, abs(jeb), minvalue=pmin, maxvalue=pmax, AXKEYWORDS=axis_format, /axes, $
	       /interpolate, multimargin=mmargin, /addcmd
      xwindow = !X.window  &  ywindow = !Y.window
    endelse 

  endelse

endelse

ENDFOR

if (eorb eq 2) then begin
  bartitle = 'log!D10!N|J!DB!N / B!Dw!N|'
endif else begin
  bartitle = 'log!D10!N|J!DE!N|'
endelse
cgColorbar, /vertical, /right, yminor=1, range=[pmin,pmax],  charsize=2.5, title=bartitle,$
	    position=[xwindow[1]-0.01, ywindow[0]+0.03, xwindow[1], ywindow[1]-0.03], /addcmd

end
