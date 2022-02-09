LOADCT, 33
GLBVAR, var

; Read files
fname = ''      ; define fname as string
READ, 'main name of file ? : ', fname
READ, 'first number of files to input ? : ', firnum
READ, 'end number of files to input ? : ', endnum
njob = endnum - firnum + 1
fnum = STRING(FIX(firnum) + INDGEN(njob))
jobname = fname + fnum
jobname = STRCOMPRESS(jobname, /REMOVE_ALL)     ; delete some spaces
PRINT, ''
PRINT, 'open files ... ', jobname


fnames = ['.ex', '.ey', '.ez', '.bx', '.by', '.bz']
fnum1 = 25
fnum2 = 26

; data skip (dt = ifdiag * dskp)
;PRINT, 'default : dskp = 1'
;READ, 'time skip (1 or even) ? : ', dskp
dskp = 4

; Read files
fchk = 0
jjt = 0L
jjjt = 0L
head = BYTARR(8)     ; --- header of file ---
FOR ijob = 0, njob-1 DO BEGIN
	prmfile = '../../dat/' + fname +'/'+ jobname(ijob)+ '.prm'          ; parameter file
	fFile1 = '../../dat/' +  fname +'/'+ jobname(ijob) + fnames(fnum1 - 21)     ; data file
	fFile2 = '../../dat/' +  fname +'/'+ jobname(ijob) + fnames(fnum2 - 21)

	READ_KEMPOPRM, prmfile, kmp
	nx = kmp.nx
	dr = kmp.dr
	dt = kmp.dt
	jobnum = kmp.jobnum
	nstep = kmp.nstep
	ifdiag = kmp.ifdiag
	cv = kmp.cv

	nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
	nt = (kmp.nstep - nprev) /kmp.ifdiag     ; total number of time to plot
	IF (fchk EQ 0) THEN BEGIN
		xcen = nx / 2
		xx = xcen	; decide region of fft
		xmin = 0
		xmax = nx - 1
		ttime = nt * njob
		print,ttime
		ret = ttime/ dskp
		start_time = (nprev + kmp.ifdiag) * kmp.dt
		field1 = complexarr(xx, ret)
		field2 = field1
		fld = Fltarr(nx)
		fld2 = fld[*]
		t = FLTARR(ret)
		tmp = FLTARR(nx)
		time = FLTARR(1)
		fchk = 1
		dk = 2.*!PI/(2*xx)/(dr/cv)
		nn = 4
		testf = sin(nn*dk*findgen(xx))
		testfft = fft(testf,-1)
	ENDIF

	PRINT, 'opening ', fFile1     ; open the file
	PRINT, 'opening ', fFile2
	OPENR, 1, fFile1
	POINT_LUN, 1, 0               ; beginning of the file
	OPENR, 2, fFile2
	POINT_LUN, 2, 0

	head = var.head
	FOR jt = 0, nt-1 DO BEGIN
		READU, 1, head  &  READU, 1, time  &  READU, 1, head
		READU, 1, head  &  READU, 1, fld   &  READU, 1, head      
		READU, 2, head  &  READU, 2, time  &  READU, 2, head
		READU, 2, head  &  READU, 2, fld2  &  READU, 2, head
		IF ((jjt MOD dskp) EQ 0) THEN BEGIN
			t(jjjt) = time
			divf = DIVFIELD(fld[*],fld2[*])
			fy_f = divf[*,0]
			fz_f = divf[*,1]
			fy_b = divf[*,2]
			fz_b = divf[*,3]
			fft_field1 = fft(fy_f[xmin:xmax-1],-1)
			fft_field2 = fft(fz_f[xmin:xmax-1],-1)
			field1(*,jjjt) = 2. * fft_field1[0:xx-1]
			field2(*,jjjt) = 2. * fft_field2[0:xx-1]
			jjjt = jjjt + 1
			ENDIF
			jjt = jjt + 1
	ENDFOR
	CLOSE, 1  &  CLOSE, 2  &  CLOSE, 3
ENDFOR     ; end of job loop
 
field1 = abs(field1)
field2 = abs(field2)

field = 2. * SQRT( field1^2 + field2^2 )
print,max(field)

x = findgen(nx)

k = FINDGEN(2*xx)
dk = 2.*!PI/ (2*xx)/ (dr/cv)
k = k(*) * dk
ku = k[0:xx-1]

; Parameters
plot_max = 2e-5
plot_min = plot_max * 0.1
dlevel = plot_max - plot_min
PRINT, 'plot_max, plot_min ', plot_max, plot_min

; Plot
cgwindow, wxsize=1024, wysize=768
ppos=[.1, .1, .9, .9]

kk = 1000
tt = 2000

xyaxes = { xrange: [ku[0], ku[kk]], yrange: [t[0], t[tt]]/10000.0, $
           xtitle: 'k [c!E-1!N!4X!X!De0!N]', ytitle:'t x10!U4!N[!4X!X!De0!N!E-1!N]'}

cgimage, field[0:kk,0:tt], ku[0:kk], t[0:tt],  $
	maxvalue=plot_max, minvalue=plot_min,  $
	/axes, position=ppos, /interpolate,    $
	axkeywords=xyaxes,/window

cgcolorbar, /vertical, /right, range=[plot_min,plot_max], $
	position=[0.91,0.1,0.92,0.9], /addcmd
	
;SHADE_SURF, field(0:kk, *), ku[0:kk], t,AX=90,AZ=0, $
;   MAX_VALUE=plot_max, MIN_VALUE=plot_min,  $
;   XST=1, YST=1, ZST=4, TICKLEN=-.02,  $
;   XTITLE='k [c!E-1!N!4X!X!De0!N]', YTITLE='t [!4X!X!De0!N!E-1!N]',  $
;   SHADES=BYTSCL(field(0:kk, *),  $
;   MAX=plot_max, MIN=plot_min, TOP=!D.TABLE_SIZE),  $
;   POS=ppos
;COLORBAR, /VERTICAL, /YLOG,  $
;   RANGE=[plot_min, plot_max], $
;   MINOR=9, TICKLEN=-.2, YTICKS=0, YTICKFORMAT='exponent',  $
;   POS=[ppos(2)+.08, ppos(1), ppos(2)+.11, ppos(3)], /addcmd

; Label
;XYOUTS, .01, .97, SYSTIME(0), SIZE=1., /NORMAL, /addcmd
;XYOUTS, .3, .97, 'filename :' + STRING(fname,'(a5)') $
;   + STRING(firnum,'(I3)') + ' - ' + STRING(endnum,'(I3)'), SIZE=1., /NORMAL, /addcmd

end
