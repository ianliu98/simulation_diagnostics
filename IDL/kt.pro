LOADCT, 39
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

READ, 'make psfile ? y[1]/ n[2] ', ps
;ps = 2
IF (ps EQ 1) THEN BEGIN
   SET_PLOT, 'ps'
   psname = './IMG/ktcont.eps'
   DEVICE, filename = psname, /COLOR, BITS=8, /LAND, /ENCAPSULATED
ENDIF

fnames = ['.ex', '.ey', '.ez', '.bx', '.by', '.bz']
;read, 'field number? (21 - 26) ', field_num
fnum1 = 25
fnum2 = 26

;xmin = 3000
;xmax = 14000
;xran = xmax - xmin + 1
;print, 'xmin, xmax : ', xmin, xmax

; data skip (dt = ifdiag * dskp)
PRINT, 'default : dskp = 1'
READ, 'time skip (1 or even) ? : ', dskp

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
      ;IF (jobnum EQ 2) THEN BEGIN
      ;   ret = ret
      ;   nt = nt
      ;ENDIF
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
   ;WHILE ~EOF(1)  DO BEGIN
      ;print,~EOF(1)
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
	 ;print,time
         fft_field1 = fft(fy_f[xmin:xmax-1],-1)
         fft_field2 = fft(fz_f[xmin:xmax-1],-1)
         field1(*,jjjt) = 2. * fft_field1[0:xx-1]
         field2(*,jjjt) = 2. * fft_field2[0:xx-1]
         jjjt = jjjt + 1
      ENDIF
      jjt = jjt + 1
      ;print,~EOF(1)
   ENDFOR
   ;ENDWHILE
   CLOSE, 1  &  CLOSE, 2  &  CLOSE, 3
ENDFOR     ; end of job loop
 
 ; for ii = 0,jjjt-1 DO BEGIN
;	ttmp1 = where(field1(*,ii) GT 1)
;	ttmp2 = where(field2(*,ii) GT 1)
;	field1(ttmp1,ii) = 0
;	field2(ttmp2,ii) = 0
 ; endfor
 
  ;ffd1 = field1
  ;ffd2 = field2
;  field1 = FFT( field1(*,*), -1, dimension = 1)
;  field2 = FFT( field2(*,*), -1, dimension = 1)
 
;  help,field1
 
;pe = 200
;field1(0:pe,*) = 0.
;field1(nx-pe:*,*) = 0.
;field2(0:pe,*) = 0.
;field2(nx-pe:*,*) = 0.

;   ; IFFT of k-t
;   field1 = FFT( field1(*,*), 1)
;   field2 = FFT( field2(*,*), 1)
   
   ;field1 = REAL_PART(field1)
   ;field2 = REAL_PART(field2)
   
   ;field1 = IMAGINARY(field1)
   ;field2 = IMAGINARY(field2)
   
   field1 = abs(field1)
   field2 = abs(field2)

   ;for ii = 0, ret DO BEGIN
;	field_1(ii,*) = field1
;	field_2(ii,*) = field2
 ;  ENDFOR	

   field = 2. * SQRT( field1^2 + field2^2 )
   print,max(field)
   ;field = ALOG10(field(*,*) > 1.e-20)



;     plot_max = MAX(field(*,*))
;     plot_max = -1.3

;     dlevel = 3
;     plot_min = FIX(plot_max - dlevel)

;     IF (ps EQ 2) THEN  WINDOW, XSIZE=1024, YSIZE=768
;     ppos=[.1, .15, .85, .9]

x = findgen(nx)
;xp = 300
;xp = nx-1
;     SHADE_SURF, field(0:xp, *), x(0:xp), t(*), /T3D,  $
;     SHADE_SURF, field(nx-xp:*, *), x(nx-xp:*), t(*), /T3D,  $
;        XST=1, YST=1, ZST=4, TICKLEN=-.02,  $
;        XTITLE='x', YTITLE='t [!4X!X!De0!N!E-1!N]',  $
;        SHADES=BYTSCL(field(0:xp, *),  $
;        SHADES=BYTSCL(field(nx-xp:*, *),  $
;      MAX=plot_max, MIN=plot_min, TOP=!D.TABLE_SIZE),  $
;       POS=ppos

;     COLORBAR, /VERTICAL, /YLOG,  $
;        RANGE=[10.^plot_min, 10.^plot_max], DIVISIONS=dlevel/10,  $
;        MINOR=9, TICKLEN=-.2, YTICKS=0, YTICKFORMAT='exponent',  $
;        POS=[ppos(2)+.08, ppos(1), ppos(2)+.11, ppos(3)]


 ;  field = SHIFT( field, nx/2-1, 0 )


; k array
     k = FINDGEN(2*xx)
  ;   k(nx/2+1:*) = k(nx/2+1:*) - nx     ; insert negative k from nx/2 + 1
   ;  k = SHIFT(k, nx/2-1)            ; shift 0 to center of array
     dk = 2.*!PI/ (2*xx)/ (dr/cv)
     k = k(*) * dk
     ku = k[0:xx-1]

; Parameters
     ;plot_max = MAX(field(*,*))
     plot_max = 2e-5
     plot_min = plot_max * 0.1
     ;plot_min = 1e-6
     dlevel = plot_max - plot_min
     ;plot_min = FIX(plot_max - dlevel)
     PRINT, 'plot_max, plot_min ', plot_max, plot_min

; Plot
     WINDOW, XSIZE=1024, YSIZE=768
     ppos=[.1, .15, .85, .9]

kk = 1000

     SHADE_SURF, field(0:kk, *), ku[0:kk], t,AX=90,AZ=0, $
        MAX_VALUE=plot_max, MIN_VALUE=plot_min,  $
;     SHADE_SURF, fiel0(0:kk, *), k(0:kk), t(*), /T3D,  $
        XST=1, YST=1, ZST=4, TICKLEN=-.02,  $
        XTITLE='k [c!E-1!N!4X!X!De0!N]', YTITLE='t [!4X!X!De0!N!E-1!N]',  $
        SHADES=BYTSCL(field(0:kk, *),  $
      MAX=plot_max, MIN=plot_min, TOP=!D.TABLE_SIZE),  $
       POS=ppos

     COLORBAR, /VERTICAL, /YLOG,  $
        RANGE=[plot_min, plot_max], $
        MINOR=9, TICKLEN=-.2, YTICKS=0, YTICKFORMAT='exponent',  $
        POS=[ppos(2)+.08, ppos(1), ppos(2)+.11, ppos(3)]

; Label
     XYOUTS, .01, .97, SYSTIME(0), SIZE=1., /NORMAL
     XYOUTS, .3, .97, 'filename :' + STRING(fname,'(a5)') $
        + STRING(firnum,'(I3)') + ' - ' + STRING(endnum,'(I3)'), SIZE=1., /NORMAL

IF (ps EQ 1) THEN BEGIN
   PRINT, 'output to', psname
   DEVICE, /CLOSE
ENDIF


end
