;+
; PURPOSE:
;      field plot (screen shot or continuous mpeg)        
; 		1 or 2 field component plot
; ------------------------------------------------------------
; read and display 'field data' obtained by the kempo code      
; Apr 12, 2007  by  M.Hikishima	
;
;-

LOADCT, 39

; Get global variables
 GLBVAR, var

 ; Components to plot
 ; Read, 'Field E[1],B[2] ? : ', wfield
 wfield = 2
 if (wfield eq 1) then  fcomp = ['.ey', '.ez']
 if (wfield eq 2) then  fcomp = ['.by', '.bz']
 if (wfield eq 3) then  fcomp = ['.je_f', '.je_b']
 if (wfield eq 4) then  fcomp = ['.jb_f', '.jb_b']

 ;;;;;;;;;;;;;;;;;;;;;;;;;
 tskp = 1     ; time skip
 ;     tskp = 512     ; time skip
 ;     tskp = 1     ; time skip
 ;;;;;;;;;;;;;;;;;;;;;;;;;

colors = ['black','red','blue','purple','green','yellow']
Read, 'Choose a position [0 ~ 32768]: ', plt_pos
Read, 'Range to plot [default 2000]: ', plt_rng
Read, 'Comparison for different cases ? Yes[0], No[1]: ', cmpr

if (cmpr eq 0) then begin
  names = ''
  Read, 'How many cases are involved? :', cases
  for i = 0,cases-1 do begin
    INPUT_FILE, jobname, prefname, firnum, endnum, njob
    
    names = [names,prefname]
    
    ; Reading files
    jfrm = 0
    jjt = 0L
    jnt = 0L
    
    time_final = [];
    field_final = [];
    for ijob = 0, njob-1 do begin    ; job loop

      fcomp1 = '../../dat/' +prefname +'/' +jobName[ijob] + fcomp[0]
      fcomp2 = '../../dat/' +prefname +'/' +jobName[ijob] + fcomp[1]

      ; Read input parameters in KEMPO from .input_idl file
      prm_file = '../../dat/' +prefname +'/'  + jobName(ijob) + '.prm'
      READ_KEMPOPRM, prm_file, kmp

      nprev = (kmp.jobnum - 1) * kmp.nstep/ kmp.jobnum
      ntime_total = (kmp.nstep - nprev)/ kmp.ifdiag
      ntime = ntime_total / tskp
      field_av = Fltarr(ntime)
      time_ar = Fltarr(ntime)

      if ( (wfield eq 3) or (wfield eq 4) ) then  $
        ntime = (kmp.nstep - nprev)/ kmp.mltstp

      if (ijob eq 0) then begin
        field = Fltarr(kmp.nx)
        field2 = field[*]
        ttmp  = Fltarr(1)
        ftmp1 = Fltarr(kmp.nx)
        ftmp2 = ftmp1[*]
        x_ax = ( Findgen(kmp.nx) - (kmp.nx-1)/2. ) * kmp.dr/ kmp.cv
        xcent = kmp.nx/ 2

        ws = 2.*!dpi/ (kmp.dt*kmp.ifdiag*tskp)
        print, 'samplinf freq: ', ws

      endif

      ;--- Open files ---
      print, 'Opening ', fcomp1   &   Openr, 1, fcomp1   &   Point_lun, 1, 0
      print, 'Opening ', fcomp2   &   Openr, 2, fcomp2   &   Point_lun, 2, 0


      ;----- Time loop in each job -----

      head = var.head

      x1 = plt_pos - plt_rng/2.
      x2 = plt_pos + plt_rng/2.
      dx = x2-x1
      dxr = 1./dx

      ;        for jt=0L, ntime-1 do begin
      for jt=0L, ntime_total-1 do begin

	if (eof(1) or eof(2)) then begin
		break
	endif
        readU, 1, head   &   readU, 1, time    &   readU, 1, head
        readU, 1, head   &   readU, 1, field   &   readU, 1, head
        readU, 2, head   &   readU, 2, time    &   readU, 2, head
        readU, 2, head   &   readU, 2, field2   &   readU, 2, head

        if (jjt Mod tskp) eq 0 then begin

          time_ar[jnt] = time

          ;div
          divf = DIVFIELD(field[*], field2[*])
          w1 = divf[*,0] ; for y
          w2 = divf[*,1] ; for z
          ;w1 = divf[*,2] ; back y
          ;w2 = divf[*,3] ; back z

          bw_x =  sqrt( w1[ x1:x2 ]^2 + w2[ x1:x2 ]^2 )
          field_av[jnt] = total( bw_x ) * dxr

          jnt = ++ jnt
        endif    ; skip

        jjt = ++ jjt

      endfor     ; time loop

      Close, 1   &   Close, 2   &   Close, 3
      field_final=[field_final,field_av];
      time_final=[time_final,time_ar];
      jnt = 0;
    endfor     ; job loop
    
    if (i eq 0) then begin
      cgWindow, woxmargin=[5,5], woymargin=[5,5], wxsize=1080, wysize=960
      cgplot, time_final, field_final, color=colors[i], ytitle='amplitude',xtitle='time',title='forward field at position'+strtrim(plt_pos)+' with a range of'+strtrim(plt_rng), /window
    endif
    cgplot, time_final, field_final, color=colors[i], /overplot, /addcmd
  endfor
  AL_Legend,names[1:cases],LineStyle=intarr(cases),color=colors[0:cases-1],/window
  
endif else begin
  INPUT_FILE, jobname, prefname, firnum, endnum, njob
  
  ; Reading files
  jfrm = 0
  jjt = 0L
  jnt = 0L

  time_final = [];
  field_final = [];
  for ijob = 0, njob-1 do begin    ; job loop

    fcomp1 = '../../dat/' +prefname +'/' +jobName[ijob] + fcomp[0]
    fcomp2 = '../../dat/' +prefname +'/' +jobName[ijob] + fcomp[1]

    ; Read input parameters in KEMPO from .input_idl file
    prm_file = '../../dat/' +prefname +'/'  + jobName(ijob) + '.prm'
    READ_KEMPOPRM, prm_file, kmp

    nprev = (kmp.jobnum - 1) * kmp.nstep/ kmp.jobnum
    ntime_total = (kmp.nstep - nprev)/ kmp.ifdiag
    ntime = ntime_total / tskp
    field_av = Fltarr(ntime)
    time_ar = Fltarr(ntime)

    if ( (wfield eq 3) or (wfield eq 4) ) then  $
      ntime = (kmp.nstep - nprev)/ kmp.mltstp

    if (ijob eq 0) then begin
      field = Fltarr(kmp.nx)
      field2 = field[*]
      ttmp  = Fltarr(1)
      ftmp1 = Fltarr(kmp.nx)
      ftmp2 = ftmp1[*]
      x_ax = ( Findgen(kmp.nx) - (kmp.nx-1)/2. ) * kmp.dr/ kmp.cv
      xcent = kmp.nx/ 2

      ws = 2.*!dpi/ (kmp.dt*kmp.ifdiag*tskp)
      print, 'samplinf freq: ', ws

    endif

    ;--- Open files ---
    print, 'Opening ', fcomp1   &   Openr, 1, fcomp1   &   Point_lun, 1, 0
    print, 'Opening ', fcomp2   &   Openr, 2, fcomp2   &   Point_lun, 2, 0


    ;----- Time loop in each job -----

    head = var.head

    x1 = plt_pos - plt_rng/2.
    x2 = plt_pos + plt_rng/2.
    dx = x2-x1
    dxr = 1./dx

    ;        for jt=0L, ntime-1 do begin
    for jt=0L, ntime_total-1 do begin

	if (eof(1) or eof(2)) then begin
		break
	endif
      readU, 1, head   &   readU, 1, time    &   readU, 1, head
      readU, 1, head   &   readU, 1, field   &   readU, 1, head
      readU, 2, head   &   readU, 2, time    &   readU, 2, head
      readU, 2, head   &   readU, 2, field2   &   readU, 2, head

      if (jjt Mod tskp) eq 0 then begin

        time_ar[jnt] = time

        ;div
        divf = DIVFIELD(field[*], field2[*])
        w1 = divf[*,0] ; for y
        w2 = divf[*,1] ; for z
        ;w1 = divf[*,2] ; back y
        ;w2 = divf[*,3] ; back z

        bw_x =  sqrt( w1[ x1:x2 ]^2 + w2[ x1:x2 ]^2 )
        field_av[jnt] = total( bw_x ) * dxr

        jnt = ++ jnt
      endif    ; skip

      jjt = ++ jjt

    endfor     ; time loop

    Close, 1   &   Close, 2   &   Close, 3
    field_final=[field_final,field_av];
    time_final=[time_final,time_ar];
    jnt = 0;
  endfor     ; job loop
  
  cgWindow, woxmargin=[5,5], woymargin=[5,5], wxsize=1080, wysize=960
  cgplot, time_final, field_final, color=colors[0], ytitle='amplitude',xtitle='time',title='forward field at position'+strtrim(plt_pos)+' with a range of'+strtrim(plt_rng), /window
  AL_Legend,prefname,LineStyle=0,color=colors[0], /window
endelse


end
