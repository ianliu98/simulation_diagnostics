;----------------------------------------------------------------------------------
;
;                         Time - Space plot of field data
;
; OPTION:
;   - Separate into Forward and Backward waves,
;   - Distribution of wave frequency 
;       zero crossing method is used by a routine. 
;
; The data obtained by the kempo code is used
; Nov 27, 2009  M. Hikishima
;
;----------------------------------------------------------------------------------

; Shrink an array
tskip = 4
xskip = 4
;xskip = 16
pltt = 1


LOADCT, 39

; Global variables
     GLBVAR, var

; label
grad = ['homogeneous','1.87x10!U-7!N','3.74x10!U-6!N','3.74x!U-5!N']
ind_lab = ['a', 'b', 'c', 'd']

cmpr = 1
read, 'how many cases are involved? :', cmpr

for cmpr_ind=0,cmpr-1 do begin

; Read files
     INPUT_FILE, jobname, prefname, firnum, endnum, njob


; Plot option
     tmp = ['ew', 'bw']
     eb = tmp[1]

     tmp = ['notseparation', 'separation']
     Read, 'Separate into Northward and Southward ?, y[1], n[0] : ', i
     sep_mode = tmp[i]

     Read, 'Make frequency distribution ?, y[1], n[0] : ', freqdis

;     read, 'make elimination? y[1] n[0] : ', elim

; ---------------------------------------------------------------------------------

; Reading of files
jjt = 0L
FOR ijob = 0, njob-1 do begin

  datjob = '../../dat/' + prefname + '/' + jobname[ijob]
  prm_file = datjob +'.prm'
  READ_KEMPOPRM, prm_file, kmp
;  kmp.ifdiag = 256
  
  if (eb eq 'ew') then begin
    fname = datjob + ['.ey', '.ez']
    coe_eb = 1.0/ kmp.cv
  endif
  if (eb eq 'bw') then begin
    fname = datjob + ['.by', '.bz']
    coe_eb = 1.0
  endif

  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum 
  ntime = (kmp.nstep - nprev) /kmp.ifdiag 
;if (cmpr_ind eq 0) then begin
;  ntime = ntime/64
;endif
  IF (ijob eq 0) then begin
    total_time = ntime * njob

    ;--- Reduced size of time and space array ----------------------
    ; time
    rnt = total_time/ tskip
    fsample =  2.d0*!dpi/ (kmp.ifdiag * kmp.dt)
    fsample_r = fsample/ tskip
    ; Spatial
    ; Calculation of wavenumber to reduce the image resolution
    ; Take care the spatial Nyquist
    omg = 0.3
    b0 = 1

    rnx = kmp.nx/ xskip
    lambda = Wavelength( omg, kmp.wp1, kmp.wp2, kmp.cv, b0)
    if xskip ge lambda/2 then begin
      print, 'xskip is too larger'
      stop
    endif

    print, '>>>'
    print, '   Wavelength: ', lambda
    print, '   Spatial skip: ', xskip 
    print, '   Then, ' 
    print, '   The output array t-x: '
    print, total_time, kmp.nx, ' ->', rnt, rnx 
    print,''  &  print, '+++ sampling frequency +++ : ', fsample, ' Omega_e0'
    print, fsample, ' ->', fsample_r, ' Omega_e0'
    print,''  &  print,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    if (fsample_r le 3.0*Abs(kmp.wc)) then begin
      print, '!!! sampling frequency is small        !!!!!'
      stop
    endif
    start_time = (nprev + kmp.ifdiag) * kmp.dt
    if sep_mode eq 'notseparation' then begin
      field = Fltarr(rnx, rnt)   
      if (freqdis eq 1) then begin
        field1 = field
      endif
    endif else begin
      forf = Fltarr(rnx, rnt)
      bacf = forf
      forf_tmp = forf
      bacf_tmp = forf
      forf1_tmp = forf
      bacf1_tmp = forf
      if (freqdis eq 1) then begin
        forf1 = forf
        bacf1 = forf1
      endif
    endelse

    head = var.head
    ttmp = Fltarr(1)
    xtmp11 = Fltarr(kmp.nx)
    xtmp22 = xtmp11[*]
    t = Fltarr(rnt)
    midd = make_array(rnt,/integer,value=0)
    line = findgen(rnt)

    jx = Lindgen(rnx, increment=xskip, start=0)
    x_ax = ( Findgen(rnx)*xskip - (kmp.nx-1.)/ 2. ) * kmp.dr/ kmp.cv

    jnt = 0L
  endif 
  
  ; Open data file
  Openr, 2, fname[0]   &   Point_lun, 2, 0     ; beginning of the file
  Openr, 3, fname[1]   &   Point_lun, 3, 0
  print, 'Opening : ', fname

  for jt = 0L, ntime-1 do begin
    Readu, 2, head   &   Readu, 2, ttmp   &   Readu, 2, head
    Readu, 2, head   &   Readu, 2, xtmp11  &   Readu, 2, head
    Readu, 3, head   &   Readu, 3, ttmp   &   Readu, 3, head
    Readu, 3, head   &   Readu, 3, xtmp22  &   Readu, 3, head

    if ( (jjt MOD tskip) eq 0 ) then begin
      t[jnt] = ttmp
      xtmp1 = xtmp11[ jx[*] ] * coe_eb
      xtmp2 = xtmp22[ jx[*] ] * coe_eb

      ; For separation
      if sep_mode eq 'separation' then begin  
        divf = DIVFIELD(xtmp1[*], xtmp2[*])
        divf = Float( divf[*,*] )
        if freqdis eq 0 then begin
          forf[*, jnt] = Alog10( Sqrt(divf[*,0]^2 + divf[*,1]^2 ) > 1.E-20 )
          bacf[*, jnt] = Alog10( Sqrt(divf[*,2]^2 + divf[*,3]^2 ) > 1.E-20 )
;          if elim eq 1 then begin
;            forf_tmp[*, jnt] = divf[*,0] ; y
;            forf1_tmp[*, jnt] = divf[*,1] ; z
;            bacf_tmp[*, jnt] = divf[*,2] ; y
;            bacf1_tmp[*, jnt] = divf[*,3] ; z
;          endif
        endif else begin
;          forf_tmp[*, jnt] = Alog10( Sqrt(divf[*,0]^2 + divf[*,1]^2 ) > 1.E-20 )
;          bacf_tmp[*, jnt] = Alog10( Sqrt(divf[*,2]^2 + divf[*,3]^2 ) > 1.E-20 )  
          forf[*, jnt] = divf[*,0] ; y
          forf1[*, jnt] = divf[*,1] ; z
          bacf[*, jnt] = divf[*,2] ; y
          bacf1[*, jnt] = divf[*,3] ; z
        endelse

      endif else begin
      ; For not separation
        if freqdis eq 0 then begin
          field[*, jnt] = Alog10( Sqrt(xtmp1[*]^2 + xtmp2[*]^2) > 1.e-20)
        endif else begin
          field[*, jnt] = xtmp1[*]
          field1[*, jnt] = xtmp2[*]
        endelse

      endelse

      jnt = ++jnt
    endif     ; mod

    jjt = ++jjt
    IF ((jjt Mod 10000) EQ 0) THEN  PRINT, 'Reading.., time : ',jjt, '/', total_time
  ENDFOR
  CLOSE, 1   &   CLOSE, 2   &   CLOSE, 3
ENDFOR

;if freqdis eq 0 then begin
;  if elim eq 1 then begin
;    dt = tskip * kmp.ifdiag * kmp.dt
;    trans_forf = Transpose(forf_tmp[*,*])
;    trans_forf1 = Transpose(forf1_tmp[*,*])
;    trans_bacf = Transpose(bacf_tmp[*,*])
;    trans_bacf1 = Transpose(bacf1_tmp[*,*])
;    forf_tmp = trans_forf ; for store
;    bacf_tmp = trans_bacf ; for store

    ; forf, bacf are overwritten
;    for jx=0, rnx-1 do begin
;      forf_tmp[*, jx] = Omg_cal( trans_forf[*,jx], trans_forf1[*,jx], dt )
;      bacf_tmp[*, jx] = Omg_cal( trans_bacf[*,jx], trans_bacf1[*,jx], dt )
;    endfor
;    forf_tmp = Transpose(forf_tmp[*,*])
;    bacf_tmp = Transpose(bacf_tmp[*,*])
;  endif
;    forf[where(abs(forf_tmp) gt 0.10)] = -20.0
;    bacf[where(abs(bacf_tmp) gt 0.10)] = -20.0
;    forf[where(abs(forf_tmp) lt 0.01)] = -20.0
;    bacf[where(abs(bacf_tmp) lt 0.01)] = -20.0
;endif

; Making the frequency distribution (zero crossing) -----------------------
if freqdis eq 1 then begin

  dt = tskip * kmp.ifdiag * kmp.dt
  if sep_mode eq 'separation' then begin  

    trans_forf = Transpose(forf[*,*])
    trans_forf1 = Transpose(forf1[*,*])
    trans_bacf = Transpose(bacf[*,*])
    trans_bacf1 = Transpose(bacf1[*,*])
    forf = trans_forf ; for store
    bacf = trans_bacf ; for store

    ; forf, bacf are overwritten
    for jx=0, rnx-1 do begin
      forf[*, jx] = Omg_cal( trans_forf[*,jx], trans_forf1[*,jx], dt )
      bacf[*, jx] = Omg_cal( trans_bacf[*,jx], trans_bacf1[*,jx], dt )
    endfor
    forf = Transpose(forf[*,*])
    bacf = Transpose(bacf[*,*])
;    if elim eq 1 then begin
;      forf[where(forf_tmp lt -4)] = 0.0
;      bacf[where(bacf_tmp lt -4)] = 0.0
;    endif

  endif else begin

    trans_field = Transpose(field[*,*])
    trans_field1 = Transpose(field1[*,*])
    field = Transpose( field )
    for jx=0, rnx-1 do begin
      field[*, jx] = Omg_cal( trans_field[*,jx], trans_field1[*,jx], dt )
    endfor
    field = Transpose( field )
  endelse

endif


; ------------------------------------------------------------------------------
;                                  Plot
; ------------------------------------------------------------------------------


; resize
xmin = kmp.nxl/ xskip
;xmax = rnx/2 -1
;xmax = rnx -1
xmax = rnx -kmp.nxr/xskip -1
jtmin = 0
jtmax = rnt -1

if sep_mode eq 'notseparation' then begin
  if freqdis eq 0 THEN begin
    field = field[xmin:xmax, jtmin:jtmax]
  endif
endif
if sep_mode eq 'separation' then begin
;  if freqdis eq 0 then begin
    forf = forf[xmin:xmax, jtmin:jtmax]
    bacf = bacf[xmin:xmax, jtmin:jtmax]
;  endif
endif
x_ax = x_ax[xmin:xmax]
t = t[jtmin: jtmax]


;IF pltt eq 1 then begin

; plot setting
c_we0 = '['+var.lomgc+ '!S!De0!N]'
c_we0m = '['+var.lomgc+ '!S!De0!R!U-1!N]'
if freqdis eq 0 then begin
  if eb eq 'ew' then   bartitle='log!D10!N E!Dw!N/cB!D0!N'
  if eb eq 'bw' then   bartitle='log!D10!N B!Dw!N/B!D0!N'
endif else begin
  bartitle= var.omgc +'' + c_we0m
endelse

xyaxes = {xrange: [min(x_ax), max(x_ax)], yrange: [min(t), max(t)]/10000.0,  $
  xstyle: 1, ystyle: 1, ticklen: -0.01,  $
  xtitle: 'h [c' +var.lomgc+ '!S!De0!R!U-1!N]',  $
  ytitle: 't     x10!U4!N[' +var.lomgc+ '!S!De0!R!U-1!N]'}

if freqdis eq 0 then begin
  plmax = -3
  dlog  = 3
  plmin = plmax - dlog
endif ELSE begin
  plmax = 0.3
  plmin = 0.0
;  forf = abs(forf)
;  bacf = abs(bacf)
endelse

;while_tmp = 0L
;while (1 GT 0) DO BEGIN

  print, 'Max, Min (log plot): ', plmax, plmin
  print, '' & print, '--------------------------------------------------------'

  ; Separated ---------------------------------------------------------------
  IF (sep_mode eq 'separation') then begin   
    if (cmpr_ind eq 0) then begin
      !P.charsize = 2
      mmargin     = 1
      topd        = !D.table_size-2
      cgwindow, woxmargin = [10, 10], woymargin = [10, 10], wxsize = 1200, wysize = 960

      ;--- first ---
      cgImage, forf, /axes,  $
        minvalue=plmin, maxvalue=plmax, axkeywords=xyaxes, multimargin=mmargin, $
        top=topd, wmulti=[0,cmpr,1], font=1, charsize=2.5, /window
      cgplot, midd, line, linestyle=2, thick=3, color='white', /overplot,/addcmd
    endif else begin
      xyaxes.ytitle = ''
      xyaxes = Create_struct(xyaxes, 'ytickformat', "(A1)") 
      ;--- others ---
      cgImage, forf[*,*], /axes,  $
        minvalue=plmin, maxvalue=plmax, axkeywords=xyaxes, multimargin=mmargin, $
        top=topd, font=1, charsize=2.5, /addcmd
      cgplot, midd, line, linestyle=2, thick=3, color='white', /overplot,/addcmd
    endelse

    xwindow = !X.window
    ywindow = !Y.window
    cgText, xwindow[1]-0.25, ywindow[1]+0.02, '('+ind_lab[cmpr_ind]+')  ' + '  a = ' + grad[cmpr_ind] + '(' + var.lomgc + '!S!De0!R   /c)!U2!N', /NORMAL, $
            font=1,charsize=1.5, /addcmd
 ;   cgText, xwindow[0]+0.3, ywindow[0]-0.55, ind_lab[cmpr_ind], font=1,charsize=2,/addcmd

  endif else begin
  ; no separation ---------------------------------------------------------------
    !P.charsize = 2
    topd        = !D.table_size-2
    cgwindow, woxmargin = [10, 10], woymargin = [10, 10], wxsize = 1080, wysize = 1080
    cgImage, field[*,*], /axes,  $
      minvalue=plmin, maxvalue=plmax, axkeywords=xyaxes, top=topd, $
      pos=[0.2, 0.2, 0.7, 0.9], /window
  endelse

endfor
  ;-----------------------------------------------------------------------
  cgColorbar, /vertical, /right,  $
    yminor=1, range=[plmin, plmax], title=bartitle, ncolors=topd+1, font=1, charsize=2.5,  $
    position=[!x.window[1]+0.01, !y.window[0], !x.window[1]+0.02, !y.window[1]], /addcmd
    
;  READ, 'maxvalue adjustment? Yes[1], NO[0]', answer
;  if (answer eq 1) then begin
;    READ, 'maximum value: ', max_adjust
;    READ, 'minimum value: ', min_adjust
;    plmax = max_adjust
;    plmin = min_adjust
;  endif else begin
;    break
;  endelse
;  while_tmp = while_tmp + 1
;endwhile

;endif
  ;file_save = './save_file/'+prefname + '.sav'
  ;save, forf, filename=file_save
         
end
