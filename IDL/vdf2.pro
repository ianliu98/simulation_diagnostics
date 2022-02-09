GLBVAR, var

; Reading files
print, 'file 1:'
INPUT_FILE, jobname, prefname, firnum, endnum, njob
print, 'file 2:'
INPUT_FILE, jobname2, prefname2, firnum2, endnum2, njob2

; Constant
dir = '../../dat'

; resolution
read, 'resolution: defult[100] ', reso

; time point
read, 'which time point: ', tpoint

; overplot contour
ocont = 0

LOADCT,  33

for ijob=0, njob-1 do begin
  fname = jobname(ijob)
  pfile = FILEPATH(fname+'.fv2', ROOT_DIR=[dir], SUBDIR=[prefname] )
  prm_file = FILEPATH(fname+'.prm', ROOT_DIR=[dir], SUBDIR=[prefname] )
  READ_KEMPOPRM, prm_file, kmp
  
  fname2 = jobname2(ijob)
  pfile2 = FILEPATH(fname2+'.fv2', ROOT_DIR=[dir], SUBDIR=[prefname2] )
  prm_file2 = FILEPATH(fname2+'.prm', ROOT_DIR=[dir], SUBDIR=[prefname2] )
  READ_KEMPOPRM, prm_file2, kmp2
  
  nprev = (kmp.jobnum - 1) * kmp.nstep/ kmp.jobnum
  nt = (kmp.nstep - nprev)/ kmp.ivdiag2

  if ijob eq 0 then begin	
    nt  = nt + 1
	
	ttmp = Dblarr(1)
    xmi  = Dblarr(1)   &   xma = xmi
    ddvx  = Lonarr(1)
    vvmax  = Lonarr(1)
    tmpe  = Dblarr(reso+1)
    tmpa  = Dblarr(reso*2+1)
    fpape = DBLARR(reso*2+1, reso+1)
	fp_pape = fpape[*,*]
    vpa_ax    = Findgen(reso*2+1) - reso
    vpe_ax    = Findgen(reso+1)
    pe_n_ax = vpe_ax[*]/ reso
	
	ttmp2 = Dblarr(1)
    xmi2  = Dblarr(1)   &   xma2 = xmi2
    ddvx2  = Lonarr(1)
    vvmax2  = Lonarr(1)
    tmpe2  = Dblarr(reso+1)
    tmpa2  = Dblarr(reso*2+1)
    fpape2 = DBLARR(reso*2+1, reso+1)
    fpape_ini1 = fpape2[*,*]
    fpape_ini2 = fpape2[*,*]
	fp_pape2 = fpape2[*,*]
    vpa_ax2    = Findgen(reso*2+1) - reso
    vpe_ax2    = Findgen(reso+1)
    pe_n_ax2 = vpe_ax2[*]/ reso
	
  endif
  
  Openr, 2, pfile   &   POINT_LUN, 2, 0 
  print, ''  &  print, 'Open   ', pfile
  Openr, 1, pfile2   &   POINT_LUN, 1, 0 
  print, ''  &  print, 'Open   ', pfile2
  
  f1 = 0
  
  for jt=0, nt-1 do begin
	readu, 2, head  &  readu, 2, ttmp, xmi, xma, ddvx, vvmax  &  readu, 2, head
	readu, 1, head  &  readu, 1, ttmp2, xmi2, xma2, ddvx2, vvmax2  &  readu, 1, head
	
	if jt eq 0 then begin
      dvx   = ddvx[0]
      vmax = vvmax[0]
      fpah1 = Dblarr(kmp.nx/dvx+1, reso*2+1)
      fpah2 = fpah1(*,*)
      fpeh  = Dblarr(kmp.nx/dvx+1, reso+1)
      x_ax   = Findgen(kmp.nx/dvx+1) * dvx
	  
	  dvx2   = ddvx2[0]
      vmax2 = vvmax2[0]
      fpah12 = Dblarr(kmp.nx/dvx2+1, reso*2+1)
      fpah22 = fpah1(*,*)
      fpeh2  = Dblarr(kmp.nx/dvx2+1, reso+1)
      x_ax2   = Findgen(kmp.nx/dvx2+1) * dvx2
    endif
  
    ;;; f(Vpara,Vperp)
    for jpe=0, reso do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fpape[*, jpe] = tmpa[*]
    endfor
    for jpe=0, reso do begin 
      readu, 1, head  &  readu, 1, tmpa2  &  readu, 1, head
      fpape2[*, jpe] = tmpa2[*]
    endfor
  
    ;;; f(Vpara,h) 1
    for jx=0, kmp.nx/dvx do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fpah1[jx, *] = tmpa[*]
    endfor
    for jx=0, kmp2.nx/dvx do begin 
      readu, 1, head  &  readu, 1, tmpa2  &  readu, 1, head
      fpah12[jx, *] = tmpa2[*]
    endfor
  
    ;;; f(Vperp,h)
    for jx=0, kmp.nx/dvx do begin 
      readu, 2, head  &  readu, 2, tmpe  &  readu, 2, head
      fpeh[jx, *] = tmpe[*]
    endfor
	for jx=0, kmp2.nx/dvx do begin 
      readu, 1, head  &  readu, 1, tmpe2  &  readu, 1, head
      fpeh2[jx, *] = tmpe2[*]
    endfor
  
    ; f(Ppara,Pperp)
    for jpe=0, reso do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fp_pape[*, jpe] = tmpa[*]
    endfor
	for jpe=0, reso do begin 
      readu, 1, head  &  readu, 1, tmpa2  &  readu, 1, head
      fp_pape2[*, jpe] = tmpa2[*]
    endfor
	
	if (jt eq 0) then begin
      fpah01  = Total( fpah1(*,*) )
      fpah02  = Total( fpah2(*,*) )
      fpeh00  = Total( fpeh(*,*) )
      fpape02 = Total( fpape(*,*) )
      fpah01_inv = 1./ fpah01
      fpape02_inv = 1./ fpape02
	  
	  fpah012  = Total( fpah12(*,*) )
      fpah022  = Total( fpah22(*,*) )
      fpeh002  = Total( fpeh2(*,*) )
      fpape022 = Total( fpape2(*,*) )
      fpah01_inv2 = 1./ fpah012
      fpape02_inv2 = 1./ fpape022
    endif
	
	fpape = fpape02_inv * fpape[*,*]
    fpape = Alog10(fpape(*,*) > 1.e-20)
    fpah1 = fpah01_inv * fpah1[*,*]
    fpah1 = Alog10(fpah1(*,*)  > 1.e-20) 
    fpah2 = Alog10(fpah2(*,*)  > 1.e-20) 
    fpeh  = Alog10(fpeh(*,*)  > 1.e-20) 
	
	fpape2 = fpape02_inv2 * fpape2[*,*]
    fpape2 = Alog10(fpape2(*,*) > 1.e-20)
    fpah12 = fpah01_inv2 * fpah12[*,*]
    fpah12 = Alog10(fpah12(*,*)  > 1.e-20) 
    fpah22 = Alog10(fpah22(*,*)  > 1.e-20) 
    fpeh2  = Alog10(fpeh2(*,*)  > 1.e-20) 
	
	if f1 eq 0 then begin
      mafpape = FIX( MAX(fpape(*,*)) )
	  mafpape2 = FIX( MAX(fpape2(*,*)) )
      f1 = 1
    endif
    mifpape = mafpape-3
	mifpape2 = mafpape2-3

    mindf = -0.05
    maxdf = 0.05

	fpape[*,0] = !Values.f_nan
	vpa_min_g = 1
	vpa_max_g = reso*2+1
	vpe_min_g = 1
	vpe_max_g = reso+1
	
	vpa_min_ = ( vpa_ax[vpa_min_g-1]-0.5)/ reso * vmax/kmp.cv
	vpa_max_ = ( vpa_ax[vpa_max_g-1]+0.5)/ reso * vmax/kmp.cv
	vpe_min_ = ( vpe_ax[vpe_min_g-1]-0.5)/ reso * vmax/kmp.cv
	vpe_max_ = ( vpe_ax[vpe_max_g-1]+0.5)/ reso * vmax/kmp.cv
	
	vpa_min2_ = ( vpa_ax2[vpa_min_g-1]-0.5)/ reso * vmax2/kmp2.cv
	vpa_max2_ = ( vpa_ax2[vpa_max_g-1]+0.5)/ reso * vmax2/kmp2.cv
	vpe_min2_ = ( vpe_ax2[vpe_min_g-1]-0.5)/ reso * vmax2/kmp2.cv
	vpe_max2_ = ( vpe_ax2[vpe_max_g-1]+0.5)/ reso * vmax2/kmp2.cv
	
	
	
	xyaxes = { xrange: [vpa_min_, vpa_max_], yrange: [vpe_min_, vpe_max_],  $
                          ystyle: 1, ticklen: -0.01,  $
                          ytitle: 'v!D' +var.pec+ '!N/c' }
	
	xyaxes2 = { xrange: [vpa_min_, vpa_max_], yrange: [vpe_min_, vpe_max_],  $
                          xstyle: 1, ystyle: 1, ticklen: -0.01,  $
                          xtitle: 'v!D||!N/c', ytitle: 'v!D' +var.pec+ '!N/c' }
						  
	
	if (jt eq 0) then begin
	  !p.charsize = 2
	  mmargin = 4
	  topd = !D.table_size-2	
	  cgwindow, woxmargin = [5, 5], woymargin = [5, 5], wxsize = 1080, wysize = 1080
	endif
    	
    if (jt eq 0) then begin
      ; 1  
      fpape_ini1 = fpape
      cgImage, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1: vpe_max_g-1],  $
          /axes, $/interpolate,  $
          minvalue=mifpape, maxvalue=mafpape, axkeywords=xyaxes2, $
                  wmulti=[0,2,2], font=1, charsize=1, multimargin=mmargin, /window
      if (ocont eq 1) then begin
        ncontour = 15
        contourLevels =  $
            cgConLevels( fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1],  $
            NLevels=ncontour, MinValue=mifpape + 1 )
        cgContour, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
            Color='white', label=0, Levels=contourLevels
      endif

      cgText, !X.window[1]-0.35, !Y.window[1]+0.01,  $
           'h=-50[c'+var.lomgc+'!S!De0!R!U-1!N] ~ 50[c'+var.lomgc+'!S!De0!R!U-1!N]     t = ' + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
           /NORMAL, $
           font=1, charsize=1, /addcmd ; for movie
   
      ; 2
      fpape_ini2 = fpape2
      cgImage, fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1: vpe_max_g-1],  $
          /axes, $/interpolate,  $
          minvalue=mifpape2, maxvalue=mafpape2, axkeywords=xyaxes2, $
                  multimargin=mmargin, font=1, charsize=1, /addcmd
      if (ocont eq 1) then begin
        ncontour = 15
        contourLevels =  $
        cgConLevels( fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1],  $
            NLevels=ncontour, MinValue=mifpape2 + 1 )
        cgContour, fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
            Color='white', label=0, Levels=contourLevels
      endif

      cgText, !X.window[1]-0.35, !Y.window[1]+0.01,  $
           'h=50[c'+var.lomgc+'!S!De0!R!U-1!N] ~ 150[c'+var.lomgc+'!S!De0!R!U-1!N]     t = ' + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
           /NORMAL, $
           font=1, charsize=1, /addcmd ; for movie
      xwindow1 = !X.window
      ywindow1 = !Y.window

      cgColorbar, /vertical, /right,  $
      range=[mifpape2, mafpape2], ncolors=!D.table_size-1,  $
      position=[xwindow1[1]+0.01, ywindow1[0], xwindow1[1]+0.02, ywindow1[1]],  $
      title='log!D10!N f', font=1, charsize=2, /addcmd ; for movie

    endif

    ; 3
    if (jt eq tpoint) then begin
        fpape = (fpape - fpape_ini1) / fpape_ini1
	cgImage, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1: vpe_max_g-1],  $
          /axes, $/interpolate,  $
          minvalue=mindf, maxvalue=maxdf, axkeywords=xyaxes2, $
		  font=1, charsize=1, multimargin=mmargin, /addcmd
	; overplot contour
      if (ocont eq 1) then begin
        ncontour = 15
        contourLevels =  $
            cgConLevels( fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1],  $
            NLevels=ncontour, MinValue=mifpape + 1 )
        cgContour, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
            Color='white', label=0, Levels=contourLevels
      endif

      cgText, !X.window[1]-0.35, !Y.window[1]+0.01,  $
           'h=-50[c'+var.lomgc+'!S!De0!R!U-1!N] ~ 50[c'+var.lomgc+'!S!De0!R!U-1!N]     t = ' + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
           /NORMAL, $
           font=1, charsize=1, /addcmd ; for movie

    ; 4	
      fpape2 = (fpape2 - fpape_ini2) / fpape_ini2
      cgImage, fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1: vpe_max_g-1],  $
          /axes, $/interpolate,  $
          minvalue=mindf, maxvalue=maxdf, axkeywords=xyaxes2, $
		  multimargin=mmargin, font=1, charsize=1, /addcmd
      xwindow2 = !X.window
      ywindow2 = !Y.window
      cgColorbar, /vertical, /right,  $ 
          range=[mindf, maxdf], ncolors=!D.table_size-1,  $
          position=[xwindow2[1]+0.01, ywindow2[0], xwindow2[1]+0.02, ywindow2[1]],  $
          title='df/f!D0!N', font=1, charsize=2, /addcmd ; for movie
      ; overplot contour
      if (ocont eq 1) then begin
        ncontour = 15
        contourLevels =  $
        cgConLevels( fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1],  $
            NLevels=ncontour, MinValue=mifpape2 + 1 )
        cgContour, fpape2[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
            Color='white', label=0, Levels=contourLevels
      endif

      cgText, !X.window[1]-0.35, !Y.window[1]+0.01,  $
           'h=50[c'+var.lomgc+'!S!De0!R!U-1!N] ~ 150[c'+var.lomgc+'!S!De0!R!U-1!N]     t = ' + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]', $
           /NORMAL, $
           font=1, charsize=1, /addcmd ; for movie

      dummy = ''
      Read, 'Push enter', dummy
    endif

  endfor
  close,1 & close, 2
  
endfor

end
