GLBVAR, var

; Reading files
INPUT_FILE, jobname, prefname, firnum, endnum, njob

; Constant
dir = '../../dat'

; resolution
read, 'resolution: defult[100] ', reso

; time point
print, 'which time point to compare? '
read, 'compute by yourself:  ', tpoint

; overplot contour
ocont = 0

LOADCT,  33

for ijob=0, njob-1 do begin
  fname = jobname(ijob)
  pfile = FILEPATH(fname+'.fv2', ROOT_DIR=[dir], SUBDIR=[prefname] )
  prm_file = FILEPATH(fname+'.prm', ROOT_DIR=[dir], SUBDIR=[prefname] )
  READ_KEMPOPRM, prm_file, kmp
  
  nprev = (kmp.jobnum - 1) * kmp.nstep/ kmp.jobnum
  nt = (kmp.nstep - nprev)/ kmp.ivdiag2

  if ijob eq 0 then begin	

    nt  = nt + 1
    ttmp  = Dblarr(1)
    xmi   = Dblarr(1)   &   xma = xmi
    ddvx  = Lonarr(1)
    vvmax = Lonarr(1)
    tmpe  = Dblarr(reso+1)
    tmpa  = Dblarr(reso*2+1)
    fpape = DBLARR(reso*2+1, reso+1)
    fp_pape = fpape[*,*]
    fpape_ini = fpape[*,*]
    vpa_ax  = Findgen(reso*2+1) - reso
    vpe_ax  = Findgen(reso+1)
    pe_n_ax = vpe_ax[*]/ reso
	
  endif
  
  Openr, 2, pfile   &   POINT_LUN, 2, 0 
  print, ''  &  print, 'Open   ', pfile
  
  f1 = 0
  
  for jt=0, nt-1 do begin

    readu, 2, head  &  readu, 2, ttmp, xmi, xma, ddvx, vvmax  &  readu, 2, head
	
    if jt eq 0 then begin
      dvx   = ddvx[0]
      vmax  = vvmax[0]
      fpah1 = Dblarr(kmp.nx/dvx+1, reso*2+1)
      fpah2 = fpah1(*,*)
      fpeh  = Dblarr(kmp.nx/dvx+1, reso+1)
      x_ax  = Findgen(kmp.nx/dvx+1) * dvx
    endif
  
    ;;; f(Vpara,Vperp)
    for jpe=0, reso do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fpape[*, jpe] = tmpa[*]
    endfor
  
    ;;; f(Vpara,h) 
    for jx=0, kmp.nx/dvx do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fpah1[jx, *] = tmpa[*]
    endfor
  
    ;;; f(Vperp,h)
    for jx=0, kmp.nx/dvx do begin 
      readu, 2, head  &  readu, 2, tmpe  &  readu, 2, head
      fpeh[jx, *] = tmpe[*]
    endfor
  
    ; f(Ppara,Pperp)
    for jpe=0, reso do begin 
      readu, 2, head  &  readu, 2, tmpa  &  readu, 2, head
      fp_pape[*, jpe] = tmpa[*]
    endfor
	

    if (jt eq 0) then begin
      fpah01  = Total( fpah1(*,*) )
      fpah02  = Total( fpah2(*,*) )
      fpeh00  = Total( fpeh(*,*) )
      fpape02 = Total( fpape(*,*) )
      fpah01_inv = 1./ fpah01
      fpape02_inv = 1./ fpape02
    endif
	
    fpape = fpape02_inv * fpape[*,*]
    fpape = Alog10(fpape(*,*) > 1.e-20)
    fpah1 = fpah01_inv * fpah1[*,*]
    fpah1 = Alog10(fpah1(*,*)  > 1.e-20) 
    fpah2 = Alog10(fpah2(*,*)  > 1.e-20) 
    fpeh  = Alog10(fpeh(*,*)  > 1.e-20) 
	
	
    if f1 eq 0 then begin
      mafpape = FIX( MAX(fpape(*,*)) )
      f1 = 1
    endif
    mifpape = mafpape-3

    mindf = -0.2
    maxdf = 0.2

    fpape[*,0] = !Values.f_nan
    vpa_min_g = 1
    vpa_max_g = reso*2+1
    vpe_min_g = 1
    vpe_max_g = reso+1

    vpa_min_ = ( vpa_ax[vpa_min_g-1]-0.5)/ reso * vmax/kmp.cv
    vpa_max_ = ( vpa_ax[vpa_max_g-1]+0.5)/ reso * vmax/kmp.cv
    vpe_min_ = ( vpe_ax[vpe_min_g-1]-0.5)/ reso * vmax/kmp.cv
    vpe_max_ = ( vpe_ax[vpe_max_g-1]+0.5)/ reso * vmax/kmp.cv
	
	
    xyaxes = { xrange: [vpa_min_, vpa_max_], yrange: [vpe_min_, vpe_max_],  $
               ystyle: 1, ticklen: -0.01,                                   $
	       xtitle: 'v!D||!N/c',                                         $
               ytitle: 'v!D' +var.pec+ '!N/c' }
	
	
    xyaxes2 = { xrange: [vpa_min_, vpa_max_], yrange: [vpe_min_, vpe_max_],  $
               ystyle: 1, ticklen: -0.01,                                   $
               xtitle: 'v!D||!N/c'}

    if (jt eq 0) and (ijob eq 0) then begin
      !p.charsize = 2
      mmargin = [5, 9, 5, 9] ;[bottom, left, top, right]
      topd = !D.table_size-2	
      cgwindow, woxmargin = [5, 5], woymargin = [5, 5], wxsize = 1440, wysize = 720
    endif
    	
    if (jt eq 0) and (ijob eq 0) then begin
      ; initial  
      fpape_ini = fpape
      cgImage, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g: vpe_max_g-1],    $
               /axes, /interpolate, minvalue=mifpape, maxvalue=mafpape,    $
               axkeywords=xyaxes, wmulti=[0,2,1], font=1, charsize=1.5,    $
               multimargin=mmargin, /window

      ; overplot contour
      if (ocont eq 1) then begin
          ncontour = 15
          contourLevels =  $
              cgConLevels( fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g:vpe_max_g-1],      $
              NLevels=ncontour, MinValue=mifpape + 1 )
          cgContour, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
              Color='white', label=0, Levels=contourLevels, /addcmd
      endif

      cgText, !X.window[1]-0.25, !Y.window[1]+0.02,                                            $
             't = '                                                                            $
             + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]',                   $
             /NORMAL, font=1, charsize=1.5, /addcmd 

      cgColorbar, /vertical, /right,                                                           $
	         range=[mifpape, mafpape], ncolors=!D.table_size-1,                            $
		 position=[!x.window[1]+0.01, !y.window[0], !x.window[1]+0.02, !y.window[1]],  $
	         title='log!D10!N f', font=1, charsize=1.5, /addcmd 
		       
    endif

    if (jt eq tpoint) and (ijob eq njob-1) then begin
      ; compare
      fpape = fpape / fpape_ini - 1.0
      cgImage, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g: vpe_max_g-1],                 $
               /axes, /interpolate, minvalue=mindf, maxvalue=maxdf,                     $
               axkeywords=xyaxes2, font=1, charsize=1.5, multimargin=mmargin,           $
	       /addcmd

      ; overplot contour
      if (ocont eq 1) then begin
        ncontour = 15
        contourLevels =  $
            cgConLevels( fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g:vpe_max_g-1],      $
            NLevels=ncontour, MinValue=mifpape + 1 )
        cgContour, fpape[vpa_min_g-1: vpa_max_g-1, vpe_min_g-1:vpe_max_g-1], /OnImage,  $
            Color='white', label=0, Levels=contourLevels, /addcmd
      endif

      cgText, !X.window[1]-0.25, !Y.window[1]+0.02,                                            $
             't = '                                                                           $ 
	     + STRING(ttmp, '(F10.2)')+ '  [' +var.lomgc+ '!S!De0!R!U-1!N]',                   $
             /NORMAL, font=1, charsize=1.5, /addcmd ; for movie

      cgColorbar, /vertical, /right,                                                           $
	         range=[mindf, maxdf], ncolors=!D.table_size-1,                                $
		 position=[!x.window[1]+0.01, !y.window[0], !x.window[1]+0.02, !y.window[1]],  $
	         title='f/f!D0!N - 1', font=1, charsize=1.5, /addcmd 
    endif

  endfor

  close, 2
  
endfor

end
