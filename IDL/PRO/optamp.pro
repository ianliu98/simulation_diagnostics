;DEVICE, DECOMPOSED=0
;DEVICE, RETAIN=2

    ; Q parameter
    we  = 1.0
;    qt = 0.5
    qt = 0.5
    bmax = 1.001
;    bmax = 1.03
;    bmax = 1.1

;    vthpa = 35.
;    vthpe = 47.   ; = vpe0
    vthpa = 28.
    vthpe = 49.   ; = vpe0

    wp1  = 5.0              ; Plasma frequency (cold)
    wp2 = 0.4           ; Plasma frequency (hot)

     cv = 100.0		;*******************

    nx  = 32768L

    w = FINDGEN(1000)*.001 + .001
    w = 0.4

    ; Local Gyrofrequency of electron
    qm  = 1.0
    wc  = 1.0
 ;--------------------------------------


; Threshold Bth   [Omura2009]

     kmp = {cv:cv, wc:wc, wp1:wp1, wp2:wp2, nx:nx}
     in = {qt:qt, omg:w, vpa:vthpa, vpe:vthpe, bmax:bmax, we:we}
     out = PLASMA_PRM(5, kmp, in=in)

     bth = out.bth[*] * qm


; Optimum amplitude   +++ At the equator +++

     tau = [.25, .5, 1., 2.]
     bopt = FLTARR(N_ELEMENTS(w), N_ELEMENTS(tau))
     kmp = {wp1:wp1, wp2:wp2, wc:wc, cv:cv}

     FOR jtau = 0, N_elements(tau)-1 do begin

        in = {qt:qt, tau:tau[jtau], omg:w, vpa:vthpa, vpe:vthpe, we:we}
        out = PLASMA_PRM(6, kmp, in=in)
        bopt[*,jtau] = out.bopt[*]
     ENDFOR


; ---------------------------------------------
; Plot

     ; Get global variables
     GLBVAR, var

     LOADCT, 39
     !p.charsize=3
     coli = (INDGEN(N_ELEMENTS(tau)-1)+1) * 254/(N_ELEMENTS(tau)-1)

     ps = 0
     SET_DISP, ps, 'chorus.eps'
     if (ps eq 1) then begin
        !P.charsize = 2
        !P.thick = 3
	    print, 'Output the EPS file'
     endif
     if (ps eq 0) then begin
        !P.charsize = 3
        !P.thick = 2
     endif

     ymin = MIN(bth[*]/10)

     ; optimum amplitude
     cgPLOT, w(*), bopt(*,0), /YLOG,  $
     ; PLOT, w(*), tran(*,0), /YLOG,  $
           xrange=[.001, .95], yrange=[ymin, 1.e-2],  $
           XGRIDSTYLE=1, YGRIDSTYLE=1, XTICKLEN=1, YTICKLEN=1, xstyle=1, ystyle=1,  $
           XTITLE= var.omgc +'     ['+ var.lomgc+ '!De0!N]', YTITLE='B!Dw!N     [B!D0eq!N]'
     ywindow = !Y.window

     ; tau0
        cgtext, 0.77, ywindow[1]-0.1, STRING(tau(0), format='(f4.2)'), /NORMAL




     FOR it=1, N_ELEMENTS(tau)-1 DO BEGIN
        cgPLOT, w(*), bopt(*,it), COLOR=coli(it-1), /overplot
        ; OPLOT, w(*), tran(*,it), COLOR=coli(it-1)
     ; tau
        XYOUTS, 0.77, ywindow[1]-0.1-.04*it, STRING(tau[it], format='(f4.2)'), COLOR=coli[it-1], /NORMAL
     ENDFOR

     ; threshold
     cgPLOT, w(*), bth(*), LINESTYLE=2, /overplot

     ; Q
     cgtext, 0.7, ywindow[1]-0.05, 'Q = ' +STRING(qt, format='(f4.2)'), /NORMAL
     cgtext, 0.7, ywindow[1]-0.1, var.tau +'  =', /NORMAL

     IF (ps EQ 1) THEN  DEVICE, /CLOSE
END
