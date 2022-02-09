;+
; PURPOSE:
;    Common setting of display
;    X window or Postscritpt
;
; INPUT   
;    ps:
;       0: X window
;       1: Postxcript (eps)
;
;    psname:
;       name of postscript file (eps)
;
;
;    M.Hikishima  Sep 24, 2012
;-


PRO set_disp, ps, psname

    ; Global variables
    GLBVAR, var

    IF (ps ne 1) THEN BEGIN
        Set_plot, 'X'
        !P.font=1          ; True-type font
;        Window, xsize=var.wsz(0), ysize=var.wsz(1)
        ; back is white, line is white
        cgdisplay, xsize=var.wsz(0), ysize=var.wsz(1)
    ENDIF 
    IF (ps eq 1) THEN BEGIN
        print, 'file overwrite may cause a strange file.'
        Set_plot, 'PS'
        !P.font=0          ; Device font
        Device, /land, /color,  /encapsulated,  /helvetica,  $
           filename=psname, bits=8             , xoffset=24
    ENDIF

END
