;
; plasma_prm.pro
;
; PURPOSE:
;     Output various plasma parameters (k, Bth, Vtr, ...) related to chorus
;
; INPUT:
;    in: (STRUCTURE)
;       vpa :
;       vpe :
;       cv : speed of light
;       nx : grid
;       position : now position
;       wp1: cold electron plasma frequency
;       wp2: hot electron plasma frequency
;       bw : Bw of wave
;       bw0 : Bw of wave at the equator
;
;    kmp:
;       (STRUCTURE)
;       should be called in current source, in advance.
;          by READ_KEMPOPRM, prm_file, kmp
;
;    ver: type of plot
;       (INTEGER)
;         = 2: Energy curves
;           5: Threshold amplitude (function of frquencies) [Omura2011]
;           6: Optimum amplitude (function of frquencies) [Omura2011]
;          10: wavenumber (function of frequency)
;          43: Resonance velocity
;                 function of frquency, Vpara and Vperp are unique value
;          44: Resonance curves
;                 function of frquencies and Vperps
;          45: Diffusion curves
;                 function of frquency and Vperps
;                 
;
; OUTPUT:
;    out: structure that the following parameters are included
;
; ------------------
; Sample:
;    out = plasma_prm(ver, kmp, in=in) 
; ------------------
;
;                                                        Nov 20, 2012 M.Hikishima @Sendai
;-

FUNCTION  plasma_prm, ver, kmp, IN=in


     ; load constant quantities
     GLBVAR, var


;     pos = -9999
;     bw  = -9999
;     bw0 = -9999
;     vpe0 = -9999
;     vpa = -9999
;     vpe = -9999

     ; Global canstant
     qm = var.qm          ;  ratio of charge to mass

     if ( (ver eq 5) or (ver eq 6) or (ver eq 10) or (ver eq 43) or (ver eq 44) ) then begin
        we0   = Abs( kmp.wc )
        we_n  = in.we/ we0
        b0    = we0/ Abs(qm)
        omg_n = in.omg/ Abs(kmp.wc)
     endif

     if ( (ver eq 6) or (ver eq 5) ) then  begin
         wp2_n  = kmp.wp2/ we0
         vpe0   = Double(in.vpe)
         vpe0_n = vpe0/ kmp.cv
     endif


;     bmax = kmp.b02

;        pos  = in.pos
;        bw   = in.bw
;        bw0  = in.bw0
;     ww  = bw  * qm      ; Gyrofrequency related to wave amplitude (local)
;     ww0 = bw0 * qm      ; Gyrofrequency related to wave amplitude (equator)


     if (ver eq 5) then begin 
        ;--- Coefficient of nonuniform B0, B = B0(1+ag*x^2)  (x = x - nx/2) ---
        ;+++ In symmetric B0 +++
        ag = (in.bmax/ b0 - 1.d0)/ (kmp.nx/ 2)^2
;        print, 'a (real) :', ag
        ; normalized
        ag_n = ag * (kmp.cv/we0)^2
        print, 'a (sim) :', ag_n
     endif


     if ( (ver eq 5) or (ver eq 6) or (ver eq 10) or (ver eq 43) or (ver eq 44) or (ver eq 45) ) then begin

        xi2    = in.omg * (in.we - in.omg)/ (kmp.wp1^2+kmp.wp2^2)

        delta2 = 1.d0/ (1.d0 + xi2)
        xi     = Sqrt( xi2 )
        delta  = Sqrt( delta2 )

        if ( (ver eq 5) or (ver eq 6) or (ver eq 43) ) then begin
           ;--- Relativistic factor ---
           gm = 1.d0/ Sqrt( (1.d0 - (in.vpa^2 + in.vpe^2)/ kmp.cv^2) )
        endif 

     endif 

;--------------------------------
;        Wavenumber
;--------------------------------
     if (ver eq 10) then begin
        wavenumber = in.omg[*]/ (kmp.cv * delta[*] * xi[*])
        wavenumber_n = wavenumber[*] * kmp.cv
        out = { k_n: wavenumber_n[*] } 
     endif


;--------------------------------
;         Phase velocity
;--------------------------------
     if ( (ver eq 5) or (ver eq 6) or (ver eq 44) or (ver eq 45) ) then begin
          vp   = kmp.cv * delta * xi
          vp_n = vp/ kmp.cv
     endif


;--------------------------------
;       Group velocity
;--------------------------------
     IF ( (ver eq 5) or (ver eq 6) ) then begin

          ; INPUT:     kmp.cv
          ;            we
          ;            in.omg

          vg   = kmp.cv * xi/ delta * ( xi2 + in.we/( 2.0*(in.we - in.omg) ) )^(-1)
          vg_n = vg/ kmp.cv

     ENDIF


;--- Resonance velocity ---
    IF ( (ver eq 5) or (ver eq 6) or (ver eq 43) ) then begin

       ; INPUT: vpa, vpe, we, omg
       ; OUTPUT: vr_n (Vpara)

       vr   = kmp.cv * delta * xi * (1.d0 - in.we/ (gm * in.omg))
;print,'delta, xi, kmp.cv, in.omg, in.we, gm, vr: ',delta, xi, kmp.cv, in.omg, in.we, gm, vr
;stop
       vr_n = vr/ kmp.cv

       out = { vr:vr_n }
    ENDIF

;--- s? parameters ---
     IF ( (ver eq 5) or (ver eq 6) ) then begin


          s0 = delta/ xi * vpe0_n
          s1 = gm * (1.0 - vr/ vg)^2
          gd = 1.0        ; Plasma density coefficient along field line
          s2 = 1.0/ (2.0* xi * delta) * ( gm * in.omg/in.we * (vpe0/kmp.cv)^2  $
               - (2.0 + gd * delta2*(in.we - gm*in.omg)/(in.we - in.omg)) * (vr_n * vp_n) )

     ENDIF


;---------------------------------
;       Thermal momenta 
;---------------------------------
     IF ( (ver eq 5) or (ver eq 6) ) then begin
          upa   = in.vpa * gm
          upe   = in.vpe * gm
          upa_n = upa/ kmp.cv
          upe_n = upe/ kmp.cv
     ENDIF


;---------------------------------
;      Resonance curves 
;---------------------------------
     IF (ver eq 44) then begin

        ; INPUT    in.omg :  1 Frequency
        ;          in.we  :  Local we
        ;          in.vpe :  Perpendicular velocities (normalized, array)
        ;
        ; OUTPUT   out.vr: array of -Vpara(s) for +k

        wwp  = omg_n^2 + we_n^2 * vp_n^2
        nume = omg_n^2 - Sqrt( omg_n^4 - wwp * (omg_n^2 - we_n^2 * (1.0 - in.vpe[*]^2)) )

print, 'replaced with "resv.pro"'
stop


        vr_n = vp_n * nume[*]/ wwp

        out = {vr_n: vr_n}

     ENDIF


;---------------------------------
;       Diffusion curves
;---------------------------------
    IF (ver eq 45) then begin
        ; INPUT    in.k:  sign of k vector, 1: +k, -1: -k
        ;          in.vpe_n_at_vp: array(s) of vpe (normalized) at phase velocity vp
        ;
        ; OUTPUT   out.vpa_n_diff: array of Vpara(s) and the number of curve(s)
 
        dv = 0.01
        nv = 1.0/dv
        vpa_n_diff = Dblarr(nv+1, N_elements(in.vpe_n_at_vp))
        vpe_n_diff = Dblarr(nv+1, N_elements(in.vpe_n_at_vp))
        
        for j=0, N_elements(in.vpe_n_at_vp)-1 do begin
           vpe_n_arr       = Dindgen(nv+1) * dv
           vpe_max_arr     = vpe_n_arr[*] * in.vpe_n_at_vp[j]   ; array that vpe_max is maximum
           vpa_n_diff[*,j] = sqrt(in.vpe_n_at_vp[j]^2 - vpe_max_arr[*]^2)
           vpe_n_diff[*,j] = vpe_max_arr[*]
        endfor

        ; If Vp is positive, then -Vpa is used, v.v.
        if (in.k eq 1) then  $
           vpa_n_diff =   vp_n - vpa_n_diff[*,*]
        if (in.k eq -1) then  $
           vpa_n_diff = - vp_n + vpa_n_diff[*,*]

        out = { vpa_n_diff:    vpa_n_diff[*,*],  $
                vpe_n_diff_ax: vpe_n_diff[*,*] }

    ENDIF

;---------------------------------
;       Energy curves 
;---------------------------------
    IF (ver eq 2) then begin

        ; INPUT    in.ke: array(s) of energies [eV] to plot
        ;          in.dv: resolution of velocity
        ;
        ; OUTPUT   ou.vpe: array(s) of Vperp(s)

        in.ke = var.qe * in.ke     ; eV to Joule
        vpa   = Findgen(2/in.dv+1) * in.dv - 1.0
        vpe   = Fltarr(2/in.dv+1, N_elements(in.ke)+1)     ; + cv curve
        for je=0, N_elements(in.ke)-1 do  $
           vpe[*, je] = Sqrt( 1.0 - ( in.ke[je]/(var.me * var.cvv^2.) + 1.0)^(-2) - vpa[*]^2 )
        vpe[*, N_elements(in.ke)] = Sqrt( 1.0 - vpa[*]^2 )     ; cv curve

        out = {vpa: vpa, vpe: vpe}

     ENDIF



;--- Threshold of amplitude Bth [Omura2011] ---
     IF (ver eq 5) then begin

          ; INPUT    in.vpa: thermal parallel velocity
          ;          in.vpe: thermal perpendicular velocity
          ;          in.omg; wave frequency
          ;          in.qt;  depth of electron hole 
          ;          kmp.cv
          ;          kmp.wp2

          bth = 100.0 * !pi^3 * gm^3 * xi/  $
                ( omg_n * wp2_n^4 * vpe0_n^5 * delta^5 )  $
                * ( ag_n * s2 * upa_n/ in.qt )^2  $
                * Exp( gm^2 * vr_n^2/ upa_n^2 )


          out = {bth:bth}
     ENDIF


;--- Optimum amplitude ---
     IF (ver eq 6) then begin

          ; INPUT    in.vpa: thermal parallel velocity
          ;          in.vpe: thermal perpendicular velocity
          ;          in.omg; wave frequency
          ;          in.qt;  depth of electron hole 
          ;          kmp.cv
          ;          kmp.wp2

          bopt = 0.81 * !pi^(-2.5) * in.qt/ in.tau  $
                    * s1 * vg_n/ (s0 * omg_n * upa_n)  $
                    * (wp2_n * vpe0_n * delta/ gm)^2  $
                    * Exp(- gm^2 * vr_n^2/ (2.0 * upa_n^2))

          out = {bopt:bopt}

     ENDIF


;; Trapping velocity   ; classic version
;
;     ;???????????????
;     vtr = 2.0 * Sqrt(in.vpe * ww * kmp.cv * dlt * xi/ in.omg)
;
;     ; normalized
;     vtr_n = vtr/ kmp.cv


;; Trapping frequency
;
;     wt  = Sqrt( in.vpe * in.omg * ww/ (kmp.cv * dlt * xi) )
;     wtr = wt * dlt * gm^(-0.5)




;; dw/dt at the equator, S = -0.4 [Omura2009]
;
;        dwdt = 0.4 * s0 * omg/ s1 * ww0
;
;        PRINT, 'dw/dt -equator- : ', dwdt



; Distance hc where nonlinear growth is sustained [Omura2009]
;
;     hc = s0 * in.omg * ww0/ (5.0 * kmp.cv * ag * s2 * we0)
;
;     ; normalized
;     hc_n = hc/ (kmp.cv/we0)


;; Inhomogeneity ratio S
;
;        ; dw/dt
;        dwdtt = 6.e-5
;
;        dwdh = 2. * ag * (pos-kmp.nx2) * we0     ; dwdh < 0 at hh < nx/2
;        a1 = s1 * dwdtt
;        a2 = kmp.cv * s2 * dwdh
;        Sv = - 1./(s0 * omg * ww) * (a1 + a2)
;
;        PRINT, 'S : ', Sv
;        PRINT, 'a1, a2 : ', a1, a2






;     print, 'Local B0 : ', b0_loc
;     print, 'Coefficient of nonuniform B0 : ', ag_n
;     print, 'Critical distance, hc : ', hc_n
;     print, 'Relativistic factor, gm : ', out.gm
;     print, 'Vt, Vtr, Wt : ', vt_n, vtr_n, wt
;     print, 'Vtr, Wt : ', vtr_n, wt
;     print, 'Vr, Vp, Vg : ', vr_n, vp_n, vg_n
;     print, 'omg, k, lambda : '
;     for i=0, N_elements(in.omg)-1 do  $
;        print, in.omg[i], in.omg[i]/ vp_n[i], 2.d0*!dpi/ (in.omg[i]/vp_n[i])
;     print, 'Threshold of Bw, Bth : ', bth


     RETURN, out

END
