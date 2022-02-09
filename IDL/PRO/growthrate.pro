;
; PURPOSE;
;     Calculate the Linear growthrate and Anisotropy [Kennenl & Petschek, 1966, JGR]
;
;          +++ a lot of number of particles will make smoothed results +++
;
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; *** CAUSION ***
;      In case of Uth is larger, the results (anisotropy, ...) may be different 
;         than expected one (relativistic effecct ?) 
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; INPUT:
;     dis: Velocity distribution function f(Vpara, Vperp) which has been normalized
;          the veloity ranges are -cv <= Vpara <= +cv and 0 <= Vperp <= +cv
;
; OUTPUT:
;     Array that contains the arraies
;        - omega
;        - anisotropy
;        - critical anisotropy
;        - growth rate
;
;                                                    M.Hikishima, 29 Dec, 2012 AT Sendai
;
;- 


FUNCTION growthrate, dis, we, kmp, model


ndis = Size(dis[*,*])
npa  = ndis[1]
npe  = ndis[2]

;-------------------------------------------------
; model distribution for check

     if (model eq 1) then begin
        print, '++++++++++ MODEL DISTRIBUTION IS USED ++++++++++'
        ; Isotropic is pa = pe (uth_pe is not Upe1(or Upe2))
        uth_pa = 0.10d
        uth_pe = 0.10d 
        for jpa=0, npa-1 do begin
           for jpe=0, npe-1 do begin
              upa = Double(jpa)/ ((npa-1)/2) - 1.d0
              upe = Double(jpe)/ (npe-1)
              dis[jpa, jpe] = 1.d0/ ( (2.d0*!dpi)^1.5 * uth_pa * uth_pe^2)  $
                              * Exp(-upa^2/(2.d0*uth_pa^2)) * Exp(-upe^2/(2.d0*uth_pe^2))
           endfor
        endfor
        dis = dis/ total(dis[*,*])
     endif
;-------------------------------------------------


; integration of the dis should be 1

     sum_dis = Total(dis[*,*])
     if ( Abs(Round(sum_dis) - sum_dis) ge 1.d-6) then begin
        print, 'integration is not 1'
        stop
     endif


print, 'Neglecting "relativity" and "displacement current"'
        print, '   and assume "Maxwellian"'

print,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print,'assumed + k only'
print,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'



;----------------------------------------------------------
;                    Calculation
;----------------------------------------------------------

; Definition of arraies

     dw    = 0.01d
     dv    = 1.d0/ (npe-1)
     nfreq = we/ dw
     omg   = Dindgen(nfreq) * dw + dw
     ani   = Make_array( nfreq, /double, value=!Values.d_nan )
     ac    = ani[*]
     gamma = ani[*]
     eta   = ani[*]
     ;dis  = Congrid(dis,101,51)
     pa_ax = Dindgen( npa )/ ( (npa-1)/2 ) - 1.d0 
     pe_ax = Dindgen( npe )/ (npe-1)


; Calculation of resonance velocity Vr

     in  = {we:we, omg:Double(omg), vpa:0.d0, vpe:0.d0}     ; non relativistic
     out = PLASMA_PRM( 43, kmp, in=in )

; Calculation


     ; F = f/ (dvpa dvpe), see Kennel
     dis = dis[*,*]/ dv^2

     for jf=0, nfreq-1  do begin

        jvpa = Fix( out.vr[jf]*100.0 + 100.0 )

        if ( (jvpa ge 1) and (jvpa le (npa - 2)) ) then begin
        ; jvpa >= 1, to make dis_l

           ; interpolating
           sl = ( out.vr[jf] - pa_ax[jvpa] )/ dv
           sr = 1.d0 - sl

           dis_vr = sr * dis[jvpa  ,*]  +  sl * dis[jvpa+1,*]
           dis_l  = sr * dis[jvpa-1,*]  +  sl * dis[jvpa  ,*]
           dis_r  = sr * dis[jvpa+1,*]  +  sl * dis[jvpa+2,*]

           ; Integration
           a_nume_int  = 0.d0
           a_denom_int = 0.d0
           eta_int     = 0.d0

           for jvpe=2, npe-2 do begin          ; * dis[*,0] = NAN

                a_nume_int  = pe_ax[jvpe] * dv  $
                              * ( out.vr[jf]   * (dis_vr[jvpe+1] - dis_vr[jvpe-1])/(2.d0*dv)  $
                                 - pe_ax[jvpe] * (dis_r[jvpe]    - dis_l[jvpe])   /(2.d0*dv) )  $
                              * pe_ax[jvpe]/ out.vr[jf]  $
                                 + a_nume_int

                a_denom_int = pe_ax[jvpe] * dv  *  dis_vr[jvpe]  $
                                 + a_denom_int

                eta_int     = pe_ax[jvpe] * dv  *  dis_vr[jvpe]  $
                                 + eta_int
           endfor

           ;eta
              eta[jf] = 2.d0*!dpi * (- out.vr[jf]) * eta_int[*]

           ; Anisotropy
              ani[jf] = a_nume_int[*]/ (2.d0 * a_denom_int[*])

           ; critical anisotropy
              ac[jf] = 1.d0 / (we/ omg[jf] - 1.d0)

           ; gamma
              gamma[jf] = !dpi * we * (1.d0 - omg[jf]/ we)^2  $
                          * eta[jf] * (ani[jf] - ac[jf])
        endif

     endfor   ; jf


     ; weighting of hot electrons
     gamma = gamma[*]  *  kmp.wp2^2/(kmp.wp1^2 + kmp.wp2^2)


     RETURN, [ [omg[*]], [ani[*]], [ac[*]], [gamma[*]] ]

END
