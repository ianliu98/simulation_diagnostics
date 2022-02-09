;                 
; PURPOSE:
;    Calculation time series of frequencies,
;       which are calculated from time difference of phase.
;
;    omega[i] = (wavephase[i+1] - wavephase[i]) / dt
;       angular frequency = radian / dt
;
; NOTE:
;    The input array is overwritten.
;
; OUTPUT:
;    1-d/ 2-d(when amplitude keyward set) array of frequencies
;
;                                         M.Hikishima  Dec 26, 2018 @ ISAS
;-

FUNCTION   OMG_CAL, $
           waveform_y, $    ;
           waveform_z, $    ;
           dt, $          ;
           AMPLITUDE=amplitude ; if set, amplitude=sqrt(y^2 + z^2)

     tlen = n_elements( waveform_y )
     freq = Fltarr(tlen, /nozero)


     pi2 = 2.*!dpi
     pi5_1 = !dpi/6
     pi5_0 = pi2 - pi5_1
     dt_inv = 1./ dt

     rad = Atan( waveform_z[*], waveform_y[*] ) 

     ; convert to 0 - 2pi, allow case of negative frequency
     rad = (rad[*] le 0) * pi2 + rad[*]

     for i=1, tlen-1 do begin

       radi = rad[i]

       ; case1. phases step across 0 deg.
       ; pi5_0 < radi < 2pi
       if (rad[i-1] ge pi5_0) then begin
         ; 0 < radi < pi5_1
         if ((radi le pi5_1) and (radi ge 0)) then begin
           radi += pi2
         endif
       endif

       ; case2. The phase is reversed (rad[i-1] > rad[i])
       if ((radi ge pi5_0) and (rad[i-1] le pi5_1)) then begin
         radi -= pi2 ; original
       endif

       freq[i] = (radi - rad[i-1]) * dt_inv 

     endfor

     freq[0] = !values.f_nan

     if KEYWORD_SET(amplitude) then begin
       amplitude = sqrt( waveform_y[*]^2 + waveform_z[*]^2 )
       amplitude[0] = !values.f_nan
       RETURN, [ [freq[*]], [amplitude[*]] ]
     endif else begin
       RETURN, freq[*]
     endelse

END
