;                 
;                                         M.Hikishima  Nov 14, 2012 @ sendai
;
; PURPOSE:
;    To divide wave field into Forward and Backward waves
;
; INPUT:
;    fy: y-axis perpendicular component of wave (e.g., Ey[0:nx-1])
;        1-d array in space
;    fz: z-axis perpendicular component of wave (e.g., Ez[0:nx-1])
;        1-d array in space
;
; OUTPUT:
;    array [nx x 4] in real component
;    [*,0]: Y component of forward wave
;    [*,1]: Z component of forward wave
;    [*,2]: Y component of backward wave
;    [*,3]: Z component of backward wave
;
;-


FUNCTION  divfield,  fy, fz

     nnx = N_elements( fy )


     ;--- FFT ---
     ify = FFT( Double(fy[*]), -1, /double )     ; complex array
     ifz = FFT( Double(fz[*]), -1, /double )     ; complex array


     ifyz = 0.5d * ( ify[*] * dComplex(0.0d, 1.0d)  +  ifz[*] )

     for_if  = ifyz[*]     ; for forward
     back_if  = ifyz[*]    ; for backward

     ;--- Nyquist component is not included ---
     ; DC(i=0) and Nyquist(i=nnx/2: nnx is even) are zero
     if ( (nnx mod 2) ne 0) then  stop
     for_if[nnx/2: *]  = 0.0d
     for_if[0]         = 0.0d
     back_if[0: nnx/2] = 0.0d

     ;--- power is doubled ---
     for_if  = 2.0d * for_if[*]
     back_if = 2.0d * back_if[*]

     ;--- IFFT ---
     for_f  = FFT(for_if, 1, /double)
     back_f = FFT(back_if, 1, /double)

     f_fbyz      = Dblarr(nnx, 4)
     f_fbyz[*,0] = Imaginary(for_f[*])      ; forward y
     f_fbyz[*,1] = Real_part(for_f[*])      ; forward z
     f_fbyz[*,2] = Imaginary(back_f[*])     ; backward y
     f_fbyz[*,3] = Real_part(back_f[*])     ; backward z

     RETURN, f_fbyz[*,*]
END
