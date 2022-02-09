;
;     calculate energy
;
;-

    cv = 3.d8			; [m]
    m0 = 9.10938291d-31		; [kg]
    qe = 1.602176565d-19	; [c]

    v = 0.01d *  cv

    gamma = 1.d0/ sqrt( 1.d0-v^2/cv^2)
    E     = m0/qe * cv^2*(gamma-1.d0)

    print, 'energy [eV] : ', E

END
