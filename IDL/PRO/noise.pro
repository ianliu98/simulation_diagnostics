;                 
;                                         M.Hikishima  Dec 12, 2012 @ sendai
;
; PURPOSE:
;    To generate a noise signale in a finite frequency bandwidth with Gaussian or flat spectrum 
;
; INPUT:
;    * The frequencies is angular frequencies
;
;    t: time
;
;    wcen: a center frequency of the bandwidth 
;
;    wband: the bandwidth
; 
;    nwbin: the number of bins in the bandwidth (odd number)
;
;    ampc: amplitude of a spectra of the center frequency
;
; KEYWORD:
;    gaussian: make spectram in frequency domain
;
; OUTPUT:
;    noise: noise amplitude in any time
;
;-


FUNCTION  noise,  t, wcen, wband, nwbin, ampc, GAUSSIAN=gaussian

     if ( (nwbin mod 2) eq 0 ) then  stop
     ; better that wbin is about the same as dw in wt.pro
     dwbin = wband/ (nwbin-1)
     wbin = Findgen(nwbin) * dwbin + (wcen - wband/2)

     amp = Make_array( nwbin, value=ampc, /float)
     if Keyword_set(gaussian) then begin
        ; make gaussian amplitude; maximum amplitude (center) =1, sigma^2=1
        amp[*] = 0.0     ; reset
        xwid = 6.0     ; width of -3sigma to +3sigma
        xx   = Findgen(nwbin) * xwid/(nwbin-1) - xwid/2
        amp = 1.0/Sqrt(2.0*!pi) * Exp(- xx^2/2.0)
        amp = amp/ Max(amp)  *  ampc
     endif

     ; initial phase with random numbers
     ; *** this random number must be fixed ***
     seed = 1001L
     iniphi = Randomu(seed, nwbin) * 2.0*!pi

     noise = Total( amp[*] * Sin(wbin[*] * t  +  iniphi[*]) )


     RETURN, noise

END
