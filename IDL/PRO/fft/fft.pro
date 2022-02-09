
cgdisplay, 1200, 900
!p.multi = [0,1,4]
!p.charsize = 2

fs = 65536
dt = 1./ fs
len = 2048
freq1 = 450 
freq2 = 50
df = fs/ len
time = findgen(len) * dt
freq = findgen(fs) * df


iplot = 3


;; FFT(one wave) -> IFFT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if iplot eq 1 then begin

y = sin( 2.*!pi* freq1 * time)
plot, y

; FFT
ffty = fft(y, -1) ; complex
plot, freq[0:30], abs( ffty[0:30] ), xtitle='Hz'

; IFFT 
iffty = FFT( ffty, 1 ) ; complex
plot, real_part(iffty)
endif


; FFT(two wave) -> IFFT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if iplot eq 2 then begin

y = sin( 2.*!pi* freq1 * time) + sin( 2.*!pi* freq2 * time)
plot, y

; FFT
ffty = fft(y, -1) ; complex
plot, freq[0:30], abs( ffty[0:30] ), xtitle='Hz'

; IFFT 
iffty = FFT( ffty, 1 ) ; complex
plot, real_part(iffty)
endif


; FFT(one wave) -> bandpass -> IFFT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if iplot eq 3 then begin

y = sin( 2.*!pi* freq1 * time)
plot, y

; FFT
ffty = fft(y, -1) ; complex
plot, freq, abs( ffty ), xtitle='Hz', xrange=[0,2000]

; bandpass (negative=0 and power*2)
ffty[0:14] = 0.
ffty[len/2:*] = 0.
ffty *= 2
plot, freq, abs( ffty ), xtitle='Hz', xrange=[0,2000]

; IFFT 
iffty = FFT( ffty, 1 ) ; complex
plot, real_part(iffty)
endif




end
