
!p.multi = [0,1,3]
!p.charsize = 2

fs = 65536
dt = 1./ fs
len = 1024
freq = 512 
df = fs/ len
time = findgen(len) * dt

;real
y = sin( 2.*!pi* freq * time)
plot, y

; FFT (complex)
ffty = fft(y, -1)
plot, abs( ffty )

;ffty[len/2:*] = 0

; IFFT (complex)
iffty = FFT( ffty, 1 )
plot, real_part(iffty)

end
