; kt -> gn, fre, dwdt, S
ssize = size(gn)
width = ssize[1]
hight = ssize[2]
nnx = width / 2
gn_kt = fltarr(nnx,hight)
ww_kt = fltarr(nnx,hight)
dwdt_kt = fltarr(nnx,hight)
ss_kt = fltarr(nnx,hight)
for ii=0,hight-1 do begin
	gn_kt_tmp = fft(gn[*,ii],-1)
	gn_kt[*,ii] = 2.d * gn_kt_tmp[0:nnx-1]

	ww_kt_tmp = fft(forf[*,ii],-1)
	ww_kt[*,ii] = 2.d * ww_kt_tmp[0:nnx-1]

	dwdt_kt_tmp = fft(dwdt[*,ii],-1)
	dwdt_kt[*,ii] = 2.d * dwdt_kt_tmp[0:nnx-1]

	ss_kt_tmp = fft(S[*,ii],-1)
	ss_kt[*,ii] = 2.d * ss_kt_tmp[0:nnx-1]
endfor



end
