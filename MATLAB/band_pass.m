% function to band-pass field 
% a "tkip=4, xkip=4, ..." is necessary !!!

function pasf = band_pass(orif, wl, wh)

% parameters
tskip = 4;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;

% fft setting
ds = size(orif);
dts = ddt;
ws = 2*pi/dts;
fft_pnt_f = ds(1);
dw_f = ws / fft_pnt_f;
fss = (1:1:ds(1)) * dw_f - dw_f/2.0;

% band index
tmp1 = find(fss<wl);
if (isempty(tmp1))
    bandl = 0;
else
    bandl = tmp1(end);
end
tmp2 = find(fss>wh);
if (isempty(tmp2))
    bandh = length(fss);
else
    bandh = tmp2(1);
end

% filtering
pasf = zeros(ds(1),ds(2));
for i=1:ds(2)
    tmp = fft(orif(:,i));
    tmp(1:bandl) = 0;  tmp(bandh+1:end) = 0;
    pasf(:,i) = 2.0 * real(ifft(tmp));
end

end