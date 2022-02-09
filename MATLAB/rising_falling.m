% check the correspondence between w-t & instantaneous fre
clc
clear

% load data
fyn100 = readmatrix('./csv/homo_forfy_hn100.csv');
fyn50  = readmatrix('./csv/homo_forfy_hn50.csv');
fy0    = readmatrix('./csv/homo_forfy_h0.csv');
fyp50  = readmatrix('./csv/homo_forfy_hp50.csv');
fyp100 = readmatrix('./csv/homo_forfy_hp100.csv');

fzn100 = readmatrix('./csv/homo_forfz_hn100.csv');
fzn50  = readmatrix('./csv/homo_forfz_hn50.csv');
fz0    = readmatrix('./csv/homo_forfz_h0.csv');
fzp50  = readmatrix('./csv/homo_forfz_hp50.csv');
fzp100 = readmatrix('./csv/homo_forfz_hp100.csv');

jbn100 = readmatrix('./csv/homo_jb_hn100.csv');
jbn50  = readmatrix('./csv/homo_jb_hn50.csv');
jb0    = readmatrix('./csv/homo_jb_h0.csv');
jbp50  = readmatrix('./csv/homo_jb_hp50.csv');
jbp100 = readmatrix('./csv/homo_jb_hp100.csv');

jen100 = readmatrix('./csv/homo_je_hn100.csv');
jen50  = readmatrix('./csv/homo_je_hn50.csv');
je0    = readmatrix('./csv/homo_je_h0.csv');
jep50  = readmatrix('./csv/homo_je_hp50.csv');
jep100 = readmatrix('./csv/homo_je_hp100.csv');

% parameters
dt = 0.004;
ifdiag = 256;
ddt = dt * ifdiag;
nsteps = 5242880*2;

tt = ddt:ddt:dt*nsteps;

%% choose a packet from dynamic spectrum
wtf = readmatrix('./csv/homo_wtf.csv');
wtf_t = readmatrix('./csv/homo_wtf_time.csv');
wtf_w = readmatrix('./csv/homo_wtf_fre.csv');

wtfs = reshape(wtf,[length(wtf_t), length(wtf_w), 5]);
wtf_tt = wtf_t * 1e4;

ind = 5;
wtf_maxv = max(max(wtfs(:,:,ind)));

figure,
colormap('jet')
mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
%mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)))
view(2)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\omega [\Omega_{e0}]')
h = colorbar;
ylabel(h,'log_{10}(B_w/B_0)')
axis tight;
caxis([wtf_maxv-1 wtf_maxv])

points = 2;
[tp, wp] = ginput(points);

wl = min(wp);
wh = max(wp);
tl = min(tp);
th = max(tp);

wavel = round(2*pi/((wh+wl)/2)/ddt);
%% bandpass filter
len = length(fy0);
dts = ddt;
ws = 2*pi/dts;
fft_pnt_f = len;
dw_f = ws / fft_pnt_f;
fss = (1:1:len) * dw_f - dw_f/2.0;

wl_ind_vec = find(fss < wl);
wh_ind_vec = find(fss > wh);
wl_ind = wl_ind_vec(end) + 1;
wh_ind = wh_ind_vec(1) -1;

% fft
fft_fyn100 = fft(fyn100);
fft_fyn50  = fft(fyn50);
fft_fy0    = fft(fy0);
fft_fyp50  = fft(fyp50);
fft_fyp100 = fft(fyp100);

fft_fzn100 = fft(fzn100);
fft_fzn50  = fft(fzn50);
fft_fz0    = fft(fz0);
fft_fzp50  = fft(fzp50);
fft_fzp100 = fft(fzp100);

fft_jbn100 = fft(jbn100);
fft_jbn50  = fft(jbn50);
fft_jb0    = fft(jb0);
fft_jbp50  = fft(jbp50);
fft_jbp100 = fft(jbp100);

fft_jen100 = fft(jen100);
fft_jen50  = fft(jen50);
fft_je0    = fft(je0);
fft_jep50  = fft(jep50);
fft_jep100 = fft(jep100);

% elimination
fft_fyn100(1:wl_ind-1) = 0;  fft_fyn100(wh_ind+1:end) = 0;
fft_fyn50(1:wl_ind-1) = 0;   fft_fyn50(wh_ind+1:end) = 0;
fft_fy0(1:wl_ind-1) = 0;     fft_fy0(wh_ind+1:end) = 0;
fft_fyp50(1:wl_ind-1) = 0;   fft_fyp50(wh_ind+1:end) = 0;
fft_fyp100(1:wl_ind-1) = 0;  fft_fyp100(wh_ind+1:end) = 0;

fft_fzn100(1:wl_ind-1) = 0;  fft_fzn100(wh_ind+1:end) = 0;
fft_fzn50(1:wl_ind-1) = 0;   fft_fzn50(wh_ind+1:end) = 0;
fft_fz0(1:wl_ind-1) = 0;     fft_fz0(wh_ind+1:end) = 0;
fft_fzp50(1:wl_ind-1) = 0;   fft_fzp50(wh_ind+1:end) = 0;
fft_fzp100(1:wl_ind-1) = 0;  fft_fzp100(wh_ind+1:end) = 0;

fft_jbn100(1:wl_ind-1) = 0;  fft_jbn100(wh_ind+1:end) = 0;
fft_jbn50(1:wl_ind-1) = 0;   fft_jbn50(wh_ind+1:end) = 0;
fft_jb0(1:wl_ind-1) = 0;     fft_jb0(wh_ind+1:end) = 0;
fft_jbp50(1:wl_ind-1) = 0;   fft_jbp50(wh_ind+1:end) = 0;
fft_jbp100(1:wl_ind-1) = 0;  fft_jbp100(wh_ind+1:end) = 0;

fft_jen100(1:wl_ind-1) = 0;  fft_jen100(wh_ind+1:end) = 0;
fft_jen50(1:wl_ind-1) = 0;   fft_jen50(wh_ind+1:end) = 0;
fft_je0(1:wl_ind-1) = 0;     fft_je0(wh_ind+1:end) = 0;
fft_jep50(1:wl_ind-1) = 0;   fft_jep50(wh_ind+1:end) = 0;
fft_jep100(1:wl_ind-1) = 0;  fft_jep100(wh_ind+1:end) = 0;

% ifft
iff_fyn100 = 2 * real(ifft(fft_fyn100));
iff_fyn50  = 2 * real(ifft(fft_fyn50));
iff_fy0    = 2 * real(ifft(fft_fy0));
iff_fyp50  = 2 * real(ifft(fft_fyp50));
iff_fyp100 = 2 * real(ifft(fft_fyp100));

iff_fzn100 = 2 * real(ifft(fft_fzn100));
iff_fzn50  = 2 * real(ifft(fft_fzn50));
iff_fz0    = 2 * real(ifft(fft_fz0));
iff_fzp50  = 2 * real(ifft(fft_fzp50));
iff_fzp100 = 2 * real(ifft(fft_fzp100));

iff_jbn100 = 2 * real(ifft(fft_jbn100));
iff_jbn50  = 2 * real(ifft(fft_jbn50));
iff_jb0    = 2 * real(ifft(fft_jb0));
iff_jbp50  = 2 * real(ifft(fft_jbp50));
iff_jbp100 = 2 * real(ifft(fft_jbp100));

iff_jen100 = 2 * real(ifft(fft_jen100));
iff_jen50  = 2 * real(ifft(fft_jen50));
iff_je0    = 2 * real(ifft(fft_je0));
iff_jep50  = 2 * real(ifft(fft_jep50));
iff_jep100 = 2 * real(ifft(fft_jep100));

% instantaneous frequency calculation -- Nogi
fre_n100 = Nogi_frequency(iff_fyn100, iff_fzn100, ddt);
fre_n50  = Nogi_frequency(iff_fyn50 , iff_fzn50, ddt);
fre_0    = Nogi_frequency(iff_fy0   , iff_fz0, ddt);
fre_p50  = Nogi_frequency(iff_fyp50 , iff_fzp50, ddt);
fre_p100 = Nogi_frequency(iff_fyp100, iff_fzp100, ddt);

% instantaneous frequency calculation -- Hikishima
hfre_n100 = Hikishima_frequency(iff_fyn100, iff_fzn100, ddt);
hfre_n50  = Hikishima_frequency(iff_fyn50 , iff_fzn50, ddt);
hfre_0    = Hikishima_frequency(iff_fy0   , iff_fz0, ddt);
hfre_p50  = Hikishima_frequency(iff_fyp50 , iff_fzp50, ddt);
hfre_p100 = Hikishima_frequency(iff_fyp100, iff_fzp100, ddt);

% fields
fld_n100 = sqrt(iff_fyn100.^2 + iff_fzn100.^2);
fld_n50  = sqrt(iff_fyn50.^2  + iff_fzn50.^2);
fld_0    = sqrt(iff_fy0.^2    + iff_fz0.^2);
fld_p50  = sqrt(iff_fyp50.^2  + iff_fzp50.^2);
fld_p100 = sqrt(iff_fyp100.^2 + iff_fzp100.^2);

%% plotfigure,
subplot(2,1,1)
plot(tt,log10(abs(fld_p100)))
title('h=100')
xlim([tl th])
subplot(2,1,2)
plot(tt,fre_p100)
xlim([tl th])
ylim([wl wh])
figure,
subplot(2,1,1)
plot(tt,log10(abs(fld_n100)))
title('h=-100')
xlim([1.9e4 2.2e4])
subplot(2,1,2)
plot(tt,fre_n100)
xlim([1.9e4 2.2e4])
ylim([wl wh])

figure,
subplot(2,1,1)
plot(tt,log10(abs(fld_n50)))
title('h=-50')
xlim([1.9e4 2.2e4])
subplot(2,1,2)
plot(tt,fre_n50)
xlim([1.9e4 2.2e4])
ylim([wl wh])

figure,
subplot(2,1,1)
plot(tt,log10(abs(fld_0)))
title('h=0')
xlim([tl th])
subplot(2,1,2)
plot(tt,fre_0)
xlim([tl th])
ylim([wl wh])

figure,
subplot(2,1,1)
plot(tt,log10(abs(fld_p50)))
title('h=50')
xlim([1.7e4 2.0e4])
subplot(2,1,2)
plot(tt,fre_p50)
xlim([1.7e4 2.0e4])
ylim([wl wh])

figure,
subplot(2,1,1)
plot(tt,log10(abs(fld_p100)))
title('h=100')
xlim([tl th])
subplot(2,1,2)
plot(tt,fre_p100)
xlim([tl th])
ylim([wl wh])

%% inhomogeneity factor S
wind_mean = 1;
wind_mean_dwdt = 75;

[s0, s1, s2] = inhomogeneity(fre_n100);
fre_n100_mean = movmean(fre_n100,wind_mean);  % only used for dwdt
dwdt_n100 = four_order_appr(fre_n100_mean,ddt);
dwdt_n100 = [dwdt_n100; 0; 0; 0; 0];
S_n100 = -1.0 * s1 .* dwdt_n100 ./ (s0 .* fre_n100 .* fld_n100);

[s0, s1, s2] = inhomogeneity(fre_n50);
fre_n50_mean = movmean(fre_n50,wind_mean);  % only used for dwdt
dwdt_n50 = four_order_appr(fre_n50_mean,ddt);
dwdt_n50 = [dwdt_n50; 0; 0; 0; 0];
S_n50 = -1.0 * s1 .* dwdt_n50 ./ (s0 .* fre_n50 .* fld_n50);

[s0, s1, s2] = inhomogeneity(fre_0);
fre_0_mean = movmean(fre_0,wind_mean);  % only used for dwdt
dwdt_0 = four_order_appr(fre_0_mean,ddt);
dwdt_0 = [dwdt_0; 0; 0; 0; 0];
dwdt_0 = movmean(dwdt_0, wind_mean_dwdt);
%fld_0 = movmean(fld_0, wind_mean_dwdt);
S_0 = -1.0 * s1 .* dwdt_0 ./ (s0 .* fre_0 .* fld_0);

[s0, s1, s2] = inhomogeneity(fre_p50);
fre_p50_mean = movmean(fre_p50,wind_mean);  % only used for dwdt
dwdt_p50 = four_order_appr(fre_p50_mean,ddt);
dwdt_p50 = [dwdt_p50; 0; 0; 0; 0];
S_p50 = -1.0 * s1 .* dwdt_p50 ./ (s0 .* fre_p50 .* fld_p50);

[s0, s1, s2] = inhomogeneity(fre_p100);
fre_p100_mean = movmean(fre_p100,wind_mean);  % only used for dwdt
dwdt_p100 = four_order_appr(fre_p100_mean,ddt);
dwdt_p100 = [dwdt_p100; 0; 0; 0; 0];
dwdt_p100 = movmean(dwdt_p100, wind_mean_dwdt);
fld_p100 = movmean(fld_p100, wind_mean_dwdt);
S_p100 = -1.0 * s1 .* dwdt_p100 ./ (s0 .* fre_p100 .* fld_p100);

%% Growth rate

% nonlinear
Gn_n100 = Nonlinear_growthrate(fld_n100, fre_n100);
Gn_n50  = Nonlinear_growthrate(fld_n50, fre_n50);
Gn_0    = Nonlinear_growthrate(fld_0, fre_0);
Gn_p50  = Nonlinear_growthrate(fld_p50, fre_p50);
Gn_p100 = Nonlinear_growthrate(fld_p100, fre_p100);

% linear (kupdap)
gl_vec = readmatrix('./csv/linear_gr_kupdap.csv');
gl_ind1 = find(gl_vec(:,2) < wl);
gl_ind2 = find(gl_vec(:,2) > wh);
gl = mean(gl_vec(gl_ind1(end)+1:gl_ind2(1)-1, 3));