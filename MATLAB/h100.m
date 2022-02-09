% at h=100
clc
clear

% load data
by = readmatrix('./csv/by_h100.csv');
bz = readmatrix('./csv/bz_h100.csv');
jy = readmatrix('./csv/jy_h100.csv');
jz = readmatrix('./csv/jz_h100.csv');
ey = readmatrix('./csv/ey_h100.csv');
ez = readmatrix('./csv/ez_h100.csv');

% parameters
dt = 0.004;
ifdiag = 256;
tskip = 1;
ddt = dt * ifdiag * tskip;
nsteps = 5242880;
tt = ddt:ddt:dt*nsteps;

% % configuration
% wl = 0.048;
% wh = 0.052;
% mw = 50;
% wavel = round(2*pi/((wh+wl)/2)/ddt);
% wavel = round(wavel / 1);

%% choose a packet from dynamic spectrum
wtf = readmatrix('./csv/homo_wtf6.csv');
wtf_t = readmatrix('./csv/homo_wtf6_time.csv');
wtf_w = readmatrix('./csv/homo_wtf6_fre.csv');

wtfs = reshape(wtf,[length(wtf_t), length(wtf_w), 5]);
wtf_tt = wtf_t * 1e4;

ind = 1;
wtf_maxv = max(max(wtfs(:,:,ind)));

figure,
colormap('jet')
mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
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

% wl = 0.04;
% wh = 0.06;

wavel = round(2*pi/((wh+wl)/2)/ddt);
%%
% band pass
pasby = band_pass2(by, wl, wh);
pasbz = band_pass2(bz, wl, wh);
pasjy = band_pass2(jy, wl, wh);
pasjz = band_pass2(jz, wl, wh);
pasey = band_pass2(ey, wl, wh);
pasez = band_pass2(ez, wl, wh);

% compute
fld = sqrt(pasby.^2 + pasbz.^2);
fre = Nogi_frequency(pasby, pasbz, ddt);
flde = sqrt(pasey.^2 + pasez.^2);
fldeb = sqrt(fld.^2 + flde.^2);

jb = (pasjy .* pasby + pasjz .* pasbz) ./ sqrt(pasby.^2 + pasbz.^2);
je = (pasjy .* pasey + pasjz .* pasez) ./ sqrt(pasey.^2 + pasez.^2);

jb = jb ./ ddt;

dwdt = four_order_appr(fre,ddt);
dwdt = [dwdt;0;0;0;0];
dwdt = movmean(dwdt,wavel);

[s0, s1, ~] = inhomogeneity(fre);
inh = -1.0 * s1 .* dwdt ./ (s0 .* fre .* fld);
inh2 = -1.0 * s1 .* dwdt ./ (s0 .* fre .* fldeb);

gn = Nonlinear_growthrate(fld, fre);
gn2 = Nonlinear_growthrate(fldeb, fre);

jbbw = jb ./ fld;

fld  = smooth(fld,wavel);
fre  = smooth(fre,wavel);
dwdt = smooth(dwdt,wavel);
jb   = smooth(jb,wavel);
je   = smooth(je,wavel);
inh  = smooth(inh,wavel);
gn   = smooth(gn,wavel);
jbbw = smooth(jbbw,wavel);
%% plot
figure,

subplot(7,1,1)
plot(tt,log10(fld))
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('log_{10}(B_w/B_0)')
axis tight
xlim([tt(1) tt(end)])
xlim([tl th])

subplot(7,1,2)
plot(tt,fre)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\omega [\Omega_{e0}]')
axis tight
ylim([wl-0.001 wh+0.001])
xlim([tt(1) tt(end)])
xlim([tl th])

subplot(7,1,3)
plot(tt,dwdt)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\partial\omega/\partial t')
axis tight
yline(0,'r--')
xlim([tt(1) tt(end)])
xlim([tl th])

subplot(7,1,4)
plot(tt,inh)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('S')
axis tight
ylim([-5 5])
yline(0,'r--')
xlim([tt(1) tt(end)])
xlim([tl th])

subplot(7,1,5)
plot(tt,je)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('J_E')
axis tight
yline(0,'r--')
xlim([tt(1) tt(end)])
xlim([tl th])

subplot(7,1,6)
plot(tt,jb)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('J_B')
axis tight
yline(0,'r--')
xlim([tt(1) tt(end)])
xlim([tl th])

% subplot(7,1,7)
% plot(tt,jbbw)
% xlabel('t [\Omega_{e0}^{-1}]')
% ylabel('J_B/B_w')
% axis tight
% xlim([tt(1) tt(end)])
% ylim([-100 100])
% yline(0,'r--')
% xlim([tl th])

subplot(7,1,7)
plot(tt,delta_w)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\Delta\omega')
axis tight
xlim([tt(1) tt(end)])
xlim([tl th])

% ----
% subplot(7,1,5)
% plot(tt,inh2)
% xlabel('t [\Omega_{e0}^{-1}]')
% ylabel('inh2')
% axis tight
% yline(0,'r--')
% ylim([-5 5])
% xlim([tt(1) tt(end)])
% xlim([tl th])
% 
% subplot(7,1,6)
% plot(tt,gn)
% xlabel('t [\Omega_{e0}^{-1}]')
% ylabel('gn')
% axis tight
% yline(0,'r--')
% xlim([tt(1) tt(end)])
% xlim([tl th])
% 
% subplot(7,1,7)
% plot(tt,gn2)
% xlabel('t [\Omega_{e0}^{-1}]')
% ylabel('gn2')
% axis tight
% xlim([tt(1) tt(end)])
% yline(0,'--')
% xlim([tl th])
%% functions

function nw = smooth(old, mw)
    nw = movmean(old,mw);
end