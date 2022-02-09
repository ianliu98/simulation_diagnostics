% dkdt process
clc
clear

%% load
ktf_all  = readmatrix('./csv/k_t_field_all.csv');
ktk  = readmatrix('./csv/k_t_k.csv');  % w 0.04 ~ 0.06

%% parameters

jebfft = 0;
jebtype = 2;  % 1 -> all; 2 -> hot; 3 -> cold

kmin = 140;
kmax = 200;

dt = 0.004 * 256;
tt = dt:dt:0.004*5242880;

dk = ktk(2) - ktk(1);
kk = kmin*dk:dk:kmax*dk;

% k calculation
we0 = 1.0;
wpc = 15.0;
wph = 0.3;
wpe = sqrt(wpc^2 + wph^2);
cv = 100;
dr = 1;
dw = 3e-6;
w = (1:1:100000) * dw;
xi2 = w .* (we0 - w) ./ wpe^2;
xi2 = abs(xi2);
xi = sqrt(xi2);
chi2 = 1 ./ (1 + xi2);
chi = sqrt(chi2);

kt = w ./ (cv .* chi .* xi);  % theoretical k
kt = kt * cv;

wkmin_tmp = find(kt<=kk(1));
wkmax_tmp = find(kt>=kk(end));
wkmin = w(wkmin_tmp(end));
wkmax = w(wkmax_tmp(1));

wks = wkmin+(wkmax-wkmin)/length(kk):(wkmax-wkmin)/length(kk):wkmax;

%% choose packet
fld = ktf_all(:,kmin:kmax);

% % figure in thesis
% maxvv = max(max(log10(fld(1:20000,:))));
% figure,
% colormap('jet')
% mesh(kk,tt(1:20000),log10(fld(1:20000,:)))
% view(2)
% xlabel('k [c^{-1}\Omega_{e0}]');
% ylabel('t [\Omega_{e0}^{-1}')
% caxis([maxvv-0.5 maxvv])
% h = colorbar();
% ylabel(h,'log_{10}B_w [B_0]');
% axis tight;

maxvv = max(max(log10(fld(1:20000,:))));
figure,
colormap('jet')
%mesh(kk,tt,fld,'EdgeColor','interp','FaceColor','interp')
% mesh(fld(1:20000,:),'EdgeColor','interp','FaceColor','interp')
mesh(log10(fld(1:20000,:)))
view(2)
colorbar()
caxis([maxvv-0.5 maxvv])
axis tight;

[kp, tp] = ginput(2);

% calculate k,t,amp from interpolation
% point 1
intk1 = kp(1) - floor(kp(1));
intk2 = 1 - intk1;
intt1 = tp(1) - floor(tp(1));
intt2 = 1 - intt1;

kpint1 = kk(floor(kp(1))) * intk2 + kk(floor(kp(1))+1) * intk1;
wpint1 = wks(floor(kp(1))) * intk2 + wks(floor(kp(1))+1) * intk1;
tpint1 = tt(floor(tp(1))) * intt2 + tt(floor(tp(1))+1) * intt1;
fldint1 = fld(floor(tp(1)),floor(kp(1))) * (intk2 * intt2) + ...
    fld(floor(tp(1))+1,floor(kp(1))+1) * (intk1 * intt1) + ...
    fld(floor(tp(1))+1,floor(kp(1))) * (intk2 * intt1) + ...
    fld(floor(tp(1)),floor(kp(1))+1) * (intk1 * intt2);

% point 2
intk1 = kp(2) - floor(kp(2));
intk2 = 1 - intk1;
intt1 = tp(2) - floor(tp(2));
intt2 = 1 - intt1;

kpint2 = kk(floor(kp(2))) * intk2 + kk(floor(kp(2))+1) * intk1;
wpint2 = wks(floor(kp(2))) * intk2 + wks(floor(kp(2))+1) * intk1;
tpint2 = tt(floor(tp(2))) * intt2 + tt(floor(tp(2))+1) * intt1;
fldint2 = fld(floor(tp(2)),floor(kp(2))) * (intk2 * intt2) + ...
    fld(floor(tp(2))+1,floor(kp(2))+1) * (intk1 * intt1) + ...
    fld(floor(tp(2))+1,floor(kp(2))) * (intk2 * intt1) + ...
    fld(floor(tp(2)),floor(kp(2))+1) * (intk1 * intt2);


tpl = min([tpint1,tpint2]);
tph = max([tpint1,tpint2]);
frel = min([wpint1,wpint2]);
freh = max([wpint1,wpint2]);

% release memory
clear chi chi2 kt ktdk_all ktf_all w xi xi2 wkmax_tmp wkmin_tmp fld

%% band-pass
tskip = 4;
xskip = 4;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');
if (jebtype==1)
    je = readmatrix('./csv/homo_je_all.csv');
    jb = readmatrix('./csv/homo_jb_all.csv');    
end
if (jebtype==2)
    je = readmatrix('./csv/homo_je_hot.csv');
    jb = readmatrix('./csv/homo_jb_hot.csv');    
end
if (jebtype==3)
    je = readmatrix('./csv/homo_je_cold.csv');
    jb = readmatrix('./csv/homo_jb_cold.csv');    
end

data_size = size(by);

len = data_size(1);
dts = ddt;
ws = 2*pi/dts;
fft_pnt_f = len;
dw_f = ws / fft_pnt_f;
fss = (1:1:len) * dw_f - dw_f/2.0;

tmp1 = find(fss<frel);
tmp2 = find(fss>freh);
bandl = tmp1(end);
bandh = tmp2(1);

pasby = zeros(data_size(1),data_size(2));
pasbz = zeros(data_size(1),data_size(2));
if (jebfft==1)
    pasje = zeros(data_size(1),data_size(2));
    pasjb = zeros(data_size(1),data_size(2));
end
for i=1:data_size(2)
    tmpy = fft(by(:,i));
    tmpz = fft(bz(:,i));
    tmpy(1:bandl) = 0;  tmpy(bandh+1:end) = 0;
    tmpz(1:bandl) = 0;  tmpz(bandh+1:end) = 0;
    pasby(:,i) = 2.0 * real(ifft(tmpy));
    pasbz(:,i) = 2.0 * real(ifft(tmpz));
    if (jebfft==1)
        tmpe = fft(je(:,i));
        tmpb = fft(jb(:,i));
        tmpe(1:bandl) = 0;  tmpe(bandh+1:end) = 0;
        tmpb(1:bandl) = 0;  tmpb(bandh+1:end) = 0;
        pasje(:,i) = 2.0 * real(ifft(tmpe));
        pasjb(:,i) = 2.0 * real(ifft(tmpb));
    end
end

clear by bz tmpy tmpz tmp1 tmp2
if (jebfft==1)
    clear tmpe tmpb je jb
end

%% frequency calculation
fre = zeros(data_size(1),data_size(2));
for j=1:data_size(2)
    fre(:,j) = Nogi_frequency(pasby(:,j), pasbz(:,j), ddt);
end

forf = sqrt(pasby.^2 + pasbz.^2);
clear pasby pasbz

%% frequency variation
dwdt = zeros(data_size(1),data_size(2));
mean_wind = 50;
for j=1:data_size(2)
    dwdt_tmp = four_order_appr(fre(:,j),ddt);
    dwdt_tmp = [dwdt_tmp;0;0;0;0];
    dwdt(:,j) = movmean(dwdt_tmp,mean_wind);
end

%% inhomogeneity factor
[s0, s1, ~] = inhomogeneity(fre);
S = -1.0 * s1 .* dwdt ./ (s0 .* fre .* forf);

clear s0 s1 s2

%% growth rates
% nonlinear
Gn = Nonlinear_growthrate(forf, fre);

% linear (kupdap)
gl_vec = readmatrix('./csv/linear_gr_kupdap.csv');
gl_ind1 = find(gl_vec(:,2) < wpint1);
gl_ind2 = find(gl_vec(:,2) > wpint2);
gl = mean(gl_vec(gl_ind1(end):gl_ind2(1), 3));

%% analysis
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

tpl_ind = find(ytick<tpl);
tpl_ind = tpl_ind(end)+1;
tph_ind = find(ytick>tph);
tph_ind = tph_ind(1)-1;

tx_fld = forf(tpl_ind:tph_ind,:);
tx_fre = fre(tpl_ind:tph_ind,:);
tx_inh = S(tpl_ind:tph_ind,:);
tx_dwt = dwdt(tpl_ind:tph_ind,:);
tx_gnn = Gn(tpl_ind:tph_ind,:);
tx_jee = je(tpl_ind:tph_ind,:);
tx_jbb = jb(tpl_ind:tph_ind,:);

ytick2 = ytick(tpl_ind:tph_ind);
maxv = max(max(tx_fld));

%clear forf fre S dwdt
%% plot

skip_x = 32;
skip_t = 20;

% field
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),forf(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('fld')
colorbar()
axis tight;

% frequency
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),fre(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('fre')
colorbar()
caxis([frel freh])
axis tight;

% dwdt
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),dwdt(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('dwdt')
colorbar()
caxis([-2e-6 2e-6])
axis tight;

% S
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),S(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('S')
colorbar()
caxis([-5 5])
axis tight;

% je
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),je(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('je')
colorbar()
axis tight;

% jb
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),jb(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
title('jb')
colorbar()
axis tight;

% jbbw
coe = -0.5 * 0.028 * 1.257e-6;
figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),jb(1:skip_t:end,1:skip_x:end)./...
    forf(1:skip_t:end,1:skip_x:end).*coe,'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([-2e-6 2e-6])
title('jbbw')
colorbar()
axis tight;

%% plot in thesis
skip_x = 32;
skip_t = 20;

% text(-0.15,1.075,'(a)','Units','normalized','FontSize',12)
edge_pntz = [1e10 1e10];

% edge_pnt1x = [92.2 142.12];
% edge_pnt1y = [10461.2 13017.1];

edge_pnt1x = [56.36 90.16];
edge_pnt1y = [13144.1 15114.2];

edge_pnt2x = [-74.2 -53.72];
edge_pnt2y = [13144.1 15114.2];


figure,
colormap(jet)

% wave amplitude
pmax = max(max(log10(tx_fld)));
subplot(3,2,1)
mesh(xtick(1:skip_x:end),ytick2,log10(tx_fld(:,1:skip_x:end)))
view(2)
h = colorbar();
caxis([pmax-0.5 pmax])
ylabel('t [\Omega_{e0}^{-1}]')
ylabel(h,'log_{10}(|B_w/B_0|)')
text(-0.15,1.075,'(a)','Units','normalized','FontSize',10)
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz, 'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off

% frequency
subplot(3,2,2)
mesh(xtick(1:skip_x:end),ytick2,tx_fre(:,1:skip_x:end))
view(2)
h = colorbar();
caxis([frel freh])
ylabel(h,'\omega [\Omega_{e0}]')
text(-0.15,1.075,'(b)','Units','normalized','FontSize',10)
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz, 'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off

% dwdt
subplot(3,2,3)
mesh(xtick(1:skip_x:end),ytick2,tx_dwt(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
h = colorbar();
caxis([-2e-6 2e-6])
ylabel('t [\Omega_{e0}^{-1}]')
ylabel(h,'\partial\omega/\partial t')
text(-0.15,1.075,'(c)','Units','normalized','FontSize',10)
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz, 'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off

% S
subplot(3,2,4)
mesh(xtick(1:skip_x:end),ytick2,tx_inh(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
h = colorbar();
caxis([-2 2])
ylabel(h,'S')
text(-0.15,1.075,'(d)','Units','normalized','FontSize',10)
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz, 'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off

% je
subplot(3,2,6)
mesh(xtick(1:skip_x:end),ytick2,tx_jbb(:,1:skip_x:end)./tx_fld(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
h = colorbar();
caxis([-100 100])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel(h,'J_B/B_w')
text(-0.15,1.075,'(f)','Units','normalized','FontSize',10)
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz, 'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off

% jb
% subplot(3,2,6)
% mesh(xtick(1:skip_x:end),ytick2,tx_jbb(:,1:skip_x:end)./tx_fld(:,1:skip_x:end));
% view(2)
% h = colorbar;
% caxis([-0.1 0.1])
% xlabel('h [c\Omega_{e0}^{-1}]')
% text(-0.15,1.075,'(f)','Units','normalized','FontSize',10)
% ylabel(h,'J_B')
% axis tight;
% hold on
% plot3(edge_pnt1x, edge_pnt1y, edge_pntz,'black--','LineWidth',2);
% plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
% % plot3(edge_pnt3x, edge_pnt3y, edge_pntz, 'black--','LineWidth',2)
% % plot3(edge_pnt4x, edge_pnt4y, edge_pntz, 'black--','LineWidth',2)
% hold off

% gn
gnm = max(max(log10(abs(tx_gnn(:,1:skip_x:end)))));
subplot(3,2,5)
mesh(xtick(1:skip_x:end),ytick2,log10(abs(tx_gnn(:,1:skip_x:end))));
view(2)
h = colorbar;
%caxis([gnm-1 gnm])
caxis([-3.5 -2.5])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}]')
text(-0.15,1.075,'(e)','Units','normalized','FontSize',10)
ylabel(h,'log_{10}(|\Gamma_N|)')
axis tight;
hold on
plot3(edge_pnt1x, edge_pnt1y, edge_pntz,'black--','LineWidth',2)
plot3(edge_pnt2x, edge_pnt2y, edge_pntz, 'black--','LineWidth',2)
hold off


%% choose a position
% h =125;

edge_h = 75;
edget = 1.41e4;
%edget = 1.2e4;

% edge_h = 125;
% edget = 1.2e4;

mw_h = 75;

hind = (16384 + edge_h * 100 - 2048) / 4;

hfld = tx_fld(:,hind);
hfre = tx_fre(:,hind);
hdwt = tx_dwt(:,hind);
hinh = tx_inh(:,hind);
hjee = tx_jee(:,hind);
hjbb = tx_jbb(:,hind);
hgnn = tx_gnn(:,hind);

hfld = movmean(hfld,mw_h);
hfre = movmean(hfre,mw_h);
hdwt = movmean(hdwt,mw_h);
hinh = movmean(hinh,mw_h);
hjee = movmean(hjee,mw_h);
hjbb = movmean(hjbb,mw_h);
hgnn = movmean(hgnn,mw_h);


figure,
subplot(5,1,1)
plot(ytick2,log10(abs(hfld)),'black')
ylabel('log_{10}(|B_w/B_0|)')
xline(edget,'b--')
text(edget,pmax-0.15,'\leftarrow approximate end of the wave packet','FontSize',10,'Color','blue')
xlim([tpl tph])
text(-0.1,1.075,'(a)','Units','normalized','FontSize',10)
grid on

subplot(5,1,2)
plot(ytick2,hfre,'black')
ylabel('\omega [\Omega_{e0}]')
xlim([tpl tph])
xline(edget,'b--')
grid on
text(-0.1,1.075,'(b)','Units','normalized','FontSize',10)

subplot(5,1,3)
plot(ytick2,hdwt,'black')
% yline(0,'r--')
xlim([tpl tph])
ylabel('\partial\omega/\partial t')
xline(edget,'b--')
grid on
text(-0.1,1.075,'(c)','Units','normalized','FontSize',10)

subplot(5,1,4)
plot(ytick2,hinh,'black')
ylim([-2 2])
xlim([tpl tph])
% yline(1,'r--')
yline(-0.4,'r--')
% text(1.075e4, 0.6, 'S = -0.4','FontSize',10,'Color','red')
text(1.33e4, -0.6, 'S = -0.4','FontSize',10,'Color','red')
% yline(0,'r--')
grid on
ylabel('S')
xline(edget,'b--')
text(-0.1,1.075,'(d)','Units','normalized','FontSize',10)

subplot(5,1,5)
plot(ytick2,log10(abs(hgnn)),'black')
yline(log10(abs(gl)),'r--')
% text(1.075e4, -3.8, 'linear growth rate','FontSize',10,'Color','red')
text(1.33e4, -3.8, 'linear growth rate','FontSize',10,'Color','red')
xlim([tpl tph])
ylim([log10(abs(gl))-0.1 max(log10(abs(hgnn)))])
ylabel('log_{10}(|\Gamma_N|)')
xline(edget,'b--')
xlabel('t [\Omega_{e0}]^{-1}')
grid on
text(-0.1,1.075,'(e)','Units','normalized','FontSize',10)



% subplot(6,1,5)
% plot(ytick2,log10(abs(hjee)))
% ylabel('log_10(|J_E|)')
% 
% subplot(6,1,6)
% plot(ytick2,log10(abs(hjbb)))
% ylabel('log10(|J_B|)')
% xlabel('t [\Omega_{e0}^{-1}]')