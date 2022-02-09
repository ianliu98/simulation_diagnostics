% k-t ~ fre ~ w-t
clc
clear

%% load
ktf_all  = readmatrix('./csv/k_t_field_all.csv');
ktdk_all = readmatrix('./csv/k_t_deltak_all.csv');
ktk  = readmatrix('./csv/k_t_k.csv');  % w 0.04 ~ 0.06

%% parameters
kmin = 50;
kmax = 400;

dt = 0.004 * 256;
tt = dt:dt:0.004*5242880;

dk = ktk(2) - ktk(1);
kk = kmin*dk:dk:kmax*dk;

%% k calculation
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

% figure,
% plot(kt, w)
% xlabel('k')
% ylabel('w')

wmin = 0.04;
wmax = 0.06;

kwmin = kt(floor(wmin/dw));
kwmax = kt(floor(wmax/dw));

wkmin_tmp = find(kt<=kk(1));
wkmax_tmp = find(kt>=kk(end));
wkmin = w(wkmin_tmp(end));
wkmax = w(wkmax_tmp(1));

wks = wkmin+(wkmax-wkmin)/length(kk):(wkmax-wkmin)/length(kk):wkmax;
%% k-t analysis
fld = ktf_all(:,kmin:kmax);
dlk = ktdk_all(:,kmin:kmax);

maxv = max(max(fld));

figure,
colormap('jet')
mesh(kk,tt,fld,'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
axis tight;

%% fre analysis
tskip = 4;
xskip = 4;
ifdiag = 256;
dt = 0.004;
cv = 100;
% load data
file_date = '0104';
fre_file = join([file_date, '_homo_fre.csv']);
fre_data = readmatrix(['./csv/',fre_file]);
fld_file = join([file_date, '_homo_field.csv']);
fld_data = readmatrix(['./csv/',fld_file]);
% data parameters
data_size = size(fre_data);
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

% packet 1: t=9936.9~13638.7  w=0.0383~0.0407 amp=5.831e-5~1.36e-4
pt1 = 9936.9;
pt2 = 13638.7;
pt1_tmp = find(ytick>pt1);
pt2_tmp = find(ytick<pt2);
pck1 = fre_data(pt1_tmp(1):pt2_tmp(end),:);
pck1_fld = 10.^fld_data(pt1_tmp(1):pt2_tmp(end),:);

pw1 = 0.038;
pw2 = 0.041;
pf1 = 5.831e-5;
pf2 = 1.36e-4;
figure,
colormap('jet')
subplot(1,2,1)
mesh(xtick,ytick(pt1_tmp(1):pt2_tmp(end)),pck1)
view(2)
colorbar()
caxis([pw1 pw2])
axis tight;
subplot(1,2,2)
mesh(xtick,ytick(pt1_tmp(1):pt2_tmp(end)),pck1_fld)
view(2)
colorbar()
caxis([pf1 pf2])
axis tight;


% packet 2: t=7348.2~14174.2  w=0.0412~0.0441 amp=7.9246e-5~1.49e-4
pt3 = 7348.2;
pt4 = 14174.2;
pt3_tmp = find(ytick>pt3);
pt4_tmp = find(ytick<pt4);
pck2 = fre_data(pt3_tmp(1):pt4_tmp(end),:);
pck2_fld = fld_data(pt3_tmp(1):pt4_tmp(end),:);

pw3 = 0.041;
pw4 = 0.045;
pf3 = log10(7.9246e-5);
pf4 = log10(1.49e-4);
figure,
colormap('hot')
subplot(1,2,1)
mesh(xtick,ytick(pt3_tmp(1):pt4_tmp(end)),pck2)
view(2)
colorbar()
caxis([pw3 pw4])
axis tight;
subplot(1,2,2)
mesh(xtick,ytick(pt3_tmp(1):pt4_tmp(end)),pck2_fld)
view(2)
colorbar()
caxis([pf3 pf4])
axis tight;