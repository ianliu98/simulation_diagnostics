clc
clear

%% value saved

%%%
% particle 1
%%%

% 1.
% frel = 0.0718;
% freh = 0.0732;
% tpl = 1550;
% tph = 1830;
% h = -42;

% 2.
% frel = 0.0561;
% freh = 0.0596;
% tpl = 18900;
% tph = 19300;

% 2.1
% frel = 0.0562;
% freh = 0.0571;
% tpl = 18930;
% tph = 19230;

% 2.2
% frel = 0.058;
% freh = 0.0588;
% tpl = 18930;
% tph = 19230;

% 2.3
% frel = 0.059;
% freh = 0.06;
% tpl = 18930;
% tph = 19230;

% 3
% frel = 0.0646
% freh = 0.0675
% tpl = 10213
% tph = 10375

% 4
% frel = 0.0594
% freh = 0.0613
% tpl = 13356
% tph = 13420

% 5
% frel = 0.0564
% freh = 0.0569
% tpl = 13912
% tph = 14003
%%

% input fre and time range
frel = 0.0381;
freh = 0.0781;
tpl = 18930;
tph = 19230;

tskip = 4;
xskip = 4;
cv = 100;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');


data_size = size(by);

len = data_size(1);
dts = ddt;
ws = 2*pi/dts;
fft_pnt_f = len;
dw_f = ws / fft_pnt_f;
fss = (1:1:len) * dw_f - dw_f/2.0;

tmp1 = find(fss<frel);
if (isempty(tmp1))
    bandl = 0;
else
    bandl = tmp1(end);
end
tmp2 = find(fss>freh);
if (isempty(tmp2))
    bandh = length(fss);
else
    bandh = tmp2(1);
end


pasby = zeros(data_size(1),data_size(2));
pasbz = zeros(data_size(1),data_size(2));

for i=1:data_size(2)
    tmpy = fft(by(:,i));
    tmpz = fft(bz(:,i));
    tmpy(1:bandl) = 0;  tmpy(bandh+1:end) = 0;
    tmpz(1:bandl) = 0;  tmpz(bandh+1:end) = 0;
    pasby(:,i) = 2.0 * real(ifft(tmpy));
    pasbz(:,i) = 2.0 * real(ifft(tmpz));
end

clear by bz tmpy tmpz tmp1 tmp2


%% j e field
ey = readmatrix('./csv/homo_ey');
ez = readmatrix('./csv/homo_ez');

jy = readmatrix('./csv/homo_jy');
jz = readmatrix('./csv/homo_jz');

pasey = band_pass(ey, frel, freh);
pasez = band_pass(ez, frel, freh);

pasjy = band_pass(jy, frel, freh);
pasjz = band_pass(jz, frel, freh);

ajb = (pasjy .* pasby + pasjz .* pasbz) ./ sqrt(pasby.^2 + pasbz.^2);
aje = (pasjy .* pasey + pasjz .* pasez) ./ sqrt(pasey.^2 + pasez.^2);

clear ey ez jy jz pasey pasez pasjy pasjz

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

%% analysis
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

tpl_ind = find(ytick<tpl);
tpl_ind = tpl_ind(end)+1 - 5;
tph_ind = find(ytick>tph);
tph_ind = tph_ind(1)-1 + 5;

tx_fld = forf(tpl_ind:tph_ind,:);
tx_fre = fre(tpl_ind:tph_ind,:);
tx_inh = S(tpl_ind:tph_ind,:);
tx_dwt = dwdt(tpl_ind:tph_ind,:);
tx_jee = aje(tpl_ind:tph_ind,:);
tx_jbb = ajb(tpl_ind:tph_ind,:);

ytick2 = ytick(tpl_ind:tph_ind);
maxv = max(max(tx_fld));

%clear forf fre S dwdt
%% plot

skip_x = 32;
skip_t = 20;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_fld(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('B_w [B_0]')
axis tight;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_fre(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
caxis([frel freh])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('\omega [\Omega_{e0}]')
axis tight;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_inh(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
caxis([-2 2])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('S')
axis tight;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_jee(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Je')
axis tight;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_jbb(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Jb')
axis tight;

figure,
colormap('jet')
mesh(xtick(1:skip_x:end),ytick2,tx_jbb(:,1:skip_x:end)./tx_fld(:,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([-200 200])
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Jb/Bw')
axis tight;


%% test
clc
skip_x = 32;
skip_t = 20;
rgt = 5960;
lft = 4083;

pnth = [90.08 80.33 73.75 60.41 50.07 39.02 20.21];
pntt = [18933.8 18973.3 19001.3 19057.5 19101.7 19148.8 19228.7];
pntz = [9000.08 8000.33 7300.75 6000.41 5000.07 3900.02 2000.21];
txt = [" A"," B"," C"," D"," E"," F"," G"];

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_fld(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([7.19e-6 6.42e-4])
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('B_w [B_0]')
axis tight;

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_fre(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
colorbar()
caxis([frel freh])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('\omega [\Omega_{e0}]')
axis tight;

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_inh(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
colorbar()
caxis([-2 2])
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('S')
axis tight;

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_jee(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Je')
axis tight;

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_jbb(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Jb')
axis tight;

figure,
colormap('jet')
mesh(xtick(lft:rgt),ytick2,tx_jbb(:,lft:rgt)./tx_fld(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
caxis([-200 200])
colorbar()
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
title('Jb/Bw')
axis tight;

%%
clc
rgt = 5960;
lft = 4083;
pnth = [90.08 80.33 73.75 60.41 50.07 39.02 20.21];
pntt = [18933.8 18973.3 19001.3 19057.5 19101.7 19148.8 19228.7];
pntz = [9000.08 8000.33 7300.75 6000.41 5000.07 3900.02 2000.21];
txt = [" A"," B"," C"," D"," E"," F"," G"];

figure,
colormap(jet)

% wave amplitude
subplot(3,2,1)
mesh(xtick(lft:rgt),ytick2,tx_fld(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([7.19e-6 6.42e-4])
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
ylabel('t [\Omega_{e0}^{-1}')
h=colorbar();
ylabel(h,'B_w [B_0]')
axis tight;
text(-0.15,1.075,'(a)','Units','normalized','FontSize',10)

% frequency
subplot(3,2,2)
mesh(xtick(lft:rgt),ytick2,tx_fre(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
h=colorbar();
caxis([frel freh])
text(-0.15,1.075,'(b)','Units','normalized','FontSize',10)
ylabel(h,'\omega [\Omega_{e0}]')
axis tight;

% S
subplot(3,2,3)
mesh(xtick(lft:rgt),ytick2,tx_inh(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'white*')
text(pnth,pntt,pntz,txt,'Color','white')
h=colorbar();
caxis([-10 10])
ylabel('t [\Omega_{e0}^{-1}')
ylabel(h,'S')
text(-0.15,1.075,'(c)','Units','normalized','FontSize',10)
axis tight;

% jbbw
subplot(3,2,4)
mesh(xtick(lft:rgt),ytick2,tx_jbb(:,lft:rgt)./tx_fld(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
caxis([-200 200])
h=colorbar();
text(-0.15,1.075,'(d)','Units','normalized','FontSize',10)
ylabel(h,'J_B/B_w')
axis tight;

% jb
subplot(3,2,5)
mesh(xtick(lft:rgt),ytick2,tx_jbb(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([-0.02 0.02])
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
h=colorbar();
xlabel('h [c\Omega_{e0}^{-1}]')
ylabel('t [\Omega_{e0}^{-1}')
text(-0.15,1.075,'(e)','Units','normalized','FontSize',10)
ylabel(h,'-J_B [q]')
axis tight;

% je
subplot(3,2,6)
mesh(xtick(lft:rgt),ytick2,-tx_jee(:,lft:rgt),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([-0.02 0.02])
hold on
plot3(pnth,pntt,pntz,'black*')
text(pnth,pntt,pntz,txt)
h=colorbar();
xlabel('h [c\Omega_{e0}^{-1}]')
text(-0.15,1.075,'(f)','Units','normalized','FontSize',10)
ylabel(h,'-J_E [q]')
axis tight;
