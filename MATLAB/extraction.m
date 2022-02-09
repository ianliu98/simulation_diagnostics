clc
clear

% --- main extraction script --- %

% parameters
frel = 0.04;
freh = 0.06;

tskip = 4;
xskip = 4;
ifdiag = 256;
dt = 0.004;
cv = 100;
ddt = dt * tskip * ifdiag;

% configuration
pkt       =     1;      % which packet
ld        =     1;      % load or draw
offn      =     0;      % size of error windown
mean_wind =     50;     % movmean window size


if (pkt==1)
    error_tw = 425;
    pkti = "1";
else
    error_tw = 350;
    pkti = "2";
end

erro_size = round(error_tw/ddt);

%% b field
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');

pasby = band_pass(by, frel, freh);
pasbz = band_pass(bz, frel, freh);

clear by bz

%% data parameters
data_size = size(pasby);
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

%% frequency calculation
afre = zeros(data_size(1),data_size(2));
for j=1:data_size(2)
    afre(:,j) = Nogi_frequency(pasby(:,j), pasbz(:,j), ddt);
end

aforf = sqrt(pasby.^2 + pasbz.^2);


%% extraction -- single

if (ld == 1)
    load("./mat/extraction_packet_"+pkti+".mat")
else
    points = 5;  % points to choose

    maxv = max(max(log10(abs(aforf))));
    figure,
    colormap(jet)
    mesh(log10(abs(aforf)))
    view(2)
    caxis([maxv-1, maxv])
    colorbar()
    axis tight;

    hold on
    [x,y] = ginput(points);

    line = [];
    for i=2:points
        k = (y(i) - y(i-1)) / (x(i) - x(i-1));
        x_left  = floor(x(i-1));
        x_right = floor(x(i));
        if (i == 2)
            x_left = 1;
        end
        if (i == points)
            x_right = data_size(2);
        end
        for j=x_left:x_right
            edge = k * j + y(i) - k * x(i);
            line = [line, edge];
        end
    end

    bot = max(line(1),0);
    top = min(line(data_size(2)),data_size(1));
    lft = max(0, floor(x(1)-y(1)/k));
    rgt = min(data_size(2), floor((data_size(1)-y(1))/k+x(1)));

    clear ponits x y x_left x_right edge i j k
end

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

clear ey ez jy jz pasey pasez pasjy pasjz %pasby pasbz

%% frequency variation
adwdt = zeros(data_size(1),data_size(2));

for j=1:data_size(2)
    dwdt_tmp = four_order_appr(afre(:,j),ddt);
    dwdt_tmp = [dwdt_tmp;0;0;0;0];
    adwdt(:,j) = movmean(dwdt_tmp,mean_wind);
end

clear dwdt_tmp

%% inhomogeneity factor
[s0, s1, ~] = inhomogeneity(afre);
aS = -1.0 * s1 .* adwdt ./ (s0 .* afre .* aforf);

clear s0 s1 s2

%% nonlinear growth rate
agn = Nonlinear_growthrate(aforf, afre);

%% linear growth rate (kupdap)
gl_vec = readmatrix('./csv/linear_gr_kupdap.csv');
gl_ind1 = find(gl_vec(:,2) < frel);
gl_ind2 = find(gl_vec(:,2) > freh);
agl = mean(gl_vec(gl_ind1(end):gl_ind2(1), 3));

clear gl_vec gl_ind1 gl_ind2

%% convective growth rate
dh = 1.0 * xskip / cv;
pdbdt = zeros(data_size(1),data_size(2));
pdbdh = zeros(data_size(1),data_size(2));
% aafld = movmean(aforf,mean_wind,1);
% aafld = movmean(aafld,mean_wind,2);
aafld = aforf;
for j=1:data_size(2)
    dbdt_tmp = four_order_appr(aafld(:,j),ddt);
    dbdt_tmp = [dbdt_tmp;0;0;0;0];
    pdbdt(:,j) = movmean(dbdt_tmp,mean_wind/2);
end
for j=1:data_size(1)
    dbdh_tmp = four_order_appr(aafld(j,:),dh);
    dbdh_tmp = [dbdh_tmp;0;0;0;0];
    pdbdh(j,:) = movmean(dbdh_tmp,mean_wind/2);
end
vg = Group_velocity(afre)./cv;
dbdt = pdbdt + vg .* pdbdh;
agc = (log(abs(1.0 + dbdt./aforf)) ./ dh) .* vg;

clear pdbdt pdbdh dbdt_tmp dbdh_tmp vg dbdt aafre

%% extracting

vg = ((rgt - lft) * xskip / cv) / ((top - bot) * ddt);

if (offn==0)
    off = 0;
else
    off = floor(linspace(erro_size/offn,erro_size,offn));
end
bfld = zeros(2*offn+1, data_size(2));
bjee = zeros(2*offn+1, data_size(2));
bfre = zeros(2*offn+1, data_size(2));
bdwt = zeros(2*offn+1, data_size(2));
binh = zeros(2*offn+1, data_size(2));
bgcc = zeros(2*offn+1, data_size(2));
bgnn = zeros(2*offn+1, data_size(2));

ew = floor(erro_size/2);

ccfld = zeros(2*ew+1, data_size(2));
ccfre = zeros(2*ew+1, data_size(2));
ccjee = zeros(2*ew+1, data_size(2));
ccjbb = zeros(2*ew+1, data_size(2));
ccdwt = zeros(2*ew+1, data_size(2));
ccinh = zeros(2*ew+1, data_size(2));
ccgcc = zeros(2*ew+1, data_size(2));
ccgnn = zeros(2*ew+1, data_size(2));
for i=1:data_size(2)
    yy = floor(line(i));
    if (off==0)
        yys = yy;
    else
        yys = [yy-off,yy,yy+off];
    end
    % field
    bfld(yys>0,i) = aforf(yys,i);  bfld(yys<=0,i) = 1e-20;
    mnfld(i) = mean(aforf(yy-ew:yy+ew,i));
    upfld(i) = aforf(yy+ew,i);
    dwfld(i) = aforf(yy-ew,i);
    mnupfld(i) = mean(aforf(yy:yy+ew,i));
    mndwfld(i) = mean(aforf(yy-ew:yy,i));
    ccfld(:,i) = aforf(yy-ew:yy+ew,i);
    % freq
    bfre(yys>0,i) = afre(yys,i);   bfre(yys<=0,i) = 1e-20;
    mnfre(i) = mean(afre(yy-ew:yy+ew,i));
    upfre(i) = afre(yy+ew,i);
    dwfre(i) = afre(yy-ew,i);
    mnupfre(i) = mean(afre(yy:yy+ew,i));
    mndwfre(i) = mean(afre(yy-ew:yy,i));
    ccfre(:,i) = afre(yy-ew:yy+ew,i);
    % je
    bjee(yys>0,i) = aje(yys,i);    bjee(yys<=0,i) = 1e-20;
    mnje(i) = mean(aje(yy-ew:yy+ew,i));
    upje(i) = aje(yy+ew,i);
    dwje(i) = aje(yy-ew,i);
    mnupje(i) = mean(aje(yy:yy+ew,i));
    mndwje(i) = mean(aje(yy-ew:yy,i));
    ccjee(:,i) = aje(yy-ew:yy+ew,i);
    % jb
    bjbb(yys>0,i) = ajb(yys,i);    bjbb(yys<=0,i) = 1e-20;
    mnjb(i) = mean(ajb(yy-ew:yy+ew,i));
    upjb(i) = ajb(yy+ew,i);
    dwjb(i) = ajb(yy-ew,i);
    mnupjb(i) = mean(ajb(yy:yy+ew,i));
    mndwjb(i) = mean(ajb(yy-ew:yy,i));
    ccjbb(:,i) = ajb(yy-ew:yy+ew,i);
    % dwdt
    bdwt(yys>0,i) = adwdt(yys,i);  bdwt(yys<=0,i) = 1e-20;
    mndwt(i) = mean(adwdt(yy-ew:yy+ew,i));
    updwt(i) = adwdt(yy+ew,i);
    dwdwt(i) = adwdt(yy-ew,i);
    mnupdwt(i) = mean(adwdt(yy:yy+ew,i));
    mndwdwt(i) = mean(adwdt(yy-ew:yy,i));
    ccdwt(:,i) = adwdt(yy-ew:yy+ew,i);
    % S
    binh(yys>0,i) = aS(yys,i);     binh(yys<=0,i) = 1e-20;
    mninh(i) = mean(aS(yy-ew:yy+ew,i));
    upinh(i) = aS(yy+ew,i);
    dwinh(i) = aS(yy-ew,i);
    mnupinh(i) = mean(aS(yy:yy+ew,i));
    mndwinh(i) = mean(aS(yy-ew:yy,i));
    ccinh(:,i) = aS(yy-ew:yy+ew,i);
    % convective gr
    bgcc(yys>0,i) = agc(yys,i);    bgcc(yys<=0,i) = 1e-20;
    mngcc(i) = mean(agc(yy-ew:yy+ew,i));
    upgcc(i) = agc(yy+ew,i);
    dwgcc(i) = agc(yy-ew,i);
    mnupgcc(i) = mean(agc(yy:yy+ew,i));
    mndwgcc(i) = mean(agc(yy-ew:yy,i));
    ccgcc(:,i) = agc(yy-ew:yy+ew,i);
    % nonlinear gr
    bgnn(yys>0,i) = agn(yys,i);    bgnn(yys<=0,i) = 1e-20;
    mngnn(i) = mean(agn(yy-ew:yy+ew,i));
    upgnn(i) = agn(yy+ew,i);
    dwgnn(i) = agn(yy-ew,i);
    mnupgnn(i) = mean(agn(yy:yy+ew,i));
    mndwgnn(i) = mean(agn(yy-ew:yy,i));
    ccgnn(:,i) = agn(yy-ew:yy+ew,i);
end

% tick transfer -> along space
htick = xtick(max(floor(lft),1):floor(rgt));
ttick = bot:(top-bot)/(rgt-lft):top-(top-bot)/(rgt-lft);
ttick = ytick(floor(ttick));

clear yy yys

%% plot 1

% fld
figure,
subplot(5,1,1)
plot(htick,log10(abs(mnfld)),'black')
hold on
plot(htick,log10(abs(mndwfld)),'red')
plot(htick,log10(abs(mnupfld)),'blue')
ylabel('log_{10}(|B_w/B_0|)')
text(-0.1,1.075,'(a)','Units','normalized','FontSize',10)
%ylim([max(log10(abs(mnfld)))-1 max(log10(abs(mnfld)))])
axis tight

subplot(5,1,2)
plot(htick,mnfre,'black')
hold on
plot(htick,mndwfre,'red')
plot(htick,mnupfre,'blue')
yline(0.05,'g--')
ylabel('\omega [\Omega_{e0}]')
text(-0.1,1.075,'(b)','Units','normalized','FontSize',10)
ylim([0.04 0.06])

subplot(5,1,3)
plot(htick,mndwt,'black')
hold on
plot(htick,mndwdwt,'red')
plot(htick,mnupdwt,'blue')
yline(0,'g--')
ylabel('\partial\omega/\partial t')
text(-0.1,1.075,'(c)','Units','normalized','FontSize',10)
ylim([-2e-5 2e-5])

subplot(5,1,4)
plot(htick,mninh,'black')
hold on
plot(htick,mndwinh,'red')
plot(htick,mnupinh,'blue')
ylim([-10 10])
yline(0,'g--')
ylabel('S')
text(-0.1,1.075,'(d)','Units','normalized','FontSize',10)

% subplot(5,1,5)
% plot(htick,log10(abs(mngcc)),'black')
% hold on
% plot(htick,log10(abs(mndwgcc)),'red')
% plot(htick,log10(abs(mnupgcc)),'blue')
% plot(htick,log10(abs(mngnn)),'black--')
% plot(htick,log10(abs(mndwgnn)),'red--')
% plot(htick,log10(abs(mnupgnn)),'blue--')
% yline(log10(abs(agl)),'g--')
% ylabel('\Gamma')
% xlabel('h [c\Omega_{e0}]^{-1}')
% text(-0.1,1.075,'(e)','Units','normalized','FontSize',10)

subplot(5,1,5)
plot(htick,log10(abs(mngcc)),'black')
hold on
plot(htick,log10(abs(mngnn)),'red')
yline(log10(abs(agl)),'g--')
ylabel('\Gamma')
xlabel('h [c\Omega_{e0}]^{-1}')
text(-0.1,1.075,'(e)','Units','normalized','FontSize',10)
ylim([-5 -1])
%% plot 2

% fld
figure,
subplot(5,1,1)
plot(htick,log10(abs(bfld)),'black')
hold on
plot(htick,log10(abs(dwfld)),'red')
plot(htick,log10(abs(upfld)),'blue')
ylabel('log_{10}(|B_w/B_0|)')
text(-0.1,1.075,'(a)','Units','normalized','FontSize',10)
%ylim([max(log10(abs(mnfld)))-1 max(log10(abs(mnfld)))])
axis tight

subplot(5,1,2)
plot(htick,bfre,'black')
hold on
plot(htick,dwfre,'red')
plot(htick,upfre,'blue')
yline(0.05,'g--')
ylabel('\omega [\Omega_{e0}]')
text(-0.1,1.075,'(b)','Units','normalized','FontSize',10)
ylim([0.04 0.06])

subplot(5,1,3)
plot(htick,bdwt,'black')
hold on
plot(htick,dwdwt,'red')
plot(htick,updwt,'blue')
yline(0,'g--')
ylabel('\partial\omega/\partial t')
text(-0.1,1.075,'(c)','Units','normalized','FontSize',10)
ylim([-5e-5 5e-5])

subplot(5,1,4)
plot(htick,binh,'black')
hold on
plot(htick,dwinh,'red')
plot(htick,upinh,'blue')
ylim([-10 10])
yline(0,'g--')
ylabel('S')
text(-0.1,1.075,'(d)','Units','normalized','FontSize',10)

subplot(5,1,5)
plot(htick,log10(abs(mngcc)),'black')
hold on
plot(htick,log10(abs(mngnn)),'red')
yline(log10(abs(agl)),'g--')
ylabel('\Gamma')
xlabel('h [c\Omega_{e0}]^{-1}')
text(-0.1,1.075,'(e)','Units','normalized','FontSize',10)
ylim([-5 -1])


%% plot 3
clc

etick = linspace(-error_tw,error_tw,2*ew+1);
htmp = length(htick);
htick_pos = [htick(round(htmp/6))  htick(round(htmp/2)) ...
     htick(round(htmp/6*5))];
htick_row1 = {'-100', '0', '100'};
htick_row2 = [ttick(round(htmp/6)),ttick(round(htmp/6*3)), ...
    ttick(round(htmp/6*5))];
labelArray = [htick_row1; compose('(%.0f)',htick_row2)]; 
htick_labels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

figure,
colormap('jet')

mccfld = max(max(log10(abs(ccfld))));
subplot(3,2,1)
mesh(htick,etick,log10(abs(ccfld)),'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([mccfld-1 mccfld])
ylabel('W_e [\Omega_{e0}^{-1}]')
h = colorbar();
ylabel(h,'log_{10}(|B_w/B_0|)')
xticks(htick_pos)
xticklabels(htick_labels)
axis tight
text(-0.1,1.1,'(a)','Units','normalized','FontSize',10)

subplot(3,2,2)
mesh(htick,etick,ccfre,'EdgeColor','interp','FaceColor','interp')
view(2)
%caxis([0.045 0.055])
xticks(htick_pos)
xticklabels(htick_labels)
h = colorbar();
ylabel(h,'\omega [\Omega_{e0}]')
axis tight
text(-0.1,1.1,'(b)','Units','normalized','FontSize',10)

subplot(3,2,3)
mesh(htick,etick,ccdwt,'EdgeColor','interp','FaceColor','interp')
view(2)
caxis([-5e-5 5e-5])
ylabel('W_e [\Omega_{e0}^{-1}]')
xticks(htick_pos)
xticklabels(htick_labels)
h = colorbar();
ylabel(h,'\partial\omega/\partial t')
axis tight
text(-0.1,1.1,'(c)','Units','normalized','FontSize',10)

subplot(3,2,4)
mesh(htick,etick,ccinh,'EdgeColor','interp','FaceColor','interp')
view(2)
xticks(htick_pos)
xticklabels(htick_labels)
h = colorbar();
ylabel(h,'S')
caxis([-10 10])
axis tight
text(-0.1,1.1,'(d)','Units','normalized','FontSize',10)


% subplot(3,2,5)
% mesh(htick,etick,ccjbb./ccfld,'EdgeColor','interp','FaceColor','interp')
% view(2)
% ylabel('W_e [\Omega_{e0}^{-1}]')
% xticks(htick_pos)
% xticklabels(htick_labels)
% text(-0.2,-0.075,'h [c\Omega_{e0}^{-1}]','Units','normalized','FontSize',9)
% text(-0.2,-0.17,'t [\Omega_{e0}^{-1}]','Units','normalized','FontSize',9)
% h = colorbar();
% caxis([-200 200])
% ylabel(h,'log_{10}(|\Gamma_N|)')
% axis tight
% text(-0.1,1.1,'(e)','Units','normalized','FontSize',10)
% 
% 
% subplot(3,2,6)
% mesh(htick,etick,-ccjee,'EdgeColor','interp','FaceColor','interp')
% view(2)
% xticks(htick_pos)
% xticklabels(htick_labels)
% h = colorbar();
% %caxis([-200 200])
% ylabel(h,'log_{10}(|\Gamma_C|)')
% axis tight
% text(-0.1,1.1,'(f)','Units','normalized','FontSize',10)


subplot(3,2,5)
mesh(htick,etick,log10(abs(ccgnn)),'EdgeColor','interp','FaceColor','interp')
view(2)
ylabel('W_e [\Omega_{e0}^{-1}]')
xticks(htick_pos)
xticklabels(htick_labels)
text(-0.2,-0.075,'h [c\Omega_{e0}^{-1}]','Units','normalized','FontSize',9)
text(-0.2,-0.17,'t [\Omega_{e0}^{-1}]','Units','normalized','FontSize',9)
h = colorbar();
caxis([-4 -2])
ylabel(h,'log_{10}(|\Gamma_N|)')
axis tight
text(-0.1,1.1,'(e)','Units','normalized','FontSize',10)
% xlabel(strtrim(sprintf('%s\\newline\n%s\n','h [c\Omega_{e0}^{-1}]', '(t [\Omega_{e0}^{-1}])')))


subplot(3,2,6)
mesh(htick,etick,log10(abs(ccgcc)),'EdgeColor','interp','FaceColor','interp')
view(2)
xticks(htick_pos)
xticklabels(htick_labels)
h = colorbar();
caxis([-4 -2])
ylabel(h,'log_{10}(|\Gamma_C|)')
axis tight
text(-0.1,1.1,'(f)','Units','normalized','FontSize',10)
% xlabel(strtrim(sprintf('%s\\newline\n%s\n','h [c\Omega_{e0}^{-1}]', '(t [\Omega_{e0}^{-1}])')))

%% choose a position
% h1 = -100;
% h2 = 0;
% h3 = 100;
% 
% error_tw2 = error_tw + 100; wtf_maxv
% 
% h1_indv = find(htick<h1);
% h1_ind = h1_indv(end)+1;
% h2_indv = find(htick<h2);
% h2_ind = h2_indv(end)+1;
% h3_indv = find(htick<h3);
% h3_ind = h3_indv(end)+1;
% 
% h1_t = ttick(h1_ind);
% h1_tl = h1_t - error_tw2;
% h1_th = h1_t + error_tw2;
% 
% h2_t = ttick(h2_ind);
% h2_tl = h2_t - error_tw2;
% h2_th = h2_t + error_tw2;
% 
% h3_t = ttick(h3_ind);
% h3_tl = h3_t - error_tw2;
% h3_th = h3_t + error_tw2;
% 
% h1_wl = min(ccfre(:,h1_ind));
% h1_wh = max(ccfre(:,h1_ind));
% h1_wl = dwfre(:,h1_ind);
% h1_wh = upfre(:,h1_ind);
% 
% h2_wl = min(ccfre(:,h2_ind));
% h2_wh = max(ccfre(:,h2_ind));
% h2_wl = dwfre(:,h2_ind);
% h2_wh = upfre(:,h2_ind);
% 
% h3_wl = min(ccfre(:,h3_ind));
% h3_wh = max(ccfre(:,h3_ind));
% h3_wl = dwfre(:,h3_ind);
% h3_wh = upfre(:,h3_ind);
% % -- spectrum --
% wtf = readmatrix('./csv/homo_wtf6.csv');
% wtf_t = readmatrix('./csv/homo_wtf6_time.csv');
% wtf_w = readmatrix('./csv/homo_wtf6_fre.csv');
% 
% wtfs = reshape(wtf,[length(wtf_t), length(wtf_w), 5]);
% wtf_tt = wtf_t * 1e4;
% 
% 
% figure,
% colormap('jet')
% 
% ind = 1;
% wtf_maxv = max(max(wtfs(:,:,ind)));
% subplot(3,2,1)
% mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% colorbar()
% axis tight;
% caxis([wtf_maxv-1 wtf_maxv])
% 
% 
% oft = 0;
% ofw = 2;
% wt1iv = find(wtf_tt<h1_tl);
% wt1i = wt1iv(end)+1;
% wt2iv = find(wtf_tt>h1_th);
% wt2i = wt2iv(1)-1;
% ww1iv = find(wtf_w<h1_wl);
% ww1i = ww1iv(end) + 1 - ofw;
% ww2iv = find(wtf_w>h1_wh);
% ww2i = ww2iv(1) -1 + ofw;
% 
% wt1i = 340;
% wt2i = 444;
% ww1i = 31;
% ww2i = 33;
% 
% subplot(3,2,2)
% wtf_maxv = max(max(max(wtfs(wt1i:wt2i,ww1i:ww2i,ind))));
% mesh(wtf_tt(wt1i:wt2i),wtf_w(ww1i:ww2i),transpose(wtfs(wt1i:wt2i,ww1i:ww2i,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% h = colorbar;
% ylabel(h,'log_{10}(B_w/B_0)')
% axis tight;
% caxis([wtf_maxv-0.5 wtf_maxv])
% 
% % --
% 
% ind = 3;
% wtf_maxv = max(max(wtfs(:,:,ind)));
% subplot(3,2,3)
% mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% h = colorbar;
% axis tight;
% caxis([wtf_maxv-1 wtf_maxv])
% 
% % ---
% 
% oft = 0;
% ofw = 2;
% wt1iv = find(wtf_tt<h2_tl);
% wt1i = wt1iv(end)+1;
% wt2iv = find(wtf_tt>h2_th);
% wt2i = wt2iv(1)-1;
% ww1iv = find(wtf_w<h2_wl);
% ww1i = ww1iv(end) + 1 - ofw;
% ww2iv = find(wtf_w>h2_wh);
% ww2i = ww2iv(1) - 1 + ofw;
% 
% 
% wt1i = 605;
% wt2i = 659;
% ww1i = 32;
% ww2i = 34;
% subplot(3,2,4)
% wtf_maxv = max(max(max(wtfs(wt1i:wt2i,ww1i:ww2i,ind))));
% mesh(wtf_tt(wt1i:wt2i),wtf_w(ww1i:ww2i),transpose(wtfs(wt1i:wt2i,ww1i:ww2i,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% h = colorbar;
% ylabel(h,'log_{10}(B_w/B_0)')
% axis tight;
% caxis([wtf_maxv-0.25 wtf_maxv])
% 
% %---
% 
% ind = 5;
% wtf_maxv = max(max(wtfs(:,:,ind)));
% subplot(3,2,5)
% mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% xlabel('t [\Omega_{e0}^{-1}]')
% h = colorbar;
% axis tight;
% caxis([wtf_maxv-1 wtf_maxv])
% 
% % --
% 
% oft = 0;
% ofw = 3;
% wt1iv = find(wtf_tt<h3_tl);
% wt1i = wt1iv(end)+1;
% wt2iv = find(wtf_tt>h3_th);
% wt2i = wt2iv(1)-1;
% ww1iv = find(wtf_w<h3_wl);
% ww1i = ww1iv(end) + 1 - ofw;
% ww2iv = find(wtf_w>h3_wh);
% ww2i = ww2iv(1) - 1 + ofw;
% 
% wt1i = 827;
% wt2i = 905;
% ww1i = 34;
% ww2i = 36;
% 
% subplot(3,2,6)
% wtf_maxv = max(max(max(wtfs(wt1i:wt2i,ww1i:ww2i,ind))));
% mesh(wtf_tt(wt1i:wt2i),wtf_w(ww1i:ww2i),transpose(wtfs(wt1i:wt2i,ww1i:ww2i,ind)),'EdgeColor','flat','FaceColor','flat')
% view(2)
% ylabel('\omega [\Omega_{e0}]')
% h = colorbar;
% ylabel(h,'log_{10}(B_w/B_0)')
% axis tight;
% caxis([wtf_maxv-0.25 wtf_maxv])