clc
clear

%%
frel = 0.08;
freh = 0.084;

tskip = 4;
xskip = 4;
ifdiag = 256;
dt = 0.004;
cv = 100;
ddt = dt * tskip * ifdiag;

mean_wind =     50;     % movmean window size

%% 
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');

ey = readmatrix('./csv/homo_ey');
ez = readmatrix('./csv/homo_ez');

jy = readmatrix('./csv/homo_jy');
jz = readmatrix('./csv/homo_jz');

%%
data_size = size(by);
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

pasby = band_pass(by, frel, freh);
pasbz = band_pass(bz, frel, freh);

pasey = band_pass(ey, frel, freh);
pasez = band_pass(ez, frel, freh);

pasjy = band_pass(jy, frel, freh);
pasjz = band_pass(jz, frel, freh);

ajb = (pasjy .* pasby + pasjz .* pasbz) ./ sqrt(pasby.^2 + pasbz.^2);
aje = (pasjy .* pasey + pasjz .* pasez) ./ sqrt(pasey.^2 + pasez.^2);

%% 
afre = zeros(data_size(1),data_size(2));
for j=1:data_size(2)
    afre(:,j) = Nogi_frequency(pasby(:,j), pasbz(:,j), ddt);
end

aforf = sqrt(pasby.^2 + pasbz.^2);

%% source velocity
vg = Group_velocity(afre)./cv;
vr = Resonance_velocity(afre);

vs = vg - vr;

%%
figure,
subplot(4,1,1)
plot(ytick,movmean(aforf(:,5586),wavel))
xlim([ytick(1) ytick(end)])

subplot(4,1,2)
plot(ytick,movmean(afre(:,5588),wavel))
xlim([ytick(1) ytick(end)])

subplot(4,1,3)
plot(ytick,movmean(ajb(:,5588),wavel))
xlim([ytick(1) ytick(end)])
yline(0)
title('jb')

subplot(4,1,4)
plot(ytick,movmean(aje(:,5586),wavel))
xlim([ytick(1) ytick(end)])
yline(0)
title('je')


%%
pos = 3530;

figure,
subplot(6,1,1)
plot(ytick,movmean(ajb(:,pos),wavel))
xlim([ytick(1) ytick(end)])
yline(0)

subplot(6,1,2)
plot(ytick,movmean(aje(:,pos),wavel))
xlim([ytick(1) ytick(end)])
yline(0)

subplot(6,1,3)
plot(ytick,movmean(aje(:,pos-1),wavel))
xlim([ytick(1) ytick(end)])
yline(0)

subplot(6,1,4)
plot(ytick,movmean(aje(:,pos-2),wavel))
xlim([ytick(1) ytick(end)])
yline(0)

subplot(6,1,5)
plot(ytick,movmean(aje(:,pos-3),wavel))
xlim([ytick(1) ytick(end)])
yline(0)

subplot(6,1,6)
plot(ytick,movmean(aje(:,pos-4),wavel))
xlim([ytick(1) ytick(end)])
yline(0)


%% 
figure,
colormap(jet)
mesh(xtick(1:8:end),ytick(1:5:end),vs(1:5:end,1:8:end),'FaceColor','interp','EdgeColor','interp')
view(2)
colorbar()
axis tight

figure,
colormap(jet)
mesh(xtick(1:8:end),ytick(1:5:end),aforf(1:5:end,1:8:end),'FaceColor','interp','EdgeColor','interp')
view(2)
colorbar()
axis tight
