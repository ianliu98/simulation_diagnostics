% select single packet from dynamic spectrum

clc
clear

wtf = readmatrix('./csv/homo_wtf.csv');
wtf_t = readmatrix('./csv/homo_wtf_time.csv');
wtf_w = readmatrix('./csv/homo_wtf_fre.csv');

wtfs = reshape(wtf,[length(wtf_t), length(wtf_w), 5]);
maxv = max(wtf);

tt = wtf_t * 1e4;

ind = 5;
figure,
colormap('jet')
mesh(tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','interp','FaceColor','interp')
view(2)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\omega [\Omega_{e0}]')
colorbar()
axis tight;
caxis([maxv-1 maxv])

points = 2;
[tp, wp] = ginput(points);

wl = min(wp);
wh = max(wp);
tl = min(tp);
th = max(tp);

