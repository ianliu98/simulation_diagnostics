clc
clear

%% main parameters
frel = 0.044;
freh = 0.046;

vrl = -0.272;
vrh = -0.267;

vtro = 0.013;
wtro = 0.02;

wl = frel - wtro;
wh = freh + wtro;
vl = vrl - 1 * vtro;
vh = vrh + 1 * vtro;

xoff = 50;

% bins
vedge = 2;

bin_vpar = 29;
if (vedge==1)
    v = linspace(0,1,bin_vpar);
    vt = 0.26;
    vshift = 0;
    f = 1.0 / (vt * sqrt(2*pi)) .* exp(-0.5 .* ((v-vshift)./vt).^2);
    f = f ./ max(f);
    area = sum(f);
    wid = f ./ area .* bin_vpar;
    edge_v = zeros(1,bin_vpar);
    for i=1:bin_vpar
        edge_v(i) = sum(wid(1:i));
    end
    edge_v = [0,edge_v];
    edge_v = edge_v ./ max(edge_v);
    edge_v = flip(edge_v) - 1.0;
else
    edge_v = linspace(vl,vh,bin_vpar+1);
    edge_v = flip(edge_v);
end

bin_zeta = 30;
edge_z = linspace(0,2*pi,bin_zeta);

cv = 100;
tskip = 1;
xskip = 1;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;
wavel = round(2*pi/((wh+wl)/2)/ddt);

xx = 16384-75:1:16384+75;

%% field
fy = readmatrix('./csv/homo_hole_forfy.csv');
fz = readmatrix('./csv/homo_hole_forfz.csv');

% fft
data_size = size(fy);
tt = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

pasby = zeros(data_size(1),data_size(2));
pasbz = zeros(data_size(1),data_size(2));
for i=1:data_size(2)
    pasby(:,i) = band_pass2(fy(:,i), wl, wh);
    pasbz(:,i) = band_pass2(fz(:,i), wl, wh);
end

pasby = movmean(pasby, wavel, 1);
pasbz = movmean(pasbz, wavel, 1);

clear fy fz
clear dt ifdiag tskip vedge
%% main loop

name = "./particles/1Jan1.wpia";

for num=1:1023

    num
    
    fp = fopen(name+string(num),'r');
    test_read = fread(fp,[6,Inf],"double");
    fclose(fp);
    parts = length(test_read);
    clear test_read

    % vals -> t, vx, vy, vz, x
    vals = zeros(5,parts);
    fp = fopen(name+string(num),'r');
    for i=1:parts
        tmp = fread(fp,1,"single");
        vals(:,i) = fread(fp,5,"double");
        tmp = fread(fp,1,"single");
    end
    fclose(fp);
    vals(2:4,:) = vals(2:4,:) / cv;

    % calculation
    if (num==1)
        step_ini = round(vals(1,1) / ddt);
        step_end = round(vals(1,end) / ddt);
        steps = step_end - step_ini + 1;
        hist_data = zeros(steps,bin_zeta,bin_vpar+1);
    end
    
    hist_data_tmp = zeros(steps,bin_zeta,bin_vpar+1);

    j = 1;
    i = 1;
    te = 1;
    while(j < parts)
        inds = find(vals(1,:) == vals(1,j));
        seg = vals(:,inds);
        seg(abs(seg(5,:)-16384) > xoff) = [];

        % field interpolation
        intl = seg(5,:) - floor(seg(5,:));
        intr = 1 - intl;
        indl = floor(seg(5,:)) - xx(1) + 1;
        segt = round(seg(1,1)/ddt);
        segy = pasby(segt,indl) .* intr + pasby(segt,indl+1) .* intl;
        segz = pasbz(segt,indl) .* intr + pasbz(segt,indl+1) .* intl;

        % compute zeta
        bdotv = segy .* seg(3,:) + segz .* seg(4,:);
        bcrsv = segz .* seg(3,:) - segy .* seg(4,:);
        absb  = sqrt(segy.^2 + segz.^2);
        absv  = sqrt(seg(3,:).^2 + seg(4,:).^2);
        absb(absb<1e-30) = 1e-30;
        absv(absv<1e-30) = 1e-30;
        cosvb = bdotv ./ (absb .* absv);
        cosvb(cosvb < -1) = -1;
        cosvb(cosvb > 1)  = 1;
        zeta = acos(cosvb);
        zeta(bcrsv<0) = -zeta(bcrsv<0) + 2 * pi;

        % counter interpolation
        for is=1:length(seg)
            % zeta
            ind_z = find(edge_z <= zeta(is));
            ciz_l = zeta(is) - edge_z(ind_z(end));
            if (ind_z(end)+1 > bin_zeta)
                continue
            end
            ciz_h = edge_z(ind_z(end)+1) - zeta(is);
            % vpara
            ind_v = find(edge_v >= seg(2,is));
            if (isempty(ind_v))
                te=te+1;
                continue
            end
            civ_l = edge_v(ind_v(end)) - seg(2,is);
            if (ind_v(end)+1 > bin_vpar+1)
                continue
            end
            civ_h = seg(2,is) - edge_v(ind_v(end)+1);
            
            ss = (edge_z(ind_z(end)+1) - edge_z(ind_z(end))) * ...
                (edge_v(ind_v(end)+1) - edge_v(ind_v(end)));
            sll = ciz_l * civ_l / ss;
            slr = ciz_h * civ_l / ss;
            shl = ciz_l * civ_h / ss;
            shr = ciz_h * civ_h / ss;
            hist_data_tmp(i,ind_z(end),ind_v(end)) = ...
                hist_data_tmp(i,ind_z(end),ind_v(end)) + shr;
            hist_data_tmp(i,ind_z(end)+1,ind_v(end)) = ...
                hist_data_tmp(i,ind_z(end)+1,ind_v(end)) + shl;
            hist_data_tmp(i,ind_z(end),ind_v(end)+1) = ...
                hist_data_tmp(i,ind_z(end),ind_v(end)+1) + slr;
            hist_data_tmp(i,ind_z(end)+1,ind_v(end)+1) = ...
                hist_data_tmp(i,ind_z(end)+1,ind_v(end)+1) + sll;
        end
        
        %[cunts, ~, ~] = histcounts2(zeta,seg(2,:),edges_zeta,edges_vpar);
        %hist_data_tmp(i,:,:) = cunts;

        j = j+length(inds);
        i = i + 1;
    end

    hist_data = hist_data + hist_data_tmp;
    
end

hist_data_perm = permute(hist_data, [3 2 1]);
%% data manupulation
hist_init_f = load('./mat/hole_init_0124_ave.mat');
hist_init = hist_init_f.data_ave;
hist_init(hist_init>-1e-16) = -1e-16;

meanwindow = 10;

hist_data_perm2 = movmean(hist_data_perm,meanwindow,3);

poten = (abs(hist_data_perm2) - abs(hist_init)) ./ abs(hist_init);

%% video
clc
v = VideoWriter('./videos/test0124_7.avi');
open(v)
figure,
colormap(jet)
for k = 1:5:steps
    k
    mesh(edge_z,edge_v,poten(:,:,k),'EdgeColor','interp','FaceColor','interp')
    view(2)
    yline(vrh,'--')
    yline(vrl,'--')
    xlabel('\zeta')
    ylabel('v_{||}')
    xticks([0 pi/2 pi pi*1.5 2*pi-0.1])
    xticklabels({'0','\pi/2', '\pi', '3\pi/2', '2\pi'})
    title("t = " + string(floor(tt(step_ini+k-1))) + " [\Omega_{e0}^{-1}]")
    colorbar()
    caxis([-0.5 0.5])
    axis tight;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);