% big circle
clc
clear

% configure
ddw = 0.002;
wstart = 0.004;
wend = 0.1;
ncir = (wend - wstart) / ddw + 1;

% parameter
tskip = 4;
xskip = 4;
cv = 100;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;

% load data
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');
ey = readmatrix('./csv/homo_ey');
ez = readmatrix('./csv/homo_ez');
jy = readmatrix('./csv/homo_jy');
jz = readmatrix('./csv/homo_jz');

data_size = size(by);
xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;
ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

for iw=1:ncir
    frel = wstart + ddw * (iw - 1);
    freh = wstart + ddw * iw;
    
    % fft
    pasby = band_pass(by, frel, freh);
    pasbz = band_pass(bz, frel, freh);
    pasey = band_pass(ey, frel, freh);
    pasez = band_pass(ez, frel, freh);
    pasjy = band_pass(jy, frel, freh);
    pasjz = band_pass(jz, frel, freh);
    
    % field
    fld = sqrt(pasby.^2 + pasbz.^2);
    
    % jb je
    jb = (pasjy .* pasby + pasjz .* pasbz) ./ sqrt(pasby.^2 + pasbz.^2);
    je = (pasjy .* pasey + pasjz .* pasez) ./ sqrt(pasey.^2 + pasez.^2);
    
    % fre
    fre = zeros(data_size(1),data_size(2));
    for j=1:data_size(2)
        fre(:,j) = Nogi_frequency(pasby(:,j), pasbz(:,j), ddt);
    end
    
    % clear
    clear pasby pasbz pasey pasez pasjy pasjz
    
    % dwdt
    dwdt = zeros(data_size(1),data_size(2));
    mean_wind = 50;
    for j=1:data_size(2)
        dwdt_tmp = four_order_appr(fre(:,j),ddt);
        dwdt_tmp = [dwdt_tmp;0;0;0;0];
        dwdt(:,j) = movmean(dwdt_tmp,mean_wind);
    end
    
    % S
    [s0, s1, ~] = inhomogeneity(fre);
    S = -1.0 * s1 .* dwdt ./ (s0 .* fre .* fld);
    
    % nonlinear
    Gn = Nonlinear_growthrate(fld, fre);
    
    % plot & save
    plots(fld,"fld",frel,[min(min(fld)) max(max(fld))])
    plots(fre,"fre",frel,[frel freh])
    plots(jb,"jb",frel,[min(min(jb)) max(max(jb))])
    plots(je,"je",frel,[min(min(je)) max(max(je))])
    plots(dwdt,"dwdt",frel,[-1e-5 1e-5])
    plots(S,"S",frel,[-5 5])
    plots(Gn,"Gn",frel,[1e-5 1e-2])
end

function plots(data,dtitle,label,drange)
    % paramters
    skip_x = 16;
    skip_t = 10;
    tskip = 4;
    xskip = 4;
    cv = 100;
    ifdiag = 256;
    dt = 0.004;
    data_size = size(data);
    xtick = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;
    ytick = (1:1:data_size(1)) .* (dt * tskip * ifdiag);
    
    figure('visible','off'),
    colormap('jet')
    mesh(xtick(1:skip_x:end),ytick(1:skip_t:end),data(1:skip_t:end,1:skip_x:end),'EdgeColor','interp','FaceColor','interp')
    view(2)
    colorbar()
    caxis(drange)
    xlabel('h [c\Omega_{e0}^{-1}]')
    ylabel('t [\Omega_{e0}^{-1}')
    title(dtitle)
    axis tight;
    print("./pdfs/"+dtitle+string(label),'-dpdf','-bestfit')
    close(gcf)
    
end