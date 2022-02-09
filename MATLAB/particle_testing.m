clc
clear

%% theoretical values
bbound = 1.0;
thop = 1;
wpe = 15.0;

[vvr,vtr,ww,wtr] = theory_vector(bbound, thop, wpe);  % modification may be needed


%% read data
pnts = 6;

vxs = readmatrix('./csv/pti1_vx.csv');
vys = readmatrix('./csv/pti1_vy.csv');
vzs = readmatrix('./csv/pti1_vz.csv');
xes = readmatrix('./csv/pti1_xe.csv');
t = readmatrix('./csv/pti1_t.csv');

cv = 100;
tskip = 4;
xskip = 4;
ifdiag = 256;
dt = 0.004;
ddt = dt * tskip * ifdiag;
by = readmatrix('./csv/0104_homo_by.csv');
bz = readmatrix('./csv/0104_homo_bz.csv');
ey = readmatrix('./csv/homo_ey');
ez = readmatrix('./csv/homo_ez');
jy = readmatrix('./csv/homo_jy');
jz = readmatrix('./csv/homo_jz');
data_size = size(by);

tickh = ((1:1:data_size(2)) - data_size(2)/2) .* xskip ./ cv;  % -150 ~ 150
tickt = (1:1:data_size(1)) .* (dt * tskip * ifdiag);

%% analysis

prti = 1;

fvx = vxs(prti, vxs(prti,:) < 0);
fvy = vys(prti, vxs(prti,:) < 0);
fvz = vzs(prti, vxs(prti,:) < 0);
fxe = xes(prti, vxs(prti,:) < 0);
ftt = t(vxs(prti,:) < 0);

fww  = zeros(1,length(fvx));
fwtr = zeros(1,length(fvx));
fxe_interp = zeros(4,length(fvx));
ftt_interp = zeros(4,length(fvx));

fxe = fxe / xskip;
dam1 = 2048 / xskip;
dam2 = 32768 / xskip - dam1;
for i=1:length(fvx)
    indv = find(vvr < fvx(i));
    if (indv(end)==9000)
        fww(i) = ww(end);
        fwtr(i) = wtr(end);
    else
        fww(i) = ww(indv(end)+1); 
        fwtr(i) = wtr(indv(end)+1);
    end
    % space
    if ((fxe(i)<=dam1) || (fxe(i)>=dam2))
        fxe_interp(:,i) = -1;
    else
        fxe_tmp = fxe(i) - dam1;
        xli = floor(fxe_tmp);
        xri = xli + 1;
        lfi = fxe_tmp - xli;
        rfi = 1 - lfi;
        if (xli==0)
            xli = 2;
        end
        fxe_interp(1,i) = xli;
        fxe_interp(2,i) = xri;
        fxe_interp(3,i) = lfi;
        fxe_interp(4,i) = rfi;
    end
    % time
    ftt_tmp_v = find(tickt < ftt(i));
    if (isempty(ftt_tmp_v))
        ftt_interp(:,i) = -1;
        continue
    end
    if (ftt_tmp_v(end) == data_size(1))
        ftt_interp(:,i) = -1;
    else
        ftt_li = ftt_tmp_v(end);  ftt_hi = ftt_li + 1;
        ftt_l = tickt(ftt_li);  ftt_h = tickt(ftt_hi);
        dfi = (ftt(i) - ftt_l)/ddt;  ofi = (ftt_h - ftt(i))/ddt;
        ftt_interp(1,i) = ftt_li;
        ftt_interp(2,i) = ftt_hi;
        ftt_interp(3,i) = dfi;
        ftt_interp(4,i) = ofi;
    end
end


%% field

wtr_w = 1; % important parameter

frel = fww - wtr_w .* fwtr;
freh = fww + wtr_w .* fwtr;
% frel = ones(1,length(fvx)) * 0.04;
% freh = ones(1,length(fvx)) * 0.06;

zetas = zeros(1,length(fvx));
jbs   = zeros(1,length(fvx));
jes   = zeros(1,length(fvx));
bfld  = zeros(1,length(fvx));
efld  = zeros(1,length(fvx));

for i=1:length(fvx)
    if ((fxe_interp(1,i)==-1) || (ftt_interp(1,i)==-1))
        continue
    end
    
    % b
    pasby1 = band_pass(by(:, fxe_interp(1,i)), frel(i), freh(i));
    pasby2 = band_pass(by(:, fxe_interp(2,i)), frel(i), freh(i));
    pasbz1 = band_pass(bz(:, fxe_interp(1,i)), frel(i), freh(i));
    pasbz2 = band_pass(bz(:, fxe_interp(2,i)), frel(i), freh(i));
    % e
    pasey1 = band_pass(ey(:, fxe_interp(1,i)), frel(i), freh(i));
    pasey2 = band_pass(ey(:, fxe_interp(2,i)), frel(i), freh(i));
    pasez1 = band_pass(ez(:, fxe_interp(1,i)), frel(i), freh(i));
    pasez2 = band_pass(ez(:, fxe_interp(2,i)), frel(i), freh(i));
    % j
    pasjy1 = band_pass(jy(:, fxe_interp(1,i)), frel(i), freh(i));
    pasjy2 = band_pass(jy(:, fxe_interp(2,i)), frel(i), freh(i));
    pasjz1 = band_pass(jz(:, fxe_interp(1,i)), frel(i), freh(i));
    pasjz2 = band_pass(jz(:, fxe_interp(2,i)), frel(i), freh(i));    
    
    
    % interpolate
    stl = fxe_interp(3,i) * ftt_interp(3,i);
    str = fxe_interp(4,i) * ftt_interp(4,i);
    sll = fxe_interp(3,i) * ftt_interp(3,i);
    slr = fxe_interp(4,i) * ftt_interp(3,i);
   
    % by
    bytl = pasby1(ftt_interp(2,i));  bytr = pasby2(ftt_interp(2,i));
    byll = pasby1(ftt_interp(1,i));  bylr = pasby2(ftt_interp(1,i));
    pby = byll * str + bylr * stl + bytl * slr + bytr * sll;
    % bz
    bztl = pasbz1(ftt_interp(2,i));  bztr = pasbz2(ftt_interp(2,i));
    bzll = pasbz1(ftt_interp(1,i));  bzlr = pasbz2(ftt_interp(1,i));
    pbz = bzll * str + bzlr * stl + bztl * slr + bztr * sll;
    % ey
    eytl = pasey1(ftt_interp(2,i));  eytr = pasey2(ftt_interp(2,i));
    eyll = pasey1(ftt_interp(1,i));  eylr = pasey2(ftt_interp(1,i));
    pey = eyll * str + eylr * stl + eytl * slr + eytr * sll;
    % ez
    eztl = pasez1(ftt_interp(2,i));  eztr = pasez2(ftt_interp(2,i));
    ezll = pasez1(ftt_interp(1,i));  ezlr = pasez2(ftt_interp(1,i));
    pez = ezll * str + ezlr * stl + eztl * slr + eztr * sll;
    % jy
    jytl = pasjy1(ftt_interp(2,i));  jytr = pasjy2(ftt_interp(2,i));
    jyll = pasjy1(ftt_interp(1,i));  jylr = pasjy2(ftt_interp(1,i));
    pjy = jyll * str + jylr * stl + jytl * slr + jytr * sll;
    % jz
    jztl = pasjz1(ftt_interp(2,i));  jztr = pasjz2(ftt_interp(2,i));
    jzll = pasjz1(ftt_interp(1,i));  jzlr = pasjz2(ftt_interp(1,i));
    pjz = jzll * str + jzlr * stl + jztl * slr + jztr * sll;
    
    % computer jb je
    jbs(i) = (pjy * pby + pjz * pbz) / sqrt(pby^2 + pbz^2);
    jes(i) = (pjy * pey + pjz * pez) / sqrt(pey^2 + pez^2);
    
    % compute bfield efield
    bfld(i) = sqrt(pby^2 + pbz^2);
    efld(i) = sqrt(pey^2 + pez^2);
    
    % compute zeta
    bdotv = pby .* fvy(i) + pbz .* fvz(i);
    bcrsv = pbz .* fvy(i) - pby .* fvz(i);
    absb  = sqrt(pby^2 + pbz^2);
    absv  = sqrt(fvy(i)^2 + fvz(i)^2);
    absb(absb<1e-30) = 1e-30;
    absv(absv<1e-30) = 1e-30;
    cosvb = bdotv ./ (absb .* absv);
    cosvb(cosvb < -1) = -1;
    cosvb(cosvb > 1)  = 1;
    zeta = acos(cosvb);
    zeta(bcrsv<0) = -zeta(bcrsv<0) + 2 * pi;
    
    zetas(i) = zeta;
end

clear by bz ey ez jy jz
%% plot
%uzetas = unwrap(zetas);
% sep_iv = find(diff(fxe)>1000);
% if (fvx(1)>0)
%     tm1s = sep_iv(1:2:end);
%     tm2s = sep_iv(2:2:end);
%     if (fvx(end)<0)
%         tm2s = [tm2s, length(fvx)];
%     end
% else
%     tm1s = sep_iv(2:2:end);
%     tm1s = [1,tm1s];
%     tm2s = sep_iv(1:2:end);
%     if (fvx(end)<0)
%         tm2s = [tm2s, length(fvx)];
%     end
% end
% 
% seps = zeros(2*length(tm1s),1);
% seps(1:2:end) = tm1s;
% seps(2:2:end) = tm2s;
% 
% for iv=2:length(seps)
%     
%     tm1 = seps(iv-1) + 2;
%     tm2 = seps(iv) - 2;
%     
    % rotate
%     maxvx_r = max(abs(fvx(tm1:tm2)));
%     v = VideoWriter("./videos/rotate"+ string(iv)+".avi");
%     open(v)
%     figure,
%     for jj = tm1:tm2
%         [xr,yr] = pol2cart(zetas(jj),abs(fvx(jj)));
%         plot([0 xr],[0 yr],'r')
%         title("t="+string(ftt(jj))+"  v_{||}="+string(round(fvx(jj)))+"  h="+string(round(fxe(jj))))
%         grid on
%         axis([-maxvx_r maxvx_r -maxvx_r maxvx_r])
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     end
%     close(v);
    
    % trace
    maxvx = max(fvx);
    minvx = min(fvx);
    v = VideoWriter("./videos/pt_small_wtr"+string(prti)+".avi");
    open(v)
    hr = [];
    vr = [];
    figure,
    for jj = 1:length(ftt)
        hr = [hr, zetas(jj)];
        vr = [vr, fvx(jj)];
        if (length(hr)>=2)
            if (abs(hr(end)-hr(end-1))<5)
                plot(hr,vr,'r-')
            else
                hr = [];
                vr = [];
                continue
            end
        end
        grid on
%         title("t="+string(ftt(jj))+"[\Omega_{e0}^{-1}]  h="+ ...
%             string(round(fxe(jj)-16384/4)/25) + "[c\Omega_{e0}^{-1}]")
        title("t="+string(ftt(jj))+" h="+ string(round(fxe(jj)-16384/4)/25))
        axis([0 2*pi minvx maxvx])
        %xlim([0 2*pi])
        %ylim([mean(fvx)-2 mean(fvx)+2])
        xlabel('\zeta')
        ylabel('V_{||}')
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    close(v);
    

    
% end

% ----

% maxvx = max(fvx(tm1:tm2));
% minvx = min(fvx(tm1:tm2));
% 
% v = VideoWriter('./videos/pti_vzeta.avi');
% open(v)
% figure,
% for jj = tm1:tm2
%     scatter(zetas(jj),fvx(jj),'r*')
%     grid on
%     title("t="+string(ftt(jj))+"  xe="+string(round(fxe(jj))))
%     axis([0 2*pi minvx maxvx])
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% close(v);


% ---
%% test
clc

pnth = [90.08 80.33 73.75 60.41 50.07 39.02 20.21];
pntt = [18933.8 18973.3 19001.3 19057.5 19101.7 19148.8 19228.7];
pntz = [9000.08 8000.33 7300.75 6000.41 5000.07 3900.02 2000.21];
txt= [" A"," B"," C"," D"," E"," F"," G"];
txt2= [" A"," B"," C"," D"," E"," F"," G"];
pind = zeros(7,1);
for i=1:7
    aa = find(ftt<pntt(i));
    pind(i) = aa(end)+1;
end



% particle 1
starts = [77   1626 3160 4561 6106 7445 8800];
ends =   [1421 2960 4468 5961 7301 8658 10052];
starts = [77   1626 3160 4561 6106 7445 9050];  % thesis A ~ G
ends =   [1421 2960 4468 5961 7301 8658 9340];
starts = [77   1626 3160 4561 6106 7445 9090];  % thesis B ~ D
ends =   [1421 2960 4468 5961 7301 8658 9173];
prts = 7;

% configure
% v_|| - zeta
maxvx = max(fvx(starts(prts):ends(prts)));
minvx = min(fvx(starts(prts):ends(prts)));
zr = [];
vr = [];

% jb je
maxjb = max(jbs(starts(prts):ends(prts)));
maxje = max(jes(starts(prts):ends(prts)));
maxj = max([maxjb maxje]);
minjb = min(jbs(starts(prts):ends(prts)));
minje = min(jes(starts(prts):ends(prts)));
minj = min([minjb minje]);

% kinetic energy
gamma = cv ./ sqrt(cv^2 - (fvx.^2 + fvy.^2 + fvz.^2));
gamma  = (gamma-1) .* 511;

% zeta rotation
maxvx_r = max(abs(fvx(starts(prts):ends(prts))));

% bf ef
fld = sqrt(bfld.^2 + efld.^2);
maxb = max(bfld(starts(prts):ends(prts)));
maxe = max(efld(starts(prts):ends(prts)));
minb = min(bfld(starts(prts):ends(prts)));
mine = min(efld(starts(prts):ends(prts)));
maxf = max([maxb maxe]);
minf = min([minb mine]);

maxf = max(fld(starts(prts):ends(prts)));
minf = min(fld(starts(prts):ends(prts)));

% video?
% v = VideoWriter("./videos/test0123"+string(prti)+".avi");
% open(v)

% fig = figure('visible','off');
fig = figure;
fig.Position = [200 500 960 720];
for jj = starts(prts):ends(prts)%length(ftt)
    subplot(2,2,1)
    zr = [zr, zetas(jj)];
    vr = [vr, fvx(jj)/100];
    if (length(zr)>=2)
        if (abs(zr(end)-zr(end-1))<5)
            plot(zr,vr,'r-')
        else
            zr = [];
            vr = [];
            continue
        end
    end
    hold on
    tmp = find(pind==jj);
    if (tmp)
        text(zetas(jj), fvx(jj)/100, txt(tmp))
    end
    grid on
    box on
    axis([0 2*pi minvx/cv maxvx/cv])

    xlabel('\zeta')
    ylabel('V_{||} [c]')
    if (jj == ends(prts))
        text(-0.1,1.1,'(a)','Units','normalized','FontSize',10)
    end
%     title("V_{||}="+string(round(fvx(jj)/100,4))+"[c] \zeta="+string(round(zetas(jj),4)))
    
    subplot(2,2,2)
    plot([ftt(jj-1) ftt(jj)],[jbs(jj-1) jbs(jj)],'r')
    hold on
    plot([ftt(jj-1) ftt(jj)],[jes(jj-1) jes(jj)],'b')
    grid on
    ylim([ minj maxj])
    tmp = find(pind==jj);
    if (tmp)
        xline(ftt(jj),'g--')
        text(ftt(jj), minj+0.005, txt(tmp))
    end
%     axis([ftt(starts(prts)) ftt(ends(prts)) minj maxj])
%     plot([0 jbs(jj)],[0 jes(jj)],'r',jbs(jj),jes(jj),'ro')
%     title("J_B="+string(round(jbs(jj),4))+"[q] J_E="+string(round(jes(jj),4))+"[q]")
%      grid on
%      axis([-maxjb maxjb -maxje maxje])
    xlabel('t [c\Omega_{e0}^{-1}]')
%     ylabel('J')
%     xlabel('J_B [q]')
%     ylabel('J_E [q]')
    ylabel('J [q]')
    legend('J_B','J_E')
    if (jj == ends(prts))
        text(-0.1,1.1,'(b)','Units','normalized','FontSize',10)
        xlim([ftt(starts(prts)) ftt(ends(prts))])
    end
    
    subplot(2,2,3)
    plot([ftt(jj-1) ftt(jj)],[gamma(jj-1) gamma(jj)],'r')
    hold on
    grid on
    %axis([ftt(starts(prts)) ftt(ends(prts)) min(gamma(starts(prts):ends(prts))) ...
    %    max(gamma(starts(prts):ends(prts)))])
    ylim([min(gamma(starts(prts):ends(prts))) max(gamma(starts(prts):ends(prts)))])
    tmp = find(pind==jj);
    if (tmp)
        xline(ftt(jj),'g--')
        text(ftt(jj), max(gamma(starts(prts):ends(prts)))-4e-3, txt(tmp))
    end
    ylabel('K_p [mc^2]')
    xlabel('[c\Omega_{e0}^{-1}]')
    if (jj == ends(prts))
        xlim([ftt(starts(prts)) ftt(ends(prts))])
        text(-0.1,1.1,'(c)','Units','normalized','FontSize',10)
    end
%     title("t="+string(round(ftt(jj),4))+" [\Omega_{e0}^{-1}]"+" h="+ ...
%         string(round((fxe(jj)-16384/4)/25,2)) + "[c\Omega_{e0}^{-1}]")
    
    subplot(2,2,4)
%     plot([ftt(jj-1) ftt(jj)],[log10(fld(jj-1)) log10(fld(jj))],'black')
    plot([ftt(jj-1) ftt(jj)],[bfld(jj-1) bfld(jj)],'r')
    hold on
%     plot([ftt(jj-1) ftt(jj)],[log10(efld(jj-1)) log10(efld(jj))],'r')
    grid on
%     axis([ftt(starts(prts)) ftt(ends(prts)) minb maxb])
    ylim([minb maxb])
    tmp = find(pind==jj);
    if (tmp)
        xline(ftt(jj),'g--')
        text(ftt(jj), maxb-5e-5, txt(tmp))
    end
    xlabel('t [c\Omega_{e0}^{-1}]')
    ylabel('B_w [B_0]')
    if (jj == ends(prts))
        xlim([ftt(starts(prts)) ftt(ends(prts))])
        text(-0.1,1.1,'(d)','Units','normalized','FontSize',10)
    end
%     ylabel('log_{10}(\Omega_{w})')
%     [xr,yr] = pol2cart(zetas(jj),fvx(jj));
%     plot([0 xr],[0 yr],'r',xr,yr,'ro')
%     grid on
%     axis([-maxvx_r maxvx_r -maxvx_r maxvx_r])
%     title("B_w="+string(bfld(jj))+" [B_0]")
%     title("pb="+string(round(log10(bfld(jj)),4)) + " pe="+string(round(log10(efld(jj)),4)))
%     legend('B_w','E_w')
    
    drawnow
%     frame = getframe(gcf);
%     writeVideo(v,frame);
%     im{jj} = frame2im(frame);
end
% close(v)

