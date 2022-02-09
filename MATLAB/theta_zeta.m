% plot theta-zeta 
clc
clear

%% parameters
% main parameters
w   = 0.05;   % frequency of wave
we  = 1.0;    % gyrofrequency at equator
wpe = 15.0;   % cold plasma frequency
c   = 1.0;    % speed of light

utpa = 0.26;  % [c] thermal parallel momentum
utpe = 0.30;  % [c]
beta = 0.3;   % shape parameter in subtracted-Maxwellian


% derived parameters
xi  = sqrt(abs( w * (we - w) / wpe^2));
chi = sqrt(abs( 1.0 / (1.0 + xi^2)));
k   = w / (c * chi * xi);

uperph  = sqrt(pi/2) * ((1 - beta^(1.5) / (1 - beta))) * utpe;
vperp  = c / sqrt(c^2 + (utpa^2 + uperph^2)) * uperph;  % perpendicular velocity
vpara  = c / sqrt(c^2 + (utpa^2 + uperph^2)) * utpa;    % parallel velocity
gamma  = 1 / sqrt(1 - (vperp^2 + vpara^2)/c^2);

% theoretical wave amplitudes
thop = 1;   % 1 -> optimum  others -> threshold
b02 = 1.0;  % homogeneous
amp = theoretical_amp(w, thop, b02);


%% calculation & plot
wt  = sqrt(abs(k * vperp * amp));
wtr = wt * chi / sqrt(abs(gamma));
vtr = 2.0 * wtr / k;

nzeta = 1000;
ncc = 20;

zeta = linspace(0,2*pi,nzeta);

S = -0.4;
cc = linspace(-5e-3,5e-3,ncc);

data_save = zeros(ncc,nzeta);

trans_ind = 0;
zeta0_tmp = 0;
figure,
for i=1:length(cc)
    theta_over_wtr_2 = cc(i) / wtr^2 - 2.0 * (cos(zeta) - S * zeta);
    % find zeta1 & zeta2
    neg_range = find(theta_over_wtr_2 < 0);
    if (max(diff(neg_range)) > 10)
        ind1 = find(diff(neg_range) == max(diff(neg_range)));
        ind2 = neg_range(ind1);
        ind3 = neg_range(ind1+1);
        trans_ind = i;
        if (zeta0_tmp == 0)
            pos_val = theta_over_wtr_2(theta_over_wtr_2 > 0);
            zeta0_1 = find(diff(pos_val) == min(diff(pos_val)));
            zeta0_2 = find(diff(pos_val) == max(diff(pos_val)));
            ind0 = round(ind2 + (zeta0_1+zeta0_2)/2);
            zeta0_tmp = 1;
        end
    end
    
    theta_over_wtr_2(theta_over_wtr_2<0) = 0;
    vert = sqrt(theta_over_wtr_2) / 2.0;
    data_save(i,:) = vert;
    %vert(theta_over_wtr_2<0) = -vert(theta_over_wtr_2<0);
    vert_neg = -vert;
    plot(zeta,vert,'black',zeta,vert_neg,'black')
    hold on
end
%xline(zeta(ind0),'--')
%xline(zeta(ind2),'--');  % zeta1
%xline(zeta(ind3),'--');  % zeta2
xline(pi)
xticks([0 pi/2 pi pi*1.5 6.28])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
vw = 0.658;
yticks([-3.5-vw -3.5 -3.5+vw 0-vw 0 0+vw 3.5-vw 3.5 3.5+vw]);
yticklabels({'V_R^{i-1}-V_{tr}^{i-1}','V_R^{i-1}','V_R^{i-1}+V_{tr}^{i-1}', ...
    'V_R^i-V_{tr}^i','V_R^i','V_R^i+V_{tr}^i', 'V_R^{i+1}-V_{tr}^{i+1}', ...
    'V_R^{i+1}','V_R^{i+1}+V_{tr}^{i+1}'});
xlabel('\zeta')
%ylabel('\theta/2\omega_{tr}')
axis tight;