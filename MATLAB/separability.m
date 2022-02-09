function [delta_w, ww]=separability(bbound,thop, wpe)

% parameter
wp = wpe;    % cold plasma frequency
nhnc = 4e-4;    % from definition of plasma frequency nh/nc = wph^2 / wpc^2
wph = wp*sqrt(nhnc);    % hot plasma frequency
Q = 0.1;   % depth of electron hole
cv = 100; % speed of light
b02 = 1.00;
b0 = bbound;
nx = 327.68;
nx = nx * cv;
a = (b02./b0 - 1) / (nx/2)^2;    % gradient
tau = 0.5;  % transition time
dw =0.0001; % resolution
fce = 1.000;    % electron gyrofrequency
utpara = 26;  % thermal velocity momentum
utperp = 30;
beta = 0.3; % parameter for subtracted maxwellian

uperp = sqrt(pi/2) * ((1 - beta^(1.5) / (1 - beta))) * utperp;  % (79) averaged perpendicular thermal velocity momentum
vperp = cv / sqrt(cv^2 + (utpara^2 + uperp^2)) * uperp; % perpendicular velocity
vpara = cv / sqrt(cv^2 + (utpara^2 + uperp^2)) * utpara; % parallel velocity
gamma2  = 1 / sqrt(1 - (vperp^2 + vpara^2)/cv^2);

for i=1:9000
    for j =1:1
        ww(i) = 0.0+dw*i;
        w = ww(i);
        xi2 = w*(fce - w) / (wp^2);    % (4) -> Omega_e = 1
        chi2 = 1 / (1+ xi2);  % (6) 
        chi = sqrt (chi2);
        cchi2(i) = chi2;
        cchi(i) = chi;
        xi = sqrt (xi2);
        xxi2(i) = xi2;
        xxi(i) = xi;
        vp = chi*xi*cv;  % (7)
        vg = (cv*xi/chi)/(xi^2 + fce/(2*(fce-w)));   % (8)
        % standardlization
        w_s = w / fce;
        vperp_s = vperp / cv;
        vp_s = vp / cv;
        vg_s = vg / cv;
        uperp_s = uperp / cv;
        wph_s = wph / fce;
        utpara_s = utpara / cv;
        a_s = a(j)*cv^2 / fce^2;
        
        vr_s = (w_s^2 - sqrt(w_s^4 + (w_s^2 + vp_s^2)*(1 - w_s^2 - vperp_s^2)))/(w_s^2/vp_s + vp_s);    % (23)
        vr = vr_s * cv;
        gamma = cv/sqrt(cv^2 - (vr^2+vperp^2)); % (12)
        vr = vp*(1 - fce/(w*gamma)); % (22)
        s0 = chi*vperp/(xi*cv);    % (38)
        ss0(i) = s0;
        s1 = gamma * (1 - vr/vg)^2;     % (39)
        ss1(i) = s1;
        s_tmp = gamma*w*((vperp/cv)^2)/fce - ( 2 + (w/fce)*chi^2*(fce - gamma*w)/(fce -w) )*vr*vp/cv^2; 
        s2 = s_tmp/(2*xi*chi);   % (40)
        ss2(i) = s2;
        
        % optimum
        bw_eps_1 = 0.8*pi^(-5/2)*Q*vp_s*vg_s*uperp_s*wph_s^2 / (tau(1)*w_s*utpara_s);
        bw_eps_2 = (1 - vr_s/vg_s)^2 * exp(-0.5 * (gamma*vr_s/utpara_s)^2);
        bw_eps(i, j) = bw_eps_1 * bw_eps_2;  % (97)
         
        % threshold
        bwth_eps_1 = 100*pi^3*gamma^4*xi*(a_s*s2*utpara_s/Q)^2 / (w_s*wph_s^4*(chi*uperp_s)^5);
        bwth_eps(i, j) = bwth_eps_1 * exp((gamma*vr_s/utpara_s)^2); % (105)
    end
end

cc = cv / 100;
if (thop == 1)
    amps = bw_eps;
else
    amps = bwth_eps;
end

k   = ww ./ (cc .* cchi .* xxi);
wt  = sqrt(abs(k .* transpose(amps) * vperp_s));
wtr = wt .* cchi / sqrt(abs(gamma2));
delta_w = 4.0 * wtr ./ (1.0 + cchi2 .* (fce ./ (gamma2 * ww) - 1.0) .* (xxi2 + fce ./ (2.0 * (fce - ww)))); 

end