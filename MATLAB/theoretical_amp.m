function amp = theoretical_amp(fre, thop, bbound)

% parameter
wp = 15;    % cold plasma frequency
nhnc = 4e-4;    % from definition of plasma frequency nh/nc = wph^2 / wpc^2
wph = wp*sqrt(nhnc);    % hot plasma frequency
Q = 0.1;   % depth of electron hole
cv = 100; % speed of light
b02 = bbound;
b0 = 1.0;
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

for i=1:9000
    for j =1:1
        ww(i) = 0.0+dw*i;
        w = ww(i);
        xi2 = w*(fce - w) / (wp^2);    % (4) -> Omega_e = 1
        delta2 = 1 / (1+ xi2);  % (6) 
        delta = sqrt (delta2);
        xi = sqrt (xi2);
        vp = delta*xi*cv;  % (7)
        vg = (cv*xi/delta)/(xi^2 + fce/(2*(fce-w)));   % (8)
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
        s0 = delta*vperp/(xi*cv);    % (38)
        ss0(i) = s0;
        s1 = gamma * (1 - vr/vg)^2;     % (39)
        ss1(i) = s1;
        s_tmp = gamma*w*((vperp/cv)^2)/fce - ( 2 + (w/fce)*delta^2*(fce - gamma*w)/(fce -w) )*vr*vp/cv^2; 
        s2 = s_tmp/(2*xi*delta);   % (40)
        ss2(i) = s2;
        
        % optimum
        bw_eps_1 = 0.8*pi^(-5/2)*Q*vp_s*vg_s*uperp_s*wph_s^2 / (tau(1)*w_s*utpara_s);
        bw_eps_2 = (1 - vr_s/vg_s)^2 * exp(-0.5 * (gamma*vr_s/utpara_s)^2);
        bw_eps(i, j) = bw_eps_1 * bw_eps_2;  % (97)
         
        % threshold
        bwth_eps_1 = 100*pi^3*gamma^4*xi*(a_s*s2*utpara_s/Q)^2 / (w_s*wph_s^4*(delta*uperp_s)^5);
        bwth_eps(i, j) = bwth_eps_1 * exp((gamma*vr_s/utpara_s)^2); % (105)
    end
end

ind = round(fre / dw);
if (thop == 1)
    amp = bw_eps(ind,j);
else
    amp = bwth_eps(ind,j);
end

end