function [s0, s1, s2] = inhomogeneity(fre)

% parameters
wpe = 15; % [\Omega_{e0}]
cv  = 1;
utpara = 0.26; % [c]
utperp = 0.30; % [c]
beta   = 0.3;
uperph  = sqrt(pi/2) * ((1 - beta^(1.5) / (1 - beta))) * utperp;
vperp  = cv / sqrt(cv^2 + (utpara^2 + uperph^2)) * uperph; % [c]
vpara  = cv / sqrt(cv^2 + (utpara^2 + uperph^2)) * utpara; % [c]
gamma  = 1 / sqrt(1 - (vperp^2 + vpara^2)/cv^2);
Omega_e = 1.00; % \Omega_{e0}

xi2 = abs(fre) .* (Omega_e - abs(fre)) ./ wpe^2;
xi2 = abs(xi2);
delta2 = 1 ./ (1 + xi2);

w_s = fre ./ Omega_e;

vp_s = sqrt(xi2) .* sqrt(delta2);
vperp_s = vperp / cv;
vr_s = (w_s.^2 - sqrt(abs(w_s.^4 + (w_s.^2 + vp_s.^2) .* (1 - w_s.^2 - vperp_s.^2)))) ./ (w_s.^2 ./ vp_s + vp_s);
vp = vp_s * cv;
vr = vr_s .* cv;
vr2 = (1 - Omega_e ./ (gamma .* fre)) .* vp;

vg = (cv .* sqrt(xi2) ./ sqrt(delta2)) ./ (xi2 + Omega_e ./ (2.0 .* (Omega_e - fre)));

% calculation
s0 = sqrt(delta2) .* vperp ./ (sqrt(xi2) .* cv);

s1 = gamma .* (1 - vr ./ vg).^2;

s2_term1 = 1 ./ (2 .* sqrt(xi2) .* sqrt(delta2));
s2_term2 = gamma .* fre .* (vperp/cv)^2 ./ Omega_e;
s2_term3 = (2 + (fre ./ Omega_e) .* delta2 .* (Omega_e - gamma .* fre) ./ (Omega_e - fre)) .* vr .* vp ./ cv^2;
s2 = s2_term1 .* (s2_term2 - s2_term3);