function Gn = Nonlinear_growthrate(fld, fre)
% parameters
Q   = 0.1;
wpe = 15; % [\Omega_{e0}]
wph = 0.3; % [\Omega_{e0}]
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

vg = (cv .* sqrt(xi2) ./ sqrt(delta2)) ./ (xi2 + Omega_e ./ (2.0 .* (Omega_e - fre)));

term1 = (Q * wph^2 * vg) / (2 * gamma * utpara);

xi2 = abs(fre .* (Omega_e - fre) ./ wpe^2);
term2 = sqrt(abs(sqrt(xi2) ./ (fre .* fld)));

delta2 = 1 ./ (1 + xi2);
term3 = (sqrt(delta2) .* uperph ./ (pi * cv)).^1.5; 

w_s = fre ./ Omega_e;
vp_s = sqrt(xi2) .* sqrt(delta2);
vperp_s = vperp / cv;
vr_s = (w_s.^2 - sqrt(w_s.^4 + (w_s.^2 + vp_s.^2) .* (1 - w_s.^2 - vperp_s^2))) ./ (w_s.^2 ./ vp_s + vp_s);
vr = abs(vr_s .* cv);
term4 = exp(-1 * gamma^2 .* vr.^2 / (2 * utperp^2));

Gn = term1 .* term2 .* term3 .* term4;