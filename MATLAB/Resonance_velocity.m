function vr = Resonance_velocity(fre)
% parameters
wpe = 15; % [\Omega_{e0}]
cv  = 1;
utpara = 0.26; % [c]
utperp = 0.30; % [c]
beta   = 0.3;
uperph  = sqrt(pi/2) * ((1 - beta^(1.5) / (1 - beta))) * utperp;
vperp  = cv / sqrt(cv^2 + (utpara^2 + uperph^2)) * uperph; % [c]
Omega_e = 1.00; % \Omega_{e0}
xi2 = abs(fre .* (Omega_e - fre) ./ wpe^2);
delta2 = 1 ./ (1 + xi2);
w_s = fre ./ Omega_e;
vp_s = sqrt(xi2) .* sqrt(delta2);
vperp_s = vperp / cv;
vr_s = (w_s.^2 - sqrt(w_s.^4 + (w_s.^2 + vp_s.^2) .* (1 - w_s.^2 - vperp_s^2))) ./ (w_s.^2 ./ vp_s + vp_s);
vr = abs(vr_s .* cv);