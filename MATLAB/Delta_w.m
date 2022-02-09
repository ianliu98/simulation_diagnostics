function delta_w=Delta_w(fre,fld)

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
chi2 = 1 ./ (1 + xi2);
xi = sqrt(xi2);
chi = sqrt(chi2);

k   = fre ./ (cv .* chi .* xi);
wt  = sqrt(abs(k .* fld * vperp));
wtr = wt .* chi / sqrt(abs(gamma));
delta_w = 4.0 * wtr ./ (1.0 + chi2 .* (Omega_e ./ (gamma * fre) - 1.0) .* (xi2 + Omega_e ./ (2.0 * (Omega_e - fre)))); 

end