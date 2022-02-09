function vg = Group_velocity(fre)
wpe = 15;
cv  = 100;
Omega_e = 1.00;
xi2 = abs(fre) .* (Omega_e - abs(fre)) ./ wpe^2;
xi2 = abs(xi2);
chi2 = 1 ./ (1 + xi2);
vg = (cv .* sqrt(xi2) ./ sqrt(chi2)) ./ (xi2 + Omega_e ./ (2.0 .* (Omega_e - fre)));
