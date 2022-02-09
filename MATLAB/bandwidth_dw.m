% delta-w: bandwidth of coherent waves

clc
clear

bbound = 1.0;
thop = 1.0;
wpes = 15;

Dw = zeros(length(wpes),9000);
for i=1:length(wpes)
    [delta_w, ww]=separability(bbound,thop, wpes(i));
    Dw(i,:) = delta_w;
end

figure,
loglog(ww,Dw,'black')
xlabel('\omega [\Omega_{e0}]')
ylabel('Bandwith (\Delta\omega)')
axis([5e-4 1e-0 1e-6 2e-2])
xline(0.05,'r--')