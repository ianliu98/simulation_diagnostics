clc
clear

ktf  = readmatrix('./csv/k_t_field.csv');
ktdk = readmatrix('./csv/k_t_deltak.csv');
ktdw = readmatrix('./csv/k_t_deltaw.csv');
ktt  = readmatrix('./csv/k_t_t.csv');
ktk  = readmatrix('./csv/k_t_k.csv');

figure,
colormap('jet')
subplot(1,2,1)
mesh(ktk,ktt,ktf,'FaceColor','interp')
view(2)
xlabel('k [c^{-1}\Omega_{e0}]')
ylabel('t [\Omega_{e0}^{-1}]')
axis tight;
h1 = colorbar();
ylabel(h1,'B_w [B_0]')

subplot(1,2,2)
mesh(ktk,ktt,ktdk,'FaceColor','interp')
view(2)
xlabel('k [c^{-1}\Omega_{e0}]')
ylabel('t [\Omega_{e0}^{-1}]')
axis tight;
h2 = colorbar();
ylabel(h2, '\Delta k')

ind1 = 2;
ind2 = 11;
ind3 = 18;

kind1 = ktk(ind1);
kind2 = ktk(ind2);
kind3 = ktk(ind3);

flds = zeros(3,length(ktt));
flds(1,:) = ktf(:,ind1);
flds(2,:) = ktf(:,ind2);
flds(3,:) = ktf(:,ind3);

dks = zeros(3,length(ktt));
dks(1,:) = ktdk(:,ind1);
dks(2,:) = ktdk(:,ind2);
dks(3,:) = ktdk(:,ind3);

figure,
subplot(2,1,1)
plot(ktt,flds)
subplot(2,1,2)
plot(ktt,dks)

knew1 = kind1 + max(dks(1,1:1000));
knew2 = kind2 + max(dks(2,1:1000));
knew3 = kind3 + max(dks(3,1:1000));


figure,
tt = 4000;
colormap('jet')
mesh(ktk,ktt(1:tt),ktf(1:tt,:))
view(2)
xlabel('k [c^{-1}\Omega_{e0}]')
ylabel('t [\Omega_{e0}^{-1}]')
title("t = "+string(ktt(tt))+" [\Omega_{e0}^{-1}]")
axis tight;
h1 = colorbar();
ylabel(h1,'B_w [B_0]')

figure,
tt = [500, 1000, 2000, 4000];
colormap('jet')
subplot(2,2,1)
mesh(ktk,ktt(1:tt(1)),ktf(1:tt(1),:))
view(2)
ylabel('t [\Omega_{e0}^{-1}]')
title("t = "+string(ktt(tt(1)))+" [\Omega_{e0}^{-1}]")
axis tight;
colorbar();
text(-0.1,1.1,'(a)','Units','normalized','FontSize',10)

subplot(2,2,2)
mesh(ktk,ktt(1:tt(2)),ktf(1:tt(2),:))
view(2)
title("t = "+string(ktt(tt(2)))+" [\Omega_{e0}^{-1}]")
axis tight;
h1 = colorbar();
ylabel(h1,'B_w [B_0]')
text(-0.1,1.1,'(b)','Units','normalized','FontSize',10)

subplot(2,2,3)
mesh(ktk,ktt(1:tt(3)),ktf(1:tt(3),:))
view(2)
xlabel('k [c^{-1}\Omega_{e0}]')
ylabel('t [\Omega_{e0}^{-1}]')
title("t = "+string(ktt(tt(3)))+" [\Omega_{e0}^{-1}]")
axis tight;
colorbar();
text(-0.1,1.1,'(c)','Units','normalized','FontSize',10)

subplot(2,2,4)
mesh(ktk,ktt(1:tt(4)),ktf(1:tt(4),:))
view(2)
xlabel('k [c^{-1}\Omega_{e0}]')
title("t = "+string(ktt(tt(4)))+" [\Omega_{e0}^{-1}]")
axis tight;
h1 = colorbar();
ylabel(h1,'B_w [B_0]')
text(-0.1,1.1,'(d)','Units','normalized','FontSize',10)