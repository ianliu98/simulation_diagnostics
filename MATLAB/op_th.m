clc
clear

data = readmatrix('./csv/bomg_case1.csv');
len = length(data(1,:));

opth = load('./mat/op_th.mat');
w  = opth.ww;
op = opth.opt;
th = opth.the;

w  = w(1:3000);
op = op(1:3000);
th = th(1:3000);

figure,
scatter(data(1,:),data(2,:),'.black')
set(gca, 'yscale', 'log')
hold on
plot(w, op, 'red')
plot(w, th, 'blue')
ylim([1e-7 1e-3])
xlabel('\omega [\Omega_{e0}]')
ylabel('B_w [B_0]')

