bt_case1 = './csv/bt_case1.csv';
bt_case2 = './csv/bt_case2.csv';
bt_case3 = './csv/bt_case3.csv';
bt_t = './csv/bt_t.csv';

bt1 = readmatrix(bt_case1);
bt2 = readmatrix(bt_case2);
bt3 = readmatrix(bt_case3);
t = readmatrix(bt_t);

figure,
plot(t,bt1, 'red')
hold on
plot(t,bt2, 'blue')
plot(t,bt3, 'green')
legend('Case 1', 'Case 2', 'Case 3')
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('B_w [B_0]')