
% hole init
data_f = load('./mat/hole_init_0124.mat');
data = data_f.hist_data_perm;
ds = size(data);

mw = ds(3);
data = movmean(data,mw,3);
data_ave = sum(data,3) / ds(3);

save('./mat/hole_init_0124_2_ave','data_ave')