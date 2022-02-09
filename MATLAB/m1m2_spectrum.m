wtf = readmatrix('./csv/homo_wtf3.csv');
wtf_t = readmatrix('./csv/homo_wtf3_time.csv');
wtf_w = readmatrix('./csv/homo_wtf3_fre.csv');

wtfs = reshape(wtf,[length(wtf_t), length(wtf_w), 5]);
wtf_tt = wtf_t * 1e4;

ind = 4;
wtf_maxv = max(max(wtfs(:,:,ind)));



if (ind == 4)
% rising tone
    tl = 13144.1;
    th = 15114.2;
    wl = 0.045;
    wh = 0.047;
else
% falling tone
    tl = 10461.2;
    th = 13017.1;
    wl = 0.038;
    wh = 0.040;
end

tliv = find(wtf_tt < tl);
thiv = find(wtf_tt > th);
tli = tliv(end)+1;
thi = thiv(1) - 1;

wliv = find(wtf_w < wl);
whiv = find(wtf_w > wh);
wli = wliv(end);
whi = whiv(1);

% ---

figure,
colormap('jet')

subplot(1,2,1)
mesh(wtf_tt,wtf_w,transpose(wtfs(:,:,ind)),'EdgeColor','flat','FaceColor','flat')
view(2)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\omega [\Omega_{e0}]')
h = colorbar;
h.Title.String = 'log_{10}(B_w/B_0)';
%ylabel(h,'log_{10}(B_w/B_0)')
rectangle('Position',[wtf_tt(tli) wtf_w(wli) wtf_tt(thi)-wtf_tt(tli) wtf_w(whi)-wtf_w(wli)], 'LineWidth', 1)
axis tight;
caxis([wtf_maxv-1 wtf_maxv])
text(-0.1,1.05,'(a)','Units','normalized','FontSize',10)


subplot(1,2,2)
mesh(wtf_tt(tli:thi),wtf_w(wli:whi),transpose(wtfs(tli:thi,wli:whi,ind)),'EdgeColor','flat','FaceColor','flat')
view(2)
xlabel('t [\Omega_{e0}^{-1}]')
ylabel('\omega [\Omega_{e0}]')
h = colorbar;
%ylabel(h,'log_{10}(B_w/B_0)')
h.Title.String = 'log_{10}(B_w/B_0)';
axis tight;
caxis([wtf_maxv-1 wtf_maxv])
text(-0.1,1.05,'(b)','Units','normalized','FontSize',10)
