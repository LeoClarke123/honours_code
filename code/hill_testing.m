L = linspace(0, 2, 100);
N = [1]; N2 = [10];
K = 0.5; K1 = 0.4; K2 = 0.6;
dt = 0.01;

%d_vals_1 = [0.5*a0, 0.75*a0, 1.5*a0, 2*a0]; % unscaled l_0
dp = dt; % max division probability 
xp = dp/10; % death probability
% 
% for n = N
%     d_vals = d_vals_1.*((dp-xp)./xp).^(1./n); % l_0
%     P = dp*L.^n./(d_vals(1).^n+L.^n);
%     plot(L, P, LineWidth=1.5)
%     hold on
% end
% 
% for n = N2
%     P = dp/2.5*L.^n./(d_vals_1(1).^n+L.^n);
%     plot(L, P, LineWidth=1.5)
%     hold on
% end

% plot([0, K1], [0, 0], '--r', [K1, K1], [0, max_prob], '--r', ...
%     [K1, 1], [max_prob, max_prob], '--r')
% hold on
% plot([0, K2], [0, 0], '--b', [K2, K2], [0, max_prob], '--b', ...
%     [K2, 1], [max_prob, max_prob], '--b')
% hold on
%plot([0, K], [0, 0], 'k', [K, K], [0, max_prob], 'k', ...
%    [K, 1], [max_prob, max_prob], 'k')

plot([0, 2], [0, 1.5*xp],LineWidth=1.5)
hold on
plot(2/1.5, xp,'*r',MarkerSize=15)
hold on
plot([0, 2], [0, 2*xp],LineWidth=1.5)
hold on
plot(2/2, xp,'*r',MarkerSize=15)
hold on
plot([0, 2], [0, 3*xp],LineWidth=1.5)
hold on
plot(2/3, xp,'*r',MarkerSize=15)
hold on
% plot([0, 0.5], [0, 0], 'r',LineWidth=1.5)
% hold on
% plot([0.5, 0.5], [0, xp*2], 'r',LineWidth=1.5)
% hold on
% plot([0.5, 1], [xp*2, xp*2], 'r',LineWidth=1.5)
% hold on
plot([0, 2], [xp, xp], 'k--', LineWidth=1.5)


legend('Linear (Stretched)','', 'Linear (Equal)','', 'Linear (Compressed)','','Death', 'latex', 'FontSize', 15)
xlabel('$$l$$', 'Interpreter','latex', 'FontSize', 15)
ylabel('$$P(l)$$', 'Interpreter','latex', 'FontSize', 15)
%ylim([-max_prob/10, max_prob + max_prob/10])
grid on
grid minor