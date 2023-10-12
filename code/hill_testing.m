L = linspace(0, 1, 1000);
N = [1, 2];
K = 0.5; K1 = 0.4; K2 = 0.6;
dt = 0.01;

d_vals_1 = [0.5*a0, 0.75*a0, 1.5*a0, 2*a0]; % unscaled l_0
xp = dt/50; % death probability
dp = xp*5; % max division probability 

for n = N
    d_vals = d_vals_1.*((dp-xp)./xp).^(1./n); % l_0
    P = dp*L.^n./(d_vals(1).^n+L.^n);
    plot(L, P)
    hold on
end

% plot([0, K1], [0, 0], '--r', [K1, K1], [0, max_prob], '--r', ...
%     [K1, 1], [max_prob, max_prob], '--r')
% hold on
% plot([0, K2], [0, 0], '--b', [K2, K2], [0, max_prob], '--b', ...
%     [K2, 1], [max_prob, max_prob], '--b')
% hold on
%plot([0, K], [0, 0], 'k', [K, K], [0, max_prob], 'k', ...
%    [K, 1], [max_prob, max_prob], 'k')
plot([0, 1], [xp, xp], 'k')
hold on
plot([0, 1], [0, 2*xp], 'b--')

legend('n=1', 'n=2', '', 'linear', 'Interpreter','latex', 'FontSize', 15)
xlabel('$$l$$', 'Interpreter','latex', 'FontSize', 15)
ylabel('$$P(l)$$', 'Interpreter','latex', 'FontSize', 15)
%ylim([-max_prob/10, max_prob + max_prob/10])
grid on
grid minor