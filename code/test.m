N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
dp = 0.01; % division probability
xp = dp/10; % death probability

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 10; % smaller values smooth out Hill function

l_H = l0.*((dp-xp)./xp).^(1./n); % minimum length for division

l_vals = linspace(0,2,100);

plot(l_vals, Hill(l_H,n,dp,l_vals))
hold on
plot([0,2],[0,xp*2], '--r')
hold on
plot([0,2],[xp,xp], '--k')