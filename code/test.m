N0 = 50; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
a = zeros(1, N0) + a0; % equilibrium length
d = l0; % minimum length for division
dp = 5e-4; % division probability
xp = 5e-5; % death probability
k_vals = linspace(0,20,10);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 2; % smaller values smooth out Hill function

x = zeros(1, N0+1); % cell bondaries
i=0;
while i<length(x)
    x(i+1) = i*l0;
    i = i+1;
end
% tic
% for i=1:1000
%     x = insert(i,x,5);
% end
% toc

tic
for i=1:1000
    x = [x(1:4) i x(5:end)];
end
toc
            