N0 = 100; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
a = zeros(1, N0) + a0; % equilibrium length
d = l0; % minimum length for division
dp = 5e-4; % division probability
xp = 5e-5; % death probability
k_vals = linspace(0,20,20);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

Hill = @(d, n, dp, L) dp*(L^n/(d^n+L^n));
n = 50; % smaller values smooth out Hill function

tmax = 1000;
dt = 0.01;
total_cells = zeros(1, length(k_vals));
p=1;
tic
for K = k_vals
    c=1;
    k = zeros(1, N0) + K; % spring const.
    x = zeros(1, N0+1); % cell bondaries
    mean_cells = zeros(1,round(tmax/4));
    i=0;
    while i<length(x)
        x(i+1) = i*l0;
        i = i+1;
    end
    t = 0;
    while t < tmax
        % simulate ODE
        i=2;
        while i<length(x)
            x(i) = x(i)+dt*(k(i)*(x(i+1)-x(i)-a(i))-k(i-1)*(x(i)-x(i-1)-a(i-1)));
            i = i+1;
        end
        
        % cell division (step function)
        j=1;
        while j<length(x)
            % if rand(1) <= Hill(d, n, dp, x(j+1)-x(j))
            %if x(j+1)-x(j)>=d && rand(1) <= dp
            if x(j+1)-x(j)>=d && rand(1) <= dp
                x = insert((x(j+1)+x(j))/2, x, j);
                k = insert(K, k, j);
                a = insert(a0, a, j);
            end
            j = j+1;
        end
        
        % cell death
        z=1;
        while z<length(x)
            if rand(1) <= xp
                if z == 1
                    x = remove(x,2);
                elseif z == length(x)-1
                    x = remove(x,z);
                else
                    x(z+1) = (x(z+1)+x(z))/2;
                    x = remove(x,z);
                end
            end
            z = z+1;
        end
        t = t+dt;
        if t>=(tmax-round(tmax/4)*dt)
            mean_cells(c) = length(x)-1;
            c=c+1;
        end   
    end
    total_cells(p) = sum(mean_cells)/length(mean_cells);
    p=p+1;
end
toc
plot(k_vals, total_cells,'o')
ylabel('Equilibrium cell numbers N')
xlabel('spring constant k')
legend('N0=50, b=0.03, d=5e-5')
grid on
grid minor

