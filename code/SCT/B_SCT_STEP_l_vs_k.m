N0 = 100; % initial cell numbers
dt = 0.005; % timestep
iter = 1e5; % iterations
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d = l0; % minimum length for division
dp = dt; % division probability
xp = dp/10; % death probability
k_vals = linspace(0,20,20);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = iter*dt;
total_cells = zeros(1, length(k_vals));
p=1;
tic
for K = k_vals
    c=1;
    k = zeros(1, N0) + K; % spring const.
    a = zeros(1, N0) + a0; % equilibrium length
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
legend('N0=100, b=5e-4, d=5e-5')
grid on
grid minor

