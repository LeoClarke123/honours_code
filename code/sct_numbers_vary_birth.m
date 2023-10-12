N0 = 50; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
a = zeros(1, N0) + a0; % equilibrium length
d = l0; % minimum length for division
dp_vals = linspace(0.01, 0.1, 20); % division probability
xp_vals = 5e-5; % death probability
k = zeros(1, N0) + 10; % spring const.

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

g = 50000; % number of vals to take mean of

tmax = 500;
dt = 0.001;
total_cells = zeros(1, length(dp_vals));
p=1;
for dp = dp_vals
    c=1;
    x = zeros(1, N0+1); % cell bondaries
    mean_cells = zeros(1,g);
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
        if t>=(tmax-g*dt)
            mean_cells(c) = length(x)-1;
            c=c+1;
        end   
    end
    total_cells(p) = sum(mean_cells)/length(mean_cells);
    p=p+1;
end

plot(dp_vals, total_cells,'-o')
hold on
ylim([50 75])
ylabel('Equilibrium cell numbers N')
xlabel('birth probability b')
legend("N0: " + N0 + ", d: 0.001, k: 4, dt: " + dt, ... 
      "N0: " + N0 + ", d: 0.001, k: 10, dt: " + dt, ...
      "N0: " + N0 + ", d: 5e-5, k: 4, dt: " + dt, ... 
      "N0: " + N0 + ", d: 5e-5, k: 10, dt: " + dt)
grid on
grid minor

