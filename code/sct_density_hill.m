N0 = 500; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
k_vals = [0,10,20]; % spring const.
a = zeros(1, N0) + a0; % equilibrium length
d= l0; % median length for division
dp = 0.05; % max division probability
xp = 5e-3; % death probability

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp*(L^n/(d^n+L^n));
n = 2; % smaller values smooth out Hill function

tmax = 10;
dt = 0.01;
total_cells = zeros(1, tmax/dt);

for K=k_vals
    k = zeros(1, N0) + K;
    t = 0;
    density = zeros(1, round(tmax/dt,0)+1);
    c=1;
    i=0;
    x = zeros(1, N0+1); % cell bondaries
    while i<length(x)
        x(i+1) = i*l0;
        i = i+1;
    end

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
            if rand(1) <= Hill(d, n, dp, x(j+1)-x(j))
                x = insert((x(j+1)+x(j))/2, x, j);
                k = insert(1, k, j);
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
        density(c) = (length(x)-1)/L;
        c=c+1;
        t = t+dt;
    end
    plot(linspace(0,tmax,length(density)), density)
    hold on
end

xlabel('t')
ylabel('density')
grid on
grid minor
legend("k=0","k=1","k=2","k=3","k=4","k=5","k=6","k=7","k=8","k=9","k=10")
title("max division prob: " + dp + ", death prob: " + xp + ", dt: " + dt + ", n: " + n)

