N0 = 100; % initial cell numbers, MUST BE EVEN
dt = 0.01; % timestep 
iter = 50000; % iterations
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
x = zeros(1, N0+1); % cell bondaries
K1 = 10; K2 = 1; % spring const for soft/hard tissue
k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % vector of spring const
a = zeros(1, N0) + a0; % equilibrium length
d = l0*1.2; % half occupation constant
dp = dt; % max division probability
xp = dp/10; % death probability
tic
% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

i=0;
while i<length(x)
    x(i+1) = i*l0;
    i = i+1;
end

t = 0; tmax = iter*dt;
total_cells = zeros(1, iter);
while t < tmax
    
    % simulate ODE
    i=2;
    while i<length(x)
        x(i) = x(i)+dt*(k(i)*(x(i+1)-x(i)-a(i))-k(i-1)*(x(i)-x(i-1)-a(i-1)));
        if x(i) > 100 
            x(i) = 100;
        elseif x(i) < 0
            x(i) = 0;
        end
        i = i+1;
    end
    
    % cell division (step function)
    j=1;
    while j<length(x)
        if x(j+1)-x(j)>=d && rand(1) <= dp
            x = insert((x(j+1)+x(j))/2, x, j);
            k = insert(k(j), k, j);
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
                k = remove(k, 1);
                a = remove(a, 1);
            elseif z == length(x)-1
                x = remove(x,z);
                k = remove(k, z);
                a = remove(a, z);
            else
                x(z+1) = (x(z+1)+x(z))/2;
                x = remove(x,z);
                k = remove(k, z);
                a = remove(a, z);
            end
        end
        z = z+1;
    end
    
    if mod(round(t,2), 10*dt) == 0
        barrier = find(circshift(k,-1)-k == (K2-K1)); % K1 > K2
        if barrier > 0
            plot(x(1:barrier), zeros(1, barrier) + t, '.r','MarkerSize', 5)
            hold on
            plot(x(barrier+1:end), zeros(1, length(x(barrier+1:end))) + t, '.b','MarkerSize', 5)
            t_lim = t;
        end
    end
    t = t+dt;
end
toc
ylim([0,t_lim])
xlim([0,100])
title('k_{1} = 10, k_{2} = 1')
xlabel('x')
ylabel('t')



