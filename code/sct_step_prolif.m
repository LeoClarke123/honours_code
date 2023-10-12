rng(0)

N0 = 10; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
x = zeros(1, N0+1); % cedll bondaries
k = zeros(1, N0) + 1; % spring const.
a = zeros(1, N0) + a0; % equilibrium length
d = l0; % minimum length for division
dp = 0.03; % division probability
xp = 5e-3; % death probability

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

i=0;
while i<length(x)
    x(i+1) = i*l0;
    i = i+1;
end

t = 0; tmax = 2;
dt = 0.05;
total_cells = zeros(1, tmax/dt);
while t < tmax
    total_cells(round((t/dt)+1,0)) = length(x)-1;
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

    if mod(round(t,1), 0.1) == 0
        plot(x, zeros(1, length(x))+tmax-t, '|k')
        xlabel('x')
        ylabel('time left in simulation')
        hold on
    end
    t = t+dt;
end

hold off
figure
plot(linspace(0,tmax, length(total_cells)), total_cells, 'k')
hold on
plot([0, tmax], [N0, N0], 'g')
ylim([N0*0.5, N0*1.5])


