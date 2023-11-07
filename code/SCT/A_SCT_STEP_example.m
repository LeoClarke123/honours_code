rng(4)

N0 = 50; % initial cell numbers
dt = 0.01; % timestep 
iter = 10000; % iterations
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
x = zeros(1, N0+1); % cedll bondaries
K=5;
k = zeros(1, N0) + K; % spring const.
a = zeros(1, N0) + a0; % equilibrium length
d = l0*1.5; % minimum length for division
dp = dt; % division probability
xp = dp/10; % death probability

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

    if mod(round(t,1), 0.1) == 0
        plot(x, zeros(1, length(x))+t, '.b', 'MarkerSize', 5)
        xlabel('$$x$$', 'Interpreter','latex', 'FontSize', 20)
        ylabel('$$t$$', 'Interpreter','latex', 'FontSize', 20)
        hold on
    end
    t = t+dt;
end

hold off
figure
plot(linspace(0,tmax, length(total_cells)), L./total_cells, 'b')
hold on
plot([0, tmax], [d, d], 'r--')
xlabel('$$t$$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$$\langle l_{i} \rangle$$', 'Interpreter','latex', 'FontSize', 20)
ylim([0.75,2])
grid on

