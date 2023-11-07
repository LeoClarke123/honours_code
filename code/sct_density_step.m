N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
K = 1; % spring const.
a = zeros(1, N0) + a0; % equilibrium length
d_vals = [0.8,0.9,1.25,1.5].*l0; % minimum length for division
dp = 0.01; % division probability
xp = dp/10; % death probability

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 20;
dt = 0.01;
total_cells = zeros(1, tmax/dt);
max_sim = 50;

for d=d_vals
    sim_num = 1;
    avg_length = zeros(1, round(tmax/dt,0));
    while sim_num<=max_sim
        k = zeros(1, N0) + K;
        t = 0;
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
            avg_length(c) = avg_length(c) + L/(length(x)-1);
            c=c+1;
            t = t+dt;
        end
        sim_num = sim_num+1;
    end
    plot(linspace(0,tmax,length(avg_length)), avg_length./max_sim, LineWidth=1.5)
    hold on
end

plot([0, tmax], [1, 1], 'r--', LineWidth=1.5)
xlabel('$$t$$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$$\langle l_{i} \rangle$$', 'Interpreter','latex', 'FontSize', 20)
grid on
grid minor
legend('$$l^{*}=0.8$$','$$l^{*}=0.9$$','$$l^{*}=1.25$$','$$l^{*}=1.5$$', 'Interpreter','latex', 'FontSize', 15)
ylim([0.5,1.75])
