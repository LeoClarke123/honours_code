N0 = 50; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d_vals = [0.8*l0, 0.9*l0, 1.1*l0, 1.2*l0]; % minimum length for division
a = a0; % equlibrium length
dp = 0.03; % division probability
xp = 5e-5; % death probability
K1 = 10; K2=10; 

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 300; t_init = 50;
dt = 0.001;

for d1 = d_vals
    a1 = a0; a2 = a0;
    d2 = l0;
    sim_num = 0; max_sim=20;
    mean_bdry = zeros(1,round(tmax/dt)+1);
    while sim_num<max_sim
        
        a = [zeros(1, N0/2) + a1, zeros(1, N0/2) + a2]; % equilibrium length
        d = [zeros(1, N0/2) + d1, zeros(1, N0/2) + d2]; % equilibrium length
        k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % spring const.
        x = zeros(1, N0+1); % cell bondaries
        cell_bdry = zeros(1,round(tmax/dt));
        
        i=0;
        while i<length(x)
            x(i+1) = i*l0;
            i = i+1;
        end
       
        c=1;
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
                if d(j) == d1 && (x(j+1)-x(j))>=d1 && rand(1) <= dp
                    x = insert((x(j+1)+x(j))/2, x, j);
                    k = insert(k(j), k, j);
                    a = insert(a(j), a, j);
                    d = insert(d(j), d, j);
                elseif d(j) == d2 && (x(j+1)-x(j))>=d2 && rand(1) <= dp
                    x = insert((x(j+1)+x(j))/2, x, j);
                    k = insert(k(j), k, j);
                    a = insert(a(j), a, j);
                    d = insert(d(j), d, j);
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
                        d = remove(d, 1);
                    elseif z == length(x)-1
                        x = remove(x,z);
                        k = remove(k, z);
                        a = remove(a, z);
                        d = remove(d, z);
                    else
                        x(z+1) = (x(z+1)+x(z))/2;
                        x = remove(x, z);
                        k = remove(k, z);
                        a = remove(a, z);
                        d = remove(d, z);
                    end
                end
                z = z+1;
            end
            t = t+dt;
            flag = 1; v = 1;
            while flag == 1
                if v == length(a)
                    cell_bdry(c) = cell_bdry(c-1);
                    mean_bdry(c) = mean_bdry(c) + cell_bdry(c-1);
                    v = 0;
                    flag = 0;
                elseif d(v)-d(v+1) ~= 0
                    cell_bdry(c) = x(v+1);
                    mean_bdry(c) = mean_bdry(c) + x(v+1);
                    flag = 0;
                else
                    v = v+1;
                end
            end   
            c=c+1;
        end
        sim_num = sim_num+1;
    end
    plot(linspace(0,tmax,length(mean_bdry)), mean_bdry/max_sim,'-')
    hold on
end
ylabel('Cell boundry position S(t)')
xlabel('time t')
title('N0=50, k = 10, b=0.03, d=5e-5')
legend('l1 = 0.8*l2', 'l1 = 0.9*l2', 'l1 = 1.1*l2', 'l1 = 1.2*l2')
ylim([0,100])
grid on
grid minor

