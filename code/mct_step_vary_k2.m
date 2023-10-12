N0 = 100; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d = 1.5*a0; % minimum length for division
dp = 0.03; % division probability
xp = 5e-5; % death probability
K1 = 1; K2_vals=[2,5,10,15,20]; 

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 1000;
dt = 0.005;

for K2=K2_vals
    sim_num = 0; max_sim=30;
    mean_bdry = zeros(1,round(tmax/dt)+1);
    while sim_num<max_sim
    
        k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % spring const.
        a = zeros(1, N0) + a0; % equilibrium length
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
                        k = remove(k,z);
                        a = remove(a,z);
                    end
                end
                z = z+1;
            end
            t = t+dt;
            flag = 1; v = 1;
            while flag == 1
                if v == length(k)
                    cell_bdry(c) = x(v+1);
                    mean_bdry(c) = mean_bdry(c) + x(v+1);
                    v = 0;
                    flag = 0;
                elseif round((k(v)-k(v+1))) ~= 0
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
title('N0=100, k1 = 1, l_{0}=1.5*a_{0}')
legend("k2=2","k2=5","k2=10","k2=15","k2=20")
grid on
grid minor

