N0 = 50; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0;
a0_vals = [0.5, 0.8, 1, 1.2, 1.5];
d = l0; % minimum length for division
dp = 0.03; % division probability
xp = 5e-5; % death probability
k_vals = linspace(0,20,20);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 1000;
dt = 0.01;

for a0_val=a0_vals
    a0=a0_val*l0;

    sim_num = 1; max_sim = 20;
    total_cells = zeros(1, length(k_vals));

    while sim_num<=max_sim
        p=1;
        for K = k_vals
            c=1;
            k = zeros(1, N0) + K; % spring const.
            x = zeros(1, N0+1); % cell bondaries
            a = zeros(1, N0) + a0; % equilibrium length
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
                if t>=(tmax-round(tmax/4))
                    mean_cells(c) = length(x)-1;
                    c=c+1;
                end   
            end
            total_cells(p) = total_cells(p) + sum(mean_cells)/length(mean_cells);
            p=p+1;
        end
        sim_num = sim_num+1;
    end
    plot(k_vals, L./(total_cells./max_sim),'o-')
    hold on
end
ylabel('Equilibrium cell numbers N')
xlabel('spring constant k')
title('N0=50, b=0.03, d=5e-5')
legend('0.5*l_{0}=a_{0}', '0.8*l_{0}=a_{0}', 'l_{0}=a_{0}', '1.2*l_{0}=a_{0}', '1.5*l_{0}=a_{0}')
grid on
grid minor

