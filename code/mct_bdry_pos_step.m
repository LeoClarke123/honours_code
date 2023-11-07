N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
tmax = 500;
dt = 0.01;

cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
d_vals = 0.8*l0; % l_0
dp = dt; % division probability
xp = dp/10 ; % death probability

K1 = 20; K2_vals=[19]; 

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tic
times = zeros(1,length(K2_vals)*length(d_vals));


for d = d_vals

for K2=K2_vals
    sim_num = 0; max_sim=30;
    mean_bdry = zeros(1,round(tmax/dt)+1);
    while sim_num<max_sim
    
        k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % spring const.
        a = zeros(1, N0) + a0; % equilibrium length
        x = zeros(1, N0+1); % cell bondaries
        
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
                        k = remove(k,z);
                        a = remove(a,z);
                    end
                end
                z = z+1;
            end
            t = t+dt;

            flag = 1; v = 1;
            while flag == 1
                if v == length(k) && k(v) == K1
                    mean_bdry(c) = mean_bdry(c) + x(v+1); % K1 invades K2
                    v = 0;
                    flag = 0;
                elseif v == length(k) && k(v) == K2
                    mean_bdry(c) = mean_bdry(c) + 0; % K2 invades K1
                    v = 0;
                    flag = 0;
                elseif round((k(v)-k(v+1))*10)/10 ~= 0
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
    
    boundary = mean_bdry/max_sim;
    plot(linspace(0,tmax,length(mean_bdry)), boundary,'-', LineWidth=1.5)
    hold on
end

end

ylabel('$$s$$', 'Interpreter','latex', 'FontSize', 20)
xlabel('$$t$$', 'Interpreter','latex', 'FontSize', 20)
legend("$$K_{1}=2$$","$$K_{1}=5$$","$$K_{1}=10$$","$$K_{1}=20$$", 'Interpreter','latex', 'FontSize', 12)
grid on
grid minor
% ylim([25,50])
toc
