N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d_vals = linspace(0.5,3,20); % scale factors for a0 and d
dp = 0.5; % division probability
xp = 1e-6; % death probability
K = 1e-8;

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 2000;
dt = 0.01;

for d_val=d_vals
    d=d_val*a0; % minimum length for division

    sim_num = 1; max_sim = 20;
    total_cells = 0;

    while sim_num<=max_sim

        c=1;
        k = zeros(1, N0) + K; % spring const.
        x = zeros(1, N0+1); % cell bondaries
        a = zeros(1, N0) + d; % equilibrium length
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

            k = k(1:length(x)-1);
            a = a(1:length(x)-1);

            t = t+dt;
            if t>=(tmax-round(tmax/4)*dt)
                mean_cells(c) = length(x)-1;
                c=c+1;
            end   
        end

        total_cells = total_cells + sum(mean_cells)/length(mean_cells);
        sim_num = sim_num+1;
        
    end
    plot(d_val.^-1, L./(total_cells./max_sim),'ob');
    hold on
end

% plot([0,1],[a0,a0],'--k')
% hold on
% plot([1,1],[0,a0],'--k')

plot(d_vals.^-1, d_vals.*a0,'--r')
hold on
plot(d_vals.^-1, d_vals.*a0.*3./4,'--k')
hold on

ylabel('Equilibrium cell length', 'Interpreter','latex')
xlabel('$$1/l_{0}$$', 'Interpreter','latex')
title('Step Function: $$N_{0}=50, b=0.01, d=10^{-3}$$', 'Interpreter','latex')
legend('$$k=0.01$$','','','','','','','','','','','','','','','','','','','', ...
       '$$k=10$$','','','','','','','','','','','','','','','','','','','', ...
       '$$l^{*}$$', 'Interpreter','latex')
grid on
grid minor

