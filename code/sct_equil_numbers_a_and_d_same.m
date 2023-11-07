N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d_vals = [1]; % scale factors for a0 and d
cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
dp = 0.01; % division probability
xp = 1e-3; % death probability
k_vals = linspace(1e-5,10,20);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

tmax = 1000;
dt = 0.01;

plot(k_vals, zeros(1,length(k_vals)) + a0,'--r')
hold on

col_ind = 1;

for d_val=d_vals
    d=d_val*a0; % minimum length for division

    sim_num = 1; max_sim = 10;
    total_cells = zeros(1, length(k_vals));

    while sim_num<=max_sim
        p=1;
        for K = k_vals
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
            total_cells(p) = total_cells(p) + sum(mean_cells)/length(mean_cells);
            p=p+1;
        end
        sim_num = sim_num+1;
    end
    plot1 = plot(k_vals, L./(total_cells./max_sim),'o-');
    hold on
    plot2 = plot(k_vals, zeros(1,length(k_vals))+d,'-');
    hold on
    set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
    col_ind = col_ind + 3;
end
hold on
spec = @(x) 2.^((N0-x.*(N0+L))./(2.*N0.*(x-1)))./(log(2))-((N0-x.*(N0+L))./(2.*N0.*(x-1)))-1/log(2);
plot(spec(linspace(0.5,0.8853,50)),linspace(0.5,0.8853,50),'--k')
hold on
spec2 = @(x) 2.^((x.*L)./(2.*(L-x.*N0)))./(log(2))-1/log(2);
plot(spec2(linspace(0.5,0.9,50)),linspace(0.5,0.9,50),'--r')

ylabel('Equilibrium cell length', 'Interpreter','latex')
xlabel('$$k$$', 'Interpreter','latex')
title('Step Function: $$N_{0}=50,\, b=0.01,\, d=10^{-3}$$', 'Interpreter','latex')
legend('','$$l_{0}=0.8a_{0}$$', '', '$$l_{0}=0.9a_{0}$$', '', '$$l_{0}=1.1a_{0}$$', '', '$$l_{0}=1.2a_{0}$$', 'Interpreter','latex')
ylim([0.6, 1.3])
grid on
grid minor

