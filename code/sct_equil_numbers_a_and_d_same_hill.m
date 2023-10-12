N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
d_vals = [0.8, 0.9, 1.1, 1.2]; % scale factors for a0 and d
cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
dp = 0.01; % division probability
xp = 1e-3; % death probability
k_vals = linspace(1,15,10);
tic
% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 1; % smaller values smooth out Hill function

tmax = 2000;
dt = 0.05;

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
                
                % cell division (Hill function)
                rand_vals = rand(1, length(x)-1);
                spaces = circshift(x,-1)-x; spaces = spaces(1,1:end-1);
                hills = Hill(d,n,dp,spaces);
                buffer = 0;
                
                for j=1:length(hills)
                     if rand_vals(j)<=hills(j)
                         x = insert((x(j+1)+x(j))/2, x, j+buffer);
                         k = insert(K, k, j);
                         a = insert(a0, a, j);
                         buffer=buffer+1;
                     end
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
    plot2 = plot(k_vals, zeros(1,length(k_vals))+((xp.*d.^n)./(dp-xp)).^(1/n),'-');
    hold on
    set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
    col_ind = col_ind + 3;
end

ylabel('Equilibrium cell length', 'Interpreter','latex')
xlabel('$$k$$', 'Interpreter','latex')
title('Hill Function: $$n=1,\, N_{0}=50,\, b=0.01,\, d=10^{-3}$$', 'Interpreter','latex')
legend('','$$l_{0}=0.8a_{0}$$', '', '$$l_{0}=0.9a_{0}$$', '', '$$l_{0}=1.1a_{0}$$', '', '$$l_{0}=1.2a_{0}$$', 'Interpreter','latex')
ylim([0.05, 0.2])
grid on
grid minor
toc
