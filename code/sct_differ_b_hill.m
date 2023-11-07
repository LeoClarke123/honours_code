N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0; 
dt = 0.01;
d_vals = linspace(0.05,2,20); % scale factors for a0 and d
cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
dp_vals = [dt/5, dt/4, dt/2, dt]; % division probability

xp = dt/10; % death probability
K = 1;
tic
% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 10; % smaller values smooth out Hill function

tmax = 20;
col_ind = 1;

for dp = dp_vals
    plotter = zeros(1,length(d_vals));
    plot_index = 1;
for d_val=d_vals
    d=d_val*a0; % minimum length for division

    sim_num = 1; max_sim = 10;
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

        total_cells = total_cells + sum(mean_cells)/length(mean_cells);
        sim_num = sim_num+1;
        
    end
    plotter(plot_index) = L./(total_cells./max_sim);
    plot_index = plot_index+1;
end
plot1 = plot(d_vals, plotter,'.-', MarkerSize=20);
set(plot1,'Color', cols(col_ind:col_ind+2))
hold on
col_ind = col_ind + 3;
end
% plot([0,1],[a0,a0],'--k')
% hold on
% plot([1,1],[0,a0],'--k')

% plot(d_vals, d_vals.*a0,'--k')
% hold on
toc
ylabel('$$\langle l_{i} \rangle$$', 'Interpreter','latex','FontSize', 15)
xlabel('$$l_{H}$$', 'Interpreter','latex','FontSize', 15)
legend('$$b_{H}=0.002$$','$$b_{H}=0.0025$$','$$b_{H}=0.005$$','$$b_{H}=0.01$$', 'Interpreter','latex' ,'FontSize', 15)
grid on
grid minor

