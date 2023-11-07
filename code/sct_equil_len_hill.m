N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
base_d_vals = linspace(0.05,2,20); % scale factors for a0 and d
cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
dp = 0.01; % division probability
xp = dp/10; % death probability
k_vals = [0.01, 1, 5, 10];

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 10; % smaller values smooth out Hill function

d_vals = base_d_vals.*((dp-xp)./xp).^(1./n); % length parameter 

tmax = 40;
dt = 0.01;
col_ind = 1;
tic
for K = k_vals
    plotter = zeros(1,length(d_vals));
    plot_index = 1;
for d=d_vals
 
    sim_num = 1; max_sim = 10;
    total_cells = 0;

    while sim_num<=max_sim

        c=1;
        x = zeros(1, N0+1); % cell bondaries
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
                x(i) = x(i)+dt*(K*(x(i+1)-x(i)-a0)-K*(x(i)-x(i-1)-a0));
                i = i+1;
            end
            
            % cell division (Hill function)
            rand_vals = rand(1, length(x)-1);
            spaces = circshift(x,-1)-x; spaces = spaces(1,1:end-1);
            hills = Hill(d,n,dp,spaces);
            
            for j=1:length(hills)
                 if rand_vals(j)<=hills(j)
                     x = insert((x(j+1)+x(j))/2, x, j);
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
plot1 = plot(base_d_vals, plotter,'.-', MarkerSize=20);
set(plot1,'Color', cols(col_ind:col_ind+2))
hold on
col_ind = col_ind + 3;
end
% plot([0,1],[a0,a0],'--k')
% hold on
% plot([1,1],[0,a0],'--k')

plot(base_d_vals, base_d_vals,'--k',LineWidth=1.3)
hold on
toc
ylabel('$$\langle l_{i} \rangle$$', 'Interpreter','latex','FontSize', 20)
xlabel('$$l^{*}$$', 'Interpreter','latex','FontSize', 20)
legend('$$k=0.01$$','$$k=1$$','$$k=5$$','$$k=10$$','$$l^{*}$$', 'Interpreter','latex' ,'FontSize', 15)
grid on
grid minor

