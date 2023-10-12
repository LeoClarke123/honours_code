N0 = 50; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = 3*L/N0; a0 = 3*L/N0;
a = zeros(1, N0) + a0; % equilibrium length
d = l0; % minimum length for division
dp = 5e-4; % division probability
xp = 5e-5; % death probability

k_vals = linspace(0,20,10);

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L1) dp.*(L1.^n./(d.^n+L1.^n));
n = 3; % smaller values smooth out Hill function

exp_equil = L/((l0^n)*xp/(dp-xp))^(1/n); % equilibrium value

tmax = 800;
dt = 0.005;

sim_num = 1; max_sim = 20;
total_cells = zeros(1, length(k_vals));
total_cells_std = zeros(1, length(k_vals));
tic
while sim_num<=max_sim
    p=1;
    for K = k_vals
        c=1;
        k = zeros(1, N0) + K; % spring const.
        x = zeros(1, N0+1); % cell bondaries
        mean_cells = zeros(1,round(tmax/4));
        i=0;
        while i<length(x)
            x(i+1) = i*(L/N0);
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
                     k = insert(K, k, j+buffer);
                     a = insert(a0, a, j+buffer);
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
            t = t+dt;
            if t>=(tmax-round(tmax/4))
                mean_cells(c) = length(x)-1;
                c=c+1;
            end   
        end
        total_cells(p) = total_cells(p) + mean(mean_cells);
        total_cells_std(p) = total_cells_std(p) + std(mean_cells);
        p=p+1;
    end
    sim_num = sim_num+1;
end
toc
errorbar(k_vals, total_cells./max_sim, (total_cells_std./max_sim)./sqrt(max_sim), 'o-')
hold on
plot([0,20],[exp_equil, exp_equil], 'k--')
ylabel('Equilibrium cell numbers N')
xlabel('spring constant k')
title('N_{0}=50, b=5e-4, d=5e-5, n=2, L=100, l_{0} = 3*L/N_{0}')
legend('l_{0}=a_{0}', 'Theoretical Equilibrium')
grid on
grid minor

