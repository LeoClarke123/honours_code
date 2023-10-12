N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
a = zeros(1, N0) + a0; % equilibrium length
d = a0*1.5; % minimum length for division
dp = 0.03; % division probability
xp = 5e-4; % death probability
K1 = 1; K2=2; 

% insert/remove new cell upon div
insert = @(m, x, n)cat(2,  x(1:n), m, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% hill function for cell division prob
Hill = @(d, n, dp, L) dp.*(L.^n./(d.^n+L.^n));
n = 1; % smaller values smooth out Hill function

tmax = 150;
dt = 0.005;
tic
sim_num = 0; max_sim=20;
mean_bdry = zeros(1,round(tmax/dt)+1);
while sim_num<max_sim

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
        
        % cell division (Hill function)
        rand_vals = rand(1, length(x)-1);
        spaces = circshift(x,-1)-x; spaces = spaces(1,1:end-1);
        hills = Hill(d,n,dp,spaces);
        buffer = 0;
        
        for j=1:length(hills)
             if rand_vals(j)<=hills(j)
                 x = insert((x(j+1)+x(j))/2, x, j+buffer);
                 k = insert(k(j), k, j);
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
                cell_bdry(c) = cell_bdry(c-1);
                mean_bdry(c) = mean_bdry(c) + cell_bdry(c-1);
                v = 0;
                flag = 0;
            elseif k(v)-k(v+1) ~= 0
                cell_bdry(c) = v;
                mean_bdry(c) = mean_bdry(c) + (x(v+1)+x(v))/2;
                flag = 0;
            else
                v = v+1;
            end
        end   
        c=c+1;
    end
    sim_num = sim_num+1;
end
toc
plot(linspace(0,tmax,length(mean_bdry)), mean_bdry/max_sim,'-')
ylabel('Cell boundry position S(t)')
xlabel('time t')
title('N0=50, k1 = 10, k2 = 15')
ylim([0,50])
grid on
grid minor

