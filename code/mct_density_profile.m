N0 = 50; % initial cell numbers, MUST BE EVEN
dt = 0.01; % timestep 
iter = 50000; % iterations
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
K1 = 1.5; K2 = 1; % spring const for soft/hard tissue
d = l0*0.8; % half occupation constant
dp = dt; % max division probability
xp = dp/10; % death probability
tic
% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

prof_size = 15; % number of points either side in profile

full_density = zeros(1, prof_size*2 + 1);

tot_sim = 30;
sim = 1;

while sim<=tot_sim
    
    x = zeros(1, N0+1); % cell bondaries
    k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % vector of spring const
    a = zeros(1, N0) + a0; % equilibrium length
    dens_count = 0;
    temp_density = zeros(1, prof_size*2 + 1);
    
    i=0;
    while i<length(x)
        x(i+1) = i*l0;
        i = i+1;
    end
    
    t = 0; tmax = iter*dt;
    total_cells = zeros(1, iter);
    while t < tmax
        
        % simulate ODE
        i=2;
        while i<length(x)
            x(i) = x(i)+dt*(k(i)*(x(i+1)-x(i)-a(i))-k(i-1)*(x(i)-x(i-1)-a(i-1)));
            if x(i) > B 
                x(i) = B;
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
                    k = remove(k, z);
                    a = remove(a, z);
                end
            end
            z = z+1;
        end
        
        if mod(round(t,2), 10*dt) == 0
            barrier = find(circshift(k,-1)-k == (K2-K1)); % K1 > K2
            if barrier > 0
                if (barrier > (prof_size)) && (barrier < (length(x)-prof_size-1))
                    dens = (x(barrier-prof_size+1:barrier+prof_size+1)-x(barrier-prof_size:barrier+prof_size)).^-1;
                    if sum(abs(dens)>5) == 0
                        temp_density = temp_density + dens;
                        dens_count = dens_count + 1;
                    end
                end
            end
        end
        t = t+dt;
    end
    full_density = full_density + temp_density./dens_count;
    sim = sim + 1;
end
toc
full_density = full_density./tot_sim;
plot(linspace(-prof_size,prof_size,length(full_density)),full_density,'.', 'MarkerSize', 15)
xlabel('$$x$$', 'Interpreter','latex','FontSize', 20)
ylabel('$$q$$', 'Interpreter','latex','FontSize', 20)
%legend("$$K_{1}-K_{2}=0.5$$", 'Interpreter','latex', 'FontSize', 15)
legend("$$K_{1}-K_{2}=3$$","$$K_{1}-K_{2}=2$$","$$K_{1}-K_{2}=1$$","$$K_{1}-K_{2}=0.5$$", 'Interpreter','latex', 'FontSize', 15)
% ylim([1.3,1.7])
grid on
grid minor
hold on


