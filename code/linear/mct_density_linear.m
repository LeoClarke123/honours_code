N0 = 50; % initial cell numbers, MUST BE EVEN
dt = 0.01; % timestep 
iter = 50000; % iterations
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
K1 = 1.5; K2 = 1; % spring const for soft/hard tissue
d = l0*1.25; % half occupation constant
xp = dt/10; % death probability
tic

% linear function for cell division prob
Lin = @(xp, d, L) xp.*L./d;

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

prof_size = 20; % number of points either side in profile

full_density = zeros(1, prof_size*2 + 1);

tot_sim = 50;
sim = 1;

while sim<=tot_sim
    
    x = zeros(1, N0+1); % cell bondaries
    k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % vector of spring const
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
            x(i) = x(i)+dt*(k(i)*(x(i+1)-x(i)-a0)-k(i-1)*(x(i)-x(i-1)-a0));
            if x(i) > B 
                x(i) = B;
            elseif x(i) < 0
                x(i) = 0;
            end
            i = i+1;
        end
        
        % cell division (Linear function)
        rand_vals = rand(1, length(x)-1);
        spaces = circshift(x,-1)-x; spaces = spaces(1,1:end-1);
        lins = Lin(xp,d,spaces);
        
        for j=1:length(lins)
             if rand_vals(j)<=lins(j)
                 x = insert((x(j+1)+x(j))/2, x, j);
                 k = insert(k(j), k, j);
             end
        end
        
        % cell death
        z=1;
        while z<length(x)
            if rand(1) <= xp
                if z == 1
                    x = remove(x,2);
                    k = remove(k, 1);
                elseif z == length(x)-1
                    x = remove(x,z);
                    k = remove(k, z);
                else
                    x(z+1) = (x(z+1)+x(z))/2;
                    x = remove(x,z);
                    k = remove(k, z);
                end
            end
            z = z+1;
        end
        
        if mod(round(t,2), 10*dt) == 0
            barrier = find(circshift(k,-1)-k == (K2-K1)); % K1 > K2
            if barrier > 0
                if (barrier > prof_size) && (barrier < (length(x)-prof_size-1))
                    dens = (x(barrier-prof_size+1:barrier+prof_size+1)-x(barrier-prof_size:barrier+prof_size)).^-1;
                    if sum(abs(dens)>15) == 0
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
plot(linspace(-prof_size,prof_size,length(full_density)),full_density,'.', 'MarkerSize', 10)
hold on
xlabel('$$x$$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$$q$$', 'Interpreter','latex', 'FontSize', 20)
legend('$$K_{1}-K_{2}=0.5$$', '$$K_{1}-K_{2}=1$$', '$$K_{1}-K_{2}=2$$', '$$K_{1}-K_{2}=5$$', 'Interpreter','latex', 'FontSize', 15)
%legend('$$K_{1}-K_{2}=0.5$$', 'Interpreter','latex', 'FontSize', 15)
grid on
grid minor
ylim([0.8, 1.2])



