N0 = 50; % initial cell numbers
A = 0; B = 50; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
tmax = 500;
dt = 0.01;

d_vals = [0.5,0.75,1.25,1.5].*l0; % l_0
xp = dt/10; % death probability

K1 = 1; K2_vals=[1.5,2,2.5,3,3.5,4]; 

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

% linear function for cell division prob
Lin = @(xp, d, L) xp.*L./d;

tic
times = zeros(1,length(K2_vals)*length(d_vals));
dist = 5;

p=1;
for d = d_vals

for K2=K2_vals
    sim_num = 0; max_sim=50;
    mean_bdry = zeros(1,round(tmax/dt)+1);
    while sim_num<max_sim
    
        k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % spring const.
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
                        k = remove(k,z);
                    end
                end
                z = z+1;
            end
            t = t+dt;

            flag = 1; v = 1;
            while flag == 1
                if v == length(k) && k(v) == K1
                    mean_bdry(c) = mean_bdry(c) + x(v+1); % K1 invades K2
                    flag = 0;
                elseif v == length(k) && k(v) == K2
                    mean_bdry(c) = mean_bdry(c) + 0; % K2 invades K1
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
    minimum = min(abs(abs(boundary-B/2)-dist));
    if boundary(end) < B/2
        times(p) = -dt*find(abs(abs(boundary-B/2)-dist) == minimum);
    elseif boundary(end) > B/2
        times(p) = dt*find(abs(abs(boundary-B/2)-dist) == minimum);
    end        
    p=p+1;
    plot(linspace(0,tmax,length(boundary)),boundary)
    hold on
end

end

toc

figure
u = length(K2_vals);

plot(K2_vals-K1 , -dist./times(1:u), '-o', LineWidth=1.3);

hold on
plot(K2_vals-K1 , -dist./times(u+1:u*2), '-o', LineWidth=1.3);

hold on
plot(K2_vals-K1 , -dist./times(u*2+1:u*3), '-o', LineWidth=1.3);

hold on
plot(K2_vals-K1 , -dist./times(u*3+1:end), '-o', LineWidth=1.3);


ylabel('$$v$$', 'Interpreter','latex', 'FontSize', 20)
xlabel('$$K_{1} - K_{2}$$', 'Interpreter','latex', 'FontSize', 15)
legend("$$l^{*}=0.5$$", "$$l^{*}=0.75$$", "$$l^{*}=1.25$$", "$$l^{*}=1.5$$", 'Interpreter','latex', 'FontSize', 15)
grid on
grid minor