N0 = 100; % initial cell numbers
A = 0; B = 100; L = B-A; % boundaries
l0 = L/N0; a0 = L/N0;
tmax = 5000;
dt = 0.01;

cols = [[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]];
d_vals = [0.5*a0, 0.75*a0, 1.5*a0, 2*a0]; % l_0
xp = dt/50; % death probability
K1 = 1; K2_vals=[1.1,1.2,1.3,1.4,1.5]; 

tic
% linear function for cell division prob
Lin = @(xp, d, L) xp.*L./d;

% insert/remove new cell upon div
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); 
remove = @(x, n)cat(2,  x(1:n-1), x(n+1:end));

times = zeros(1,length(K2_vals)*length(d_vals));
p=1;
for d = d_vals

for K2=K2_vals
    sim_num = 0; max_sim=50;
    mean_bdry = zeros(1,round(tmax/dt));
    while sim_num<max_sim
    
        k = [zeros(1, N0/2) + K1, zeros(1, N0/2) + K2]; % spring const.
        a = zeros(1, N0) + a0; % equilibrium length
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
                x(i) = x(i)+dt*(k(i)*(x(i+1)-x(i)-a(i))-k(i-1)*(x(i)-x(i-1)-a(i-1)));
                if x(i) > 100 
                    x(i) = 100;
                elseif x(i) < 0
                    x(i) = 0;
                end
                i = i+1;
            end

            % cell division (Linear function)
            rand_vals = rand(1, length(x)-1);
            spaces = circshift(x,-1)-x; spaces = spaces(1,1:end-1);
            lins = Lin(xp,d,spaces);
            buffer = 0;
            
            for j=1:length(lins)
                 if rand_vals(j)<=lins(j)
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
                if v == length(k) && k(v) == K1
                    mean_bdry(c) = mean_bdry(c) + x(v+1); % K1 invades K2
                    v = 0;
                    flag = 0;
                elseif v == length(k) && k(v) == K2
                    mean_bdry(c) = mean_bdry(c) + 0; % K2 invades K1
                    v = 0;
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
    plot(linspace(0,tmax,length(mean_bdry)), boundary,'-')
    hold on
    minimum = min(abs(abs(boundary-50)-10));
    if boundary(end) < 50
        times(p) = -dt*find(abs(abs(boundary-50)-10) == minimum);
    elseif boundary(end) > 50
        times(p) = dt*find(abs(abs(boundary-50)-10) == minimum);
    end        
    p=p+1;
end

end

ylabel('Cell boundry position')
xlabel('time t')
grid on
grid minor
toc

K2_plot = linspace(K2_vals(1),K2_vals(end),100);
theory = @(l0) -(l0-a0).*sqrt(xp/dt).*(K1-K2_plot).*sqrt(K1.*K2_plot)./(K1.^(3/2)+K2_plot.^(3/2));

figure
u = length(K2_vals);
col_ind = 1;

plot1 = plot(K2_vals-K1 , 10./times(1:u), '-o');
hold on
plot2 = plot(K2_plot-K1, theory(d_vals(1)), '--');
set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
col_ind = col_ind + 3;

hold on
plot1 = plot(K2_vals-K1 , 10./times(u+1:u*2), '-o');
hold on
plot2 = plot(K2_plot-K1, theory(d_vals(2)), '--');
set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
col_ind = col_ind + 3;

hold on
plot1 = plot(K2_vals-K1 , 10./times(u*2+1:u*3), '-o');
hold on
plot2 = plot(K2_plot-K1, theory(d_vals(3)), '--');
set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
col_ind = col_ind + 3;

hold on
plot1 = plot(K2_vals-K1 , 10./times(u*3+1:end), '-o');
hold on
plot2 = plot(K2_plot-K1, theory(d_vals(4)), '--');
set([plot1 plot2],'Color', cols(col_ind:col_ind+2));
col_ind = col_ind + 3;


ylabel('$$v$$', 'Interpreter','latex', 'FontSize', 15)
xlabel('$$k_{2} - k_{1}$$', 'Interpreter','latex', 'FontSize', 15)
legend("$$l_{0}=0.5$$","","$$l_{0}=0.75$$","","$$l_{0}=1.5$$","","$$l_{0}=2$$","", 'Interpreter','latex', 'FontSize', 15)
title('$$Linear: d=2 \cdot 10^{-4}$$', 'Interpreter','latex', 'FontSize', 15)
grid on
grid minor