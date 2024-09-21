clc
clear
close all

addpath(genpath('..\..\..\Research\ResearchTools'));
addpath('P:\gurobi903\win64\matlab');
addpath(genpath('P:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64'));
gurobi_setup;

T = 30;

m = 4; n = 4;

% Tot states
tot_states = 3*n*m + 2;

X0 = 0*ones(tot_states,1);
X0(1) = 1;

[t,y] = ode45(@my_ode,[0 T],X0);

% Plots
figure;
hold on; box on;

x4 = y(:,end);
x5 = y(:,1);

subplot(2,1,1);
plth = plot(t,x4,'k');
title('Robots fraction ending at home after seeding');
xlabel('Time (s)');
ylabel('Population Fractions');
subplot(2,1,2);
plts = plot(t,x5,'k');
title('Robots fraction that are ready to leave');
xlabel('Time (s)');
ylabel('Population Fractions');

figure;
hold on; box on;

for i = 1:m*n
    x = y(:,i+1) + y(:,i+m*n+1) + y(:,i+2*m*n+1);
    subplot(m,n,i)
    plot(t,x);
    xlabel('Time (s)');
    ylabel('Population Fractions');
    
end

pause;

figure; hold on; box on; % ('WindowState','maximized')


mult = 0;
sub_ind = 1;

tim = [0, 1, 5, 30];

for tind = 1:length(t)
    
    if (t(tind) == 0) || ((t(tind) < 1.03) && (t(tind) > 1)) || ...
            ((t(tind) < 5.05) && (t(tind) > 5)) || (tind == length(t))
        mult = mult + 1;
        x_ax = 0:4;
        y_ax = 0:6;
        
        heat_lvls = 1e-9*ones(5,7);
        
        % Home heat
        heat_lvls(2:3,1) = y(tind,1) + y(tind,end);
        
        for r = 0:m-1
            row_patch = y(tind,(r*n+1:r*n+n)+1) + ...
                y(tind,m*n+1+(r*n+1:r*n+n)) + ...
                y(tind,(r*n+1:r*n+n)+2*m*n+1);
            heat_lvls(1:end-1,r+3) = row_patch;
        end
        
        subplot(2,2,sub_ind);
        pcolor(x_ax, y_ax, heat_lvls');
        colorbar; set(gca,'colorscale','log');
        
        xlabel('P_x'); ylabel('P_y');
        title(['Robot Distribution at t = ' num2str(tim(sub_ind))]);
        
        sub_ind = sub_ind + 1;
        
    end
end

close all;

figure('WindowState','maximized'); hold on; box on; % 

xlabel('P_x'); ylabel('P_y');

title('Robot fractions on the field');

mult = 0;

for tind = 1:length(t)
    
    if (t(tind) >= 0.1*mult)
        mult = mult + 1;
        x_ax = 0:4;
        y_ax = 0:6;
        
        heat_lvls = 1e-9*ones(5,7);
        
        % Home heat
        heat_lvls(2:3,1) = y(tind,1) + y(tind,end);
        
        for r = 0:m-1
            row_patch = y(tind,(r*n+1:r*n+n)+1) + ...
                y(tind,m*n+1+(r*n+1:r*n+n)) + ...
                y(tind,(r*n+1:r*n+n)+2*m*n+1);

%             row_patch = 3*row_patch;
            heat_lvls(1:end-1,r+3) = row_patch;
        end
        
        
        
        pcolor(x_ax, y_ax, heat_lvls');
        colorbar;
        
        pause(0.1);
    end
    
    if (t(tind) > 20)
        break;
    elseif (tind == 1)
        set(gca,'colorscale','log');
        pause(2);
    end
end

% Properties of model:
check_props(X0);

            
%%

function [x_dot] = my_ode(t,x)

m = 4; n = 4;

a = ones(1,m); a(2) = 1.4;
b = 0.2*ones(1,n);

k = 2; l = 2;
K_const = 1; d = 0.1;

K = zeros(3*n*m+2);

K(1,1) = -n*a(1);
K(1,end) = d;
K(end,1+3*m*n-m+1:end) = [K_const*ones(1,n), -d];

ind = 2;

% For the U, S, D states
for inda = 0:m*n-1
    Ka_temp = zeros(1,3*n*m+2);
    Ka_temp(1+inda+1) = -k;
    Z = zeros(1,n*m);
    if (inda < n)
        Ka_temp(1) = a(1);
    else
        Z(inda+1-n) = a(fix(inda/n)+1);
    end
    if (mod(inda,n)+1) >= 2
        Z(inda) = b(fix(inda/n)+1);
    end
    if (mod(inda,n)+1) <= n-1
        Z(inda+2) = b(fix(inda/n)+1);
    end
    Ka_temp(1+2*m*n+1:end-1) = Z;
    
    K(ind,:) = Ka_temp;
    
    Kb_temp = zeros(1,3*n*m+2);
    Kb_temp(1+inda+1) = k;
    Kb_temp(1+m*n+inda+1) = -l;
    
    K(ind+n*m,:) = Kb_temp;
    
    Kc_temp = zeros(1,3*n*m+2);
    Kc_temp(1+m*n+inda+1) = l;
    Z = zeros(1,n*m);
    if (fix(inda/n)+1) <= (m-1)
        Z(inda+1) = -a(fix(inda/n)+2);
    else
        Z(inda+1) = -K_const;
    end
    if (mod(inda,n)+1) >= 2
        Z(inda+1) = Z(inda+1) - b(fix(inda/n)+1);
    end
    if (mod(inda,n)+1) <= n-1
        Z(inda+1) = Z(inda+1) - b(fix(inda/n)+1);
    end
    Kc_temp(1+2*m*n+1:end-1) = Z;
    
    K(ind+2*n*m,:) = Kc_temp;
    ind = ind + 1;
    
end
   
M = eye(3*n*m+2);
 
y_x = x;
   
x_dot = M*K*y_x;
end

%% Functions used to check properties:

function [outtie] = check_props(x_samp)

m = 4; n = 4;

a = ones(1,m); a(2) = 1.4;
b = 0.2*ones(1,n);

k = 2; l = 2;
K_const = 1; d = 0.1;

K = zeros(3*n*m+2);

K(1,1) = -n*a(1);
K(1,end) = d;
K(end,1+3*m*n-m+1:end) = [K_const*ones(1,n), -d];

ind = 2;

    % For the U, S, D states
    for inda = 0:m*n-1
        Ka_temp = zeros(1,3*n*m+2);
        Ka_temp(1+inda+1) = -k;
        Z = zeros(1,n*m);
        if (inda < n)
            Ka_temp(1) = a(1);
        else
            Z(inda+1-n) = a(fix(inda/n)+1);
        end
        if (mod(inda,n)+1) >= 2
            Z(inda) = b(fix(inda/n)+1);
        end
        if (mod(inda,n)+1) <= n-1
            Z(inda+2) = b(fix(inda/n)+1);
        end
        Ka_temp(1+2*m*n+1:end-1) = Z;

        K(ind,:) = Ka_temp;

        Kb_temp = zeros(1,3*n*m+2);
        Kb_temp(1+inda+1) = k;
        Kb_temp(1+m*n+inda+1) = -l;

        K(ind+n*m,:) = Kb_temp;

        Kc_temp = zeros(1,3*n*m+2);
        Kc_temp(1+m*n+inda+1) = l;
        Z = zeros(1,n*m);
        if (fix(inda/n)+1) <= (m-1)
            Z(inda+1) = -a(fix(inda/n)+2);
        else
            Z(inda+1) = -K_const;
        end
        if (mod(inda,n)+1) >= 2
            Z(inda+1) = Z(inda+1) - b(fix(inda/n)+1);
        end
        if (mod(inda,n)+1) <= n-1
            Z(inda+1) = Z(inda+1) - b(fix(inda/n)+1);
        end
        Kc_temp(1+2*m*n+1:end-1) = Z;

        K(ind+2*n*m,:) = Kc_temp;
        ind = ind + 1;

    end
    
    fprintf('Eigenvalues of K:\n');
    lambda_K = eig(K)
    
    % 1: Find model equilibrium state
    x_e = sdpvar(size(x_samp,1),1, 'full');
    constr = [(K*x_e == 0); (sum(x_e) == 1)];
    ops = sdpsettings('verbose',0,'debug',0, 'solver','gurobi');
    optim1 = optimize(constr, [], ops);
    
    if optim1.problem ~= 0
        warning(['Equilibrium not found.']);
    else
        fprintf('\nEquilibrium value:\n');
        x_e = value(x_e)
    end
    
    % 2: Find if model is expontentially stable
    P = sdpvar(3*n*m+2,3*n*m+2);
    constr = [P*K + (K')*P <= -1*eye(3*n*m+2), P >= 1*eye(3*n*m+2)]; % neg def
    ops = sdpsettings('verbose',5,'debug',5, 'solver','csdp');
    optim1 = optimize(constr, [], ops);
    
    if optim1.problem ~= 0
        constr = [P*K + (K')*P <= zeros(3*n*m+2), P >= 1*eye(3*n*m+2)]; % neg def
        optim1 = optimize(constr, [], ops);
        if optim1.problem ~= 0
            warning('Netwotk not stable');
        else
            fprintf('\nNegative semidefinite derivative of Lyapunov found, System stable.\n');
        end
    else
        fprintf('\nNegative definite derivative of Lyapunov found, System asymptotically stable.\n');
        x_e = value(x_e)
    end

    outtie = 1;
end