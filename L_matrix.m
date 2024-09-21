%%

T = 30;

m = 4; n = 4;

% Tot states
tot_states = 3*n*m + 2;

X0 = 0*ones(tot_states,1);
X0(1) = 1;

%% K-matrix

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

%% L matrix

x = sym('x',[tot_states 1]);
m = 1;
F = sym('F',[tot_states 1]);

for i = 1:tot_states
    for j = 1:tot_states
        if (i==j); F(i,j) = subs(F(i,j),0);continue;end
        F(i,j) = -K(i,j)*(x(i)-x(j));
        s(i)   = sum(F(i,:));
    end
end

L = -equationsToMatrix(s(:));
L = double(L);
%% 

func = @(t,x,L) -L*x;

[t,y] = ode45(@(t,x) func(t,x,L),[0 T],X0);
%%
plot(t,y);
xlabel('Time (s)');
ylabel('Information states');
Ls= 1/2*(L+L');
Eigenvalues = eig(Ls) ;

%% R matrix for zero deficiency Theorem 
m = 4; n = 4;

M = eye(3*n*m+2); % M-matrix

Edge_List = {[1,2];[1,5];[1,8];[1,11];[2,3];[3,4];[4,5];...
             [5,6];[6,7];[7,8];[8,9];[9,10];[10,11];[11,12];...
             [12,13];[4,14];[14,15];[15,16];[16,17];[17,18];...
             [18,19];[19,20];[20,21];[21,22];[22,23];[23,24];...
             [24,25];[16,26];[7,17];[10,20];[13,23];[26,27];...
             [27,28];[28,29];[29,30];[34,35];[30,31];[31,32];...
             [32,33];[33,34];[19,29];[22,32];[25,35];[34;35];...
             [35,36];[36,37];[28,38];[38,39];[39,40];[40;41];...
             [41,42];[42,43];[43,44];[44,45];[45,46];[31;41];...
             [34,44];[46,47];[47,48];[49,50];[37,47];[40;50];...
             [43,50];[46,50]};
        

for i = 1:length(Edge_List)
    R{i,1} = M(:,Edge_List{i}(1))'-M(:,Edge_List{i}(2))';
end

R = cell2mat(R);

Rank = rank(R);