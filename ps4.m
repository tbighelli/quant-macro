%% mia prima
clear all
clc
tic
beta   = 0.988;              % subjective discount factor 
delta  = 0.987;               % 1 - depreciation rate
theta=0.679;
%discretization
kss=(((1/(1-theta))*((1/beta)+0.013-1)))^(-1/theta);
for i =1:100
    k(i)= 0.05+(i-1)*((2*kss-0.05)/(100-1));
end

k=k';

v=zeros(100,1);



m=zeros(100,100);

for i=1:100
    for j=1:100
        m(i,j)=log(1-k(j,1)+delta*k(i,1));
       if k(j,1)>(1+delta*k(i,1));
           m(i,j)=-100;
       end
    end
end
dif=10000;
iter=0;
while dif>0.000001      
for i=1:100
    for j=1:100
            chi(i,j)=m(i,j)+beta*v(j);
            [v1(i,1),pos_idx(i,1)]=(max(chi(i,:)));
            [row_idx,col_idx]=ind2sub(size(chi),pos_idx);

             
    end
 end
dif=norm(v1-v);
v=v1;
iter=iter+1;
end
gk=k(row_idx);
figure
hold on
plot(v1) %ex 1.a
xlabel('k')
ylabel('v(k)')
hold off
figure
hold on
plot(k,gk);
xlabel('k');
ylabel('g(k)');

hold off
toc



%% mia prima 1.b
clear all
clc
tic
beta   = 0.988;              % subjective discount factor 
delta  = 0.987;               % 1 - depreciation rate
theta=0.679;
%discretization
kss=(((1/(1-theta))*((1/beta)+0.013-1)))^(-1/theta);
for i =1:100
    k(i)= 0.05+(i-1)*((2*kss-0.05)/(100-1));
end

k=k';

v=zeros(100,1);



m=zeros(100,100);

for i=1:100
    for j=1:100
        m(i,j)=log(1-k(j,1)+delta*k(i,1));
       if k(j,1)>(1+delta*k(i,1));
           m(i,j)=-100;
       end
    end
end
dif=10000;
iter=0;
gk(1)=0.05;
while dif>0.000001      
for i=1:100
    for j=1:100
        if  log(1-k(j,1)+delta*k(i,1))>gk(1);
            chi(i,j)=m(i,j)+beta*v(j);
            [v1(i,1),pos_idx(i,1)]=(max(chi(i,:)));
            [row_idx,col_idx]=ind2sub(size(chi),pos_idx);
            gk=k(row_idx);
        end
             
    end
 end
dif=norm(v1-v);
v=v1;
iter=iter+1;
end

figure
hold on
plot(v1) %ex 1.a
xlabel('k')
ylabel('v(k)')
hold off
figure
hold on
plot(k,gk);

hold off
toc
%% 1.c
n_k = 100;                

kss = [(1-theta)/(beta^(-1)-(1-delta))]^(1/theta);
kgrid = linspace(0.05*kss,2*kss,n_k);                      
k = kgrid';  

v0 = zeros(n_k,1);    
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
kmat = repmat(k,1,n_k);               
kpmat = repmat(kgrid,n_k,1);

c = kmat.^(1-theta)+(1-delta)*kmat-kpmat;               
neg = find(c<0);                                       
c(neg) = 0;   
utility = log(c);

v1 = zeros(n_k,1); 
k1 = zeros(n_k,1); 

tic
 
while epsilon > 0.00001  
    
    M = utility+beta*(v0mat);
    
    for i =1:n_k   
        for j = 1:n_k
        if M(i,j) > M(i,j+1)
            v1(i) = M(i,j);
            k1(i) = k(j);
            break
        end
        end
    end
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1  
    
end

toc

set(0,'defaultAxesFontName', 'Times New Roman')
subplot(2,2,1)
plot(k,k1);
title('Policy Function for Capital')
xlabel('$k$','Interpreter','latex'); ylabel('$g(k)$','Interpreter','latex');

subplot(2,2,2)
plot(k,v0);
title('Value Function')
xlabel('$k$','Interpreter','latex'); ylabel('$v(k)$','Interpreter','latex');
%% 1.d
clear all;
clc;

theta = .679;
beta = .988; 
delta = .013;
h = 1; 

n_k = 100;                

kss = [(1-theta)/(beta^(-1)-(1-delta))]^(1/theta);
kgrid = linspace(0.9*kss,1.1*kss,n_k);                      
k = kgrid';  

v0 = zeros(n_k,1);         
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
kmat = repmat(k,1,n_k);               
kpmat = repmat(kgrid,n_k,1);

c = kmat.^(1-theta)+(1-delta)*kmat-kpmat;               
neg = find(c<0);                                       
c(neg) = 0;   
utility = log(c);

tic
 
while epsilon > 0.00001 
    
    M = utility+beta*(v0mat);
    
    if its < 1400
    
    [v1,indexes]=max(M,[],2);
    k1 = k(indexes);
    
    % local search
    else 
        
        for i = 1:n_k
        [vp] = max(M(i,indexes(i,1)-1:indexes(i,1)+1),[],2);
        ind = find(M(i,:) == vp);
        indexes(i) = ind;
        v1(i) = vp;
        k1(i) = k(ind);
        end    
    end
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1 
end

toc

set(0,'defaultAxesFontName', 'Times New Roman')
subplot(2,2,1)
plot(k,k1);
title('Policy Function for Capital')
xlabel('$k$','Interpreter','latex'); ylabel('$g(k)$','Interpreter','latex');

subplot(2,2,2)
plot(k,v0);
title('Value Function')
xlabel('$k$','Interpreter','latex'); ylabel('$v(k)$','Interpreter','latex');

%% Ex 1.e

n_k = 100;                

kss = [(1-theta)/(beta^(-1)-(1-delta))]^(1/theta);
kgrid = linspace(0.05*kss,2*kss,n_k);                      
k = kgrid';  

v0 = zeros(n_k,1);    
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
kmat = repmat(k,1,n_k);               
kpmat = repmat(kgrid,n_k,1);

c = kmat.^(1-theta)+(1-delta)*kmat-kpmat;               
neg = find(c<0);                                       
c(neg) = 0;   
utility = log(c);

v1 = zeros(n_k,1); 
k1 = zeros(n_k,1); 

utility = triu(utility);
v0mat = triu(v0mat);

tic

 
while epsilon > 0.00001 
    
    M = utility+beta*(v0mat);
    
    for i =1:n_k   
        for j = 1:n_k - 1
        if M(i,j) > M(i,j+1)
            v1(i) = M(i,j);
            k1(i) = k(j);
            break
        end
        end
    end
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1  
    
end

toc

set(0,'defaultAxesFontName', 'Times New Roman')
subplot(2,2,1)
plot(k,k1);
title('Policy Function for Capital')
xlabel('$k$','Interpreter','latex'); ylabel('$g(k)$','Interpreter','latex');

subplot(2,2,2)
plot(k,v0);
title('Value Function')
xlabel('$k$','Interpreter','latex'); ylabel('$v(k)$','Interpreter','latex');
theta = .679;
beta = .988; 
delta = .013;
h = 1; 
%% 1.f,g
n_k = 100;                

kss = [(1-theta)/(beta^(-1)-(1-delta))]^(1/theta);
kgrid = linspace(0.05*kss,2*kss,n_k);                      
k = kgrid';  

v0 = zeros(n_k,1);         
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
kmat = repmat(k,1,n_k);               
kpmat = repmat(kgrid,n_k,1);

c = kmat.^(1-theta)+(1-delta)*kmat-kpmat;               
neg = find(c<0);                                       
c(neg) = 0;   
utility = log(c); 

tic
 
while epsilon > 0.00001
    
    [v1,indexes]=max(utility+beta*(v0mat),[],2);
    k1 = k(indexes);
    
    if  its == 400
    c1 = k.^(1-theta)+(1-delta)*k-k1;
    neg = find(c1<0);                                       
    c1(neg) = 0;
    u = log(c1);
    v{1} = v1;
    
    for j =1:9
        
    v{j+1} = u + beta*v{j};
    v1 = [v{j+1}];

    end
    end
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1  
    
end

toc

set(0,'defaultAxesFontName', 'Times New Roman')
subplot(2,2,1)
plot(k,k1);
title('Policy Function for Capital')
xlabel('$k$','Interpreter','latex'); ylabel('$g(k)$','Interpreter','latex');

subplot(2,2,2)
plot(k,v0);
title('Value Function')
xlabel('$k$','Interpreter','latex'); ylabel('$v(k)$','Interpreter','latex');

%% 1.2 labor
clear all; 
close all;
clc;

theta = .679;
beta = .988; 
delta = .013;
kappa = 5.24; 
nu = 2;

n_k = 100;                

% grid 
kgrid = linspace(20,48,n_k);                      
k = kgrid';  
hgrid = linspace(0, 1,n_k);
h = hgrid';

% Initialize matricies
v0 = zeros(n_k,1);         
v0mat = repmat(v0,1,n_k)';
v0mat = v0mat(:)';
v0mat = repmat(v0mat,n_k,1);
iter = 0;
epsilon = 1;
kmat = repmat(k,1,n_k);               
hmat = repmat(hgrid,n_k,1);

% compute utility
c = cell(1,n_k);
u = cell(1,n_k);
for i = 1:n_k
kpmat = repmat(k(i,1),n_k,n_k);
c{1,i} = kmat.^(1-theta).*hmat.^theta+(1-delta)*kmat-kpmat;   
neg = find(c{1,i}<0);                                       
c{1,i}(neg) = 0; 
u{1,i} = log(c{1,i})-kappa*(hmat.^(1+1/nu))./(1+(1/nu));
end 

utility = cell2mat(u);

% VFI
while epsilon > 0.000001
    
    [v1,indexes]=max(utility+beta*(v0mat),[],2);
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0,1,n_k)';
    v0mat = v0mat(:)';
    v0mat = repmat(v0mat,n_k,1);
    iter = iter + 1  
    
end

% find the decision rules on the grid
k1 = zeros(n_k,1);
hp = zeros(n_k,1);

for i = 1:n_k
    for j =2:n_k
        if indexes(i) <= n_k
        k1(i) = k(1);
        hp(i) = h(indexes(i));
        elseif indexes(i)>n_k*(j-1) & indexes(i) <= j*n_k  
        ind = indexes(i) - n_k*(j-1);
        k1(i) = k(j);
        hp(i) = h(ind);
        end   
    end
 end    


figure
plot(k,k1);

xlabel('k'); 
ylabel('g(k)');

figure
plot(k,hp);

xlabel('k'); 
ylabel('g(k)');

figure
plot(k,v0);

xlabel('k'); 
ylabel('v(k)');

%%
clear all;
close all;
clc;

theta = .679;
beta = .988; 
delta = .013;
h = 1; 

% select the capital nodes
n_k = 30;                

kss = [(1-theta)/(beta^(-1)-(1-delta))]^(1/theta);
kgrid = linspace(0.9*kss,1.1*kss,n_k);                      
k = kgrid';  

% guess a value function 
n=5;
Theta = ones(n,1);

syms phi(kp) v0(kp);
phi(kp) = [chebyshevT([0,1,2,3,4],kp)];
v0(kp) = phi(kp)*Theta;

k1 = zeros(n_k,1);
it = 0;
count = 0;
epsilon = 1;

tic

while epsilon > 1e-6

% find the policy given the guess 

for i = 1:n_k
 
 %Tv =@(kp) - log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-kp) - beta*([chebyshevT([0,1,2,3,4],kp)]*Theta);
 %k1(i) = fmincon(Tv(i),41,[],[],[],[],k(1,1),k(n_k,1));
 
 Tv =@(kp) log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-kp) + beta*([chebyshevT([0,1,2,3,4],kp)]*Theta);
 syms foc(kp);
 foc(kp) = diff(Tv,kp);
 x = fsolve(matlabFunction(foc),40);
 k1(i,1) = x;
 count = count +1

end

% update the value and the Chebyshev coefficients by LS n_k > n
v1 = zeros(n_k,1);
for i = 1: n_k

v1(i) = log(k(i,1)^(1-theta)+(1-delta)*k(i,1)-k1(i,1)) + beta*v0(k1(i,1));

end

Phi = zeros(n_k,n);

for j = 1:n_k
    
    Phi(j,:) = phi(k(j,1));
    
end

Theta1 = (Phi'*Phi)^(-1)*Phi'*v1;

% check convergence otherwise iterate
epsilon = norm(Theta1-Theta);
Theta = Theta1;
it = it + 1

end

toc


