%% Code for solving a simplified version of K&S (1998)
% Quantitative Macroeconomics - IDEA programme



%% Initial values and parameters

%%%%%%%%%%%% Finding the transition matrix for the state %%%%%%%%%%%%%%

% The system of equations 
A= [ 1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 ; ...
     0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 ; ...
     0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ; ...
     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0 5.6 0 -1  0  0  0  0  0 ; ...
    -1 0 28/3 0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
  .02 .48 .05 .45 0 0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0 0 .02 .48 .05 .45 0 0 0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 ];
  
 
b= [7/8; 7/8; 7/8; 7/8; 1/8; 1/8; 1/8; 1/8; 7/24; 21/40; 0; 0; 0.02; 0.005; 0.05; 0.02];


pize = reshape(A^-1*b,4,4);
pie_g0 = pize(1:2,1:2);
pie_b0 = pize(3:4,3:4);   
for i=1:2
    for j=1:2
        pie_g(i,j)=pie_g0(i,j)/sum(pie_g0(i,:));
    end
end
for i=1:2
    for j=1:2
        pie_b(i,j)=pie_b0(i,j)/sum(pie_b0(i,:));
    end
end


% transtion matrix aggregate state

piZ = [ 7/8  1/8;...
        1/8  7/8];
    
    
%%%%%%%%%%%%  Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%

betta=0.95;
delta=0.0025;
z=[0,1];
alpha=0.36;
kappa = .975;
psi = 1.5;
sigma = 2.5;
gamma = 2.5;

%L=[0.96, 0.9];
Xi = (1/betta - 1 + delta)/alpha;
pi    = zeros(2,2);  
pi = pie_g;
r0  = (1/betta-1)-[10^(-12),10^(-4)];
w0 = (1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));
h_grid  = 20;   
k_grid  = 201; 
khratio = ((1/betta - 1 + delta)/alpha)^(1/(alpha-1));
chratio = khratio^alpha - delta*khratio;
kgrid =[0.1:0.1:5,5.3:0.3:50.3]';
hgrid = [1/20:1/20:1]';
for i=1:20
 for j = 1:201
     for m=1:2
         for n=1:201
           c(i,j,m,n)= (1+r0(m))*kgrid(n) +  w0(m)*hgrid(i)*z(m) + (1-delta)*kgrid(n) - kgrid(j);
           if c(i,j,m,n)<0
               c(i,j,m,n)=0;
           end
         u(i,j,m,n) = c(i,j,m,n)^(1-sigma)/(1-sigma) - psi*hgrid(i)^(1+1/nu)/(1+1/nu);
         end
     end
 end
end

v = zeros(201,2);
tol=0.00001; 
diff = 1; 
iter = 0;
while diff>tol
    diff = 0;
    for m= 1:2
        for n=1:201
            objfn(:,:,m,n) = u(:,:,m,n) + betta *(pi(m,1)*(v(:,1)*ones(1,20))' + pi(m,2)*(v(:,2)*ones(1,20))');
            Tv(n,m)=max(max(objfn(:,:,m,n)));
        end
    end
    diff = norm(v-Tv);
    v = Tv;
    iter = iter + 1;
end
for m=1:2
    for n = 1:201
        objfn(:,:,m,n)= u(:,:,m,n) + betta *(pi(m,1)*(v(:,1)*ones(1,20))' + pi(m,2)*(v(:,2)*ones(1,20))');
        [tmp1,x1] = max(objfn(:,:,m,n),[],1);
        [tmp2,x2] = max(tmp1,[],2);
        kgridrule(m,n)=x2;
        hgridrule(m,n)=x1(x2);
        
        kdecrule(m,n) = kgrid(kgridrule(m,n));
        hdecrule(m,n) = hgrid(hgridrule(m,n));
        
        cdecrule(m,n) = (1+r0(m))*kgrid(n) +  w0(m)*hdecrule(m,n)*z(m) + (1-delta)*kgrid(n) - kdecrule(m,n);
    end
end
figure
hold on;
plot(kgrid,kdecrule)
title('policy function of assets');
hold off
figure
hold on
plot(kgrid, hdecrule)
title('optimal labor decision');
hold off


%% Initial values and parameters

%%%%%%%%%%%% Finding the transition matrix for the state %%%%%%%%%%%%%%

% The system of equations 
A= [ 1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 ; ...
     0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0 ; ...
     0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ; ...
     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0  0 5.6 0 -1  0  0  0  0  0 ; ...
    -1 0 28/3 0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
  .02 .48 .05 .45 0 0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
     0  0  0  0  0  0  0 0 .02 .48 .05 .45 0 0 0  0 ; ...
     0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 ];
  
 
b= [7/8; 7/8; 7/8; 7/8; 1/8; 1/8; 1/8; 1/8; 7/24; 21/40; 0; 0; 0.02; 0.005; 0.05; 0.02];


pize = reshape(A^-1*b,4,4);
pie_g0 = pize(1:2,1:2);
pie_b0 = pize(3:4,3:4);   
for i=1:2
    for j=1:2
        pie_g(i,j)=pie_g0(i,j)/sum(pie_g0(i,:));
    end
end
for i=1:2
    for j=1:2
        pie_b(i,j)=pie_b0(i,j)/sum(pie_b0(i,:));
    end
end


% transtion matrix aggregate state

piZ = [ 7/8  1/8;...
        1/8  7/8];
    
    
%%%%%%%%%%%%  Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%

betta=0.95;
delta=0.0025;
z=[0,1];
alpha=0.36;
kappa = .975;
psi = 1.5;
sigma = 2.5;
gamma = 2.5;

%L=[0.96, 0.9];
Xi = (1/betta - 1 + delta)/alpha;
pi    = zeros(2,2);  
pi = pie_b;
r0  = (1/betta-1)-[10^(-12),10^(-4)];
w0 = (1-alpha)*(alpha./(r0+delta)).^(alpha/(1-alpha));
h_grid  = 20;   
k_grid  = 201; 
khratio = ((1/betta - 1 + delta)/alpha)^(1/(alpha-1));
chratio = khratio^alpha - delta*khratio;
kgrid =[0.1:0.1:5,5.3:0.3:50.3]';
hgrid = [1/20:1/20:1]';
for i=1:20
 for j = 1:201
     for m=1:2
         for n=1:201
           c(i,j,m,n)= (1+r0(m))*kgrid(n) +  w0(m)*hgrid(i)*z(m) + (1-delta)*kgrid(n) - kgrid(j);
           if c(i,j,m,n)<0
               c(i,j,m,n)=0;
           end
         u(i,j,m,n) = c(i,j,m,n)^(1-sigma)/(1-sigma) - psi*hgrid(i)^(1+1/nu)/(1+1/nu);
         end
     end
 end
end

v = zeros(201,2);
tol=0.00001; 
diff = 1; 
iter = 0;
while diff>tol
    diff = 0;
    for m= 1:2
        for n=1:201
            objfn(:,:,m,n) = u(:,:,m,n) + betta *(pi(m,1)*(v(:,1)*ones(1,20))' + pi(m,2)*(v(:,2)*ones(1,20))');
            Tv(n,m)=max(max(objfn(:,:,m,n)));
        end
    end
    diff = norm(v-Tv);
    v = Tv;
    iter = iter + 1;
end
for m=1:2
    for n = 1:201
        objfn(:,:,m,n)= u(:,:,m,n) + betta *(pi(m,1)*(v(:,1)*ones(1,20))' + pi(m,2)*(v(:,2)*ones(1,20))');
        [tmp1,x1] = max(objfn(:,:,m,n),[],1);
        [tmp2,x2] = max(tmp1,[],2);
        kgridrule(m,n)=x2;
        hgridrule(m,n)=x1(x2);
        
        kdecrule(m,n) = kgrid(kgridrule(m,n));
        hdecrule(m,n) = hgrid(hgridrule(m,n));
        
        cdecrule(m,n) = (1+r0(m))*kgrid(n) +  w0(m)*hdecrule(m,n)*z(m) + (1-delta)*kgrid(n) - kdecrule(m,n);
    end
end
figure
hold on;
plot(kgrid,kdecrule)
title('policy function of assets');
hold off
figure
hold on
plot(kgrid, hdecrule)
title('optimal labor decision');
hold off
%%
