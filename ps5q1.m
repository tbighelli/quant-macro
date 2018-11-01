
clear all;
clc;
close all;
n_z = 2;
sigma = 2;
cbar = 100; 
r = 0.04;
rho = 0.06;
w = 1; 
z = [1-sigma_y,1+sigma_y];
Pi = [(1+gamma)/2 (1-gamma)/2; (1-gamma)/2 (1+gamma)/2]; 
beta = 1/(1+rho);
    syms qu(x,y) utility(x,y); 
    qu(x,y) = -.5*(x-y)^2;
    utility(x,y) = qu(x,y);
    y = cbar;
    
    a = zeros(n_k,1);
    a_min = -(1-sigma_y)/r+2; 
    a_max = 16; 
    for i = 1:n_k
    a(i) = a_min+(i-1)*(a_max-a_min)/(n_k-1);
    end

   syms crra(x,y) utility(x,y); 
   crra(x,y) = (x^(1-y)-1)/(1-y); 
   utility(x,y) = crra(x,y);    
   y = sigma;
      a = zeros(n_k,1);
   a_min = -(1-sigma_y)/r+2; 
   a_max = 16; 
   for i = 1:n_k
   a(i) = a_min+(i-1)*(a_max-a_min)/(n_k-1);
   end
end 

v0 = zeros(n_k,1);         
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);

c = (1+r)*amat + ones(n_k,n_k)*w - apmat;               
neg = find(c<0);    
u = double(utility(c,y));
u(neg) = -100000000;

while epsilon > 1e-6
    
    [v1,indexes]=max(u+beta*v0mat,[],2);
    a1 = a(indexes);  
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1  
    
end    


v = v0;
ap = a1;
c = (1+r)*a + ones(n_k,1)*w - ap;
    
else
    
% Plots
figure
plot(a,ap);
xlabel('a'); 
ylabel('a(a)');

figure
plot(a,c);
xlabel('a'); 
ylabel('c(a)');
figure
plot(a,v);
xlabel('a');
ylabel('v(a)');

%%
z = [1-sigma_y,1+sigma_y];
Pi = [(1+gamma)/2 (1-gamma)/2; (1-gamma)/2 (1+gamma)/2]; 
beta = 1/(1+rho);


    syms qu(x,y) utility(x,y); 
    qu(x,y) = -.5*(x-y)^2;
    utility(x,y) = qu(x,y);
    y = cbar;

   syms crra(x,y) utility(x,y); 
   crra(x,y) = x^(1-y)/(1-y); 
   utility(x,y) = crra(x,y);    
   y = sigma;


% grid 
a = zeros(n_k,1);
a_min = 0;
a_max = 16; 
for i = 1:n_k
    a(i) = a_min+(i-1)*(a_max-a_min)/(n_k-1);
end

amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);
c = (1+r)*amat + ones(n_k,n_k)*w - apmat;               
neg = find(c<0);                                        
u = double(utility(c,y));
u(neg) = -1000000;

v = cell(1,T); 
ap = cell(1,T);
c = cell(1,T);

[v1,indexes]=max(u,[],2);
ap{T} = a(indexes);
c{T} = (1+r)*a+ 1*w - ap{T};
v{T} = v1;


for t = T-1:-1:1
    vmat = repmat(v{t+1}',n_k,1);
    [v{t},indexes]=max(u + beta*vmat,[],2);
    ap{t} = a(indexes);
    c{t} = (1+r)*a+ 1*w - ap{t};
end

amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);
c1 = z(1)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c1<0);                                          
u1 = double(utility(c1,y));
u1(neg) = -100000000;
c2 =  z(2)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c2<0);                                          
u2 = double(utility(c2,y));
u2(neg) = -100000000;
u = [u1; u2];

v = cell(1,T); 
ap = cell(1,T);
c = cell(1,T);

[v1,indexes]=max(u,[],2);
ap{T} = a(indexes);
ap_1 = ap{T}(1:n_k,1);
ap_2 = ap{T}(n_k+1:end,1);
ap{T} = [ap_1 ap_2];
c_1 = z(1)*ones(n_k,1)*w + (1+r)*a-ap_1;
c_2 = z(2)*ones(n_k,1)*w + (1+r)*a-ap_2;
c{T} = [c_1 c_2];
v{T} = v1;


for t = T-1:-1:1
    v0_one = v{t+1}(1:n_k,1)'; 
    v0_two = v{t+1}(n_k+1:end,1)';
    v0 = [v0_one; v0_two];
    W = Pi*v0;
    W1 = repmat(W(1,:),n_k,1);
    W2 = repmat(W(2,:),n_k,1);
    [v{t},indexes]=max(u+beta*([W1; W2]),[],2);
    ap{t} = a(indexes);
    ap_1 = ap{t}(1:n_k,1);
    ap_2 = ap{t}(n_k+1:end,1);
    ap{t} = [ap_1 ap_2];
    c_1 = z(1)*ones(n_k,1)*w + (1+r)*a-ap_1;
    c_2 = z(2)*ones(n_k,1)*w + (1+r)*a-ap_2;
    c{t} = [c_1 c_2];
end



c5 = c{5};
c40 = c{40};
ap5 = ap{5};
ap40 = ap{40};
v5 = v{5};
v40 = v{40};
figure
plot(a,ap);
xlabel('a'); 
ylabel('a(a)');

figure
plot(a,c);
xlabel('a'); 
ylabel('c(a)');
figure
plot(a,v);
xlabel('a');
ylabel('v(a)');

%%
v0_one = zeros(n_k,1)'; 
v0_two = zeros(n_k,1)';
v0 = [v0_one; v0_two];
its = 0;
epsilon = 1;
amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);


c1 = z(1)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c1<0);                                          
u1 = double(utility(c1,y));
u1(neg) = -100000000;
c2 =  z(2)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c2<0);                                          
u2 = double(utility(c2,y));
u2(neg) = -100000000;
u = [u1; u2];

while epsilon > 1e-6
    
    W = Pi*v0;
    W1 = repmat(W(1,:),n_k,1);
    W2 = repmat(W(2,:),n_k,1);
    [v1,indexes]=max(u+beta*([W1; W2]),[],2);
    a1 = a(indexes);
    v0 = [v0(1,:)'; v0(2,:)'];
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0_one = v1(1:n_k,1)'; 
    v0_two = v1(n_k+1:end,1)';
    v0 = [v0_one; v0_two];
    its = its + 1  
    
end


 v_1 = v1(1:n_k,1);
 v_2 = v1(n_k+1:end,1);
 v = [v_1 v_2];
 ap_1 = a1(1:n_k,1);
 ap_2 = a1(n_k+1:end,1);
 ap = [ap_1 ap_2];
 c_1 = z(1)*ones(n_k,1)*w +(1+r)*a-ap_1;
 c_2 = z(2)*ones(n_k,1)*w +(1+r)*a-ap_2;
 c = [c_1 c_2];
 

end

% Discrete VFI with uncertainty ----------------------------------------------------------------------------------------------------
v0_one = zeros(n_k,1)'; 
v0_two = zeros(n_k,1)';
v0 = [v0_one; v0_two];
its = 0;
epsilon = 1;
amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);


c1 = z(1)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c1<0);                                          
u1 = double(utility(c1,y));
u1(neg) = -100000000;
c2 =  z(2)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c2<0);                                          
u2 = double(utility(c2,y));
u2(neg) = -100000000;
u = [u1; u2];

while epsilon > 1e-6
    
    W = Pi*v0;
    W1 = repmat(W(1,:),n_k,1);
    W2 = repmat(W(2,:),n_k,1);
    [v1,indexes]=max(u+beta*([W1; W2]),[],2);
    a1 = a(indexes);
    v0 = [v0(1,:)'; v0(2,:)'];
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0_one = v1(1:n_k,1)'; 
    v0_two = v1(n_k+1:end,1)';
    v0 = [v0_one; v0_two];
    its = its + 1  
    
end

% output UMP with sigma_y > 0
 v_1 = v1(1:n_k,1);
 v_2 = v1(n_k+1:end,1);
 v = [v_1 v_2];
 ap_1 = a1(1:n_k,1);
 ap_2 = a1(n_k+1:end,1);
 ap = [ap_1 ap_2];
 c_1 = z(1)*ones(n_k,1)*w +(1+r)*a-ap_1;
 c_2 = z(2)*ones(n_k,1)*w +(1+r)*a-ap_2;
 c = [c_1 c_2];
 
end
end



%crra
v0 = zeros(n_k,1);         
v0mat = repmat(v0',n_k,1);
its = 0;
epsilon = 1;
amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);

c = (1+r)*amat + ones(n_k,n_k)*w - apmat;               
neg = find(c<0);    
u = double(utility(c,y));
u(neg) = -100000000;

 
while epsilon > 1e-6
    
    [v1,indexes]=max(u+beta*v0mat,[],2);
    a1 = a(indexes);  
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0mat = repmat(v0',n_k,1);
    its = its + 1  
    
end    


v = v0;
ap = a1;
c = (1+r)*a + ones(n_k,1)*w - ap;
    
else
    

v0_one = zeros(n_k,1)'; 
v0_two = zeros(n_k,1)';
v0 = [v0_one; v0_two];
its = 0;
epsilon = 1;
amat = repmat(a,1,n_k);               
apmat = repmat(a',n_k,1);

% compute utility
c1 = z(1)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c1<0);                                          
u1 = double(utility(c1,y));
u1(neg) = -100000000;
c2 =  z(2)*ones(n_k,n_k)*w +(1+r)*amat-apmat;               
neg = find(c2<0);                                          
u2 = double(utility(c2,y));
u2(neg) = -100000000;
u = [u1; u2];

while epsilon > 1e-6
    
    W = Pi*v0;
    W1 = repmat(W(1,:),n_k,1);
    W2 = repmat(W(2,:),n_k,1);
    [v1,indexes]=max(u+beta*([W1; W2]),[],2);
    a1 = a(indexes);
    v0 = [v0(1,:)'; v0(2,:)'];
    
    epsilon = norm(v1-v0);
    
    v0 = v1; 
    v0_one = v1(1:n_k,1)'; 
    v0_two = v1(n_k+1:end,1)';
    v0 = [v0_one; v0_two];
    its = its + 1  
    
end


 v_1 = v1(1:n_k,1);
 v_2 = v1(n_k+1:end,1);
 v = [v_1 v_2];
 ap_1 = a1(1:n_k,1);
 ap_2 = a1(n_k+1:end,1);
 ap = [ap_1 ap_2];
 c_1 = z(1)*ones(n_k,1)*w +(1+r)*a-ap_1;
 c_2 = z(2)*ones(n_k,1)*w +(1+r)*a-ap_2;
 c = [c_1 c_2];
 
end

figure
plot(a,ap);
xlabel('a'); 
ylabel('a(a)');

figure
plot(a,c);
xlabel('a'); 
ylabel('c(a)');
figure
plot(a,v);
xlabel('a');
ylabel('v(a)');







