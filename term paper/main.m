
 
clc
clear all

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


piZ = [ 7/8  1/8;...
        1/8  7/8];
    
%%%%%%%%%%%%  Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.95;
delta=0.0025;
h=0.8;
ex=[1.2 0.8];
alfa=0.3;
cf=ones(1,201);
sigma=0.2;
L=[0.96, 0.9]; 
%%%%%%%%%%%%% Starting values for V %%%%%%%%%%%%%%%%%%%     

v1g = @(k,K) (((-ex(1)*cf+ ((alfa*(K/L(1))^(alfa-1))+ex(1)-1).*k+ h*(1-alfa-0.1)*(K/L(1))^(alfa) -delta.*k )).^((sigma-1)/sigma)+(cf*ex(1)).^((sigma-1)/sigma)).^(sigma/(sigma-1)/(1-betta));
v1b = @(k,K) (((-ex(2)*cf+ ((alfa*(K/L(2))^(alfa-1))+ex(2)-1).*k+ h*(1-alfa-0.1)*(K/L(2))^(alfa) -delta.*k )).^((sigma-1)/sigma)+(cf*ex(2)).^((sigma-1)/sigma)).^(sigma/(sigma-1)/(1-betta));
v0g = @(k,K) (((-ex(1)*cf+ ((alfa*(K/L(1))^(alfa-1))+ex(1)-1).*k -delta.*k )).^((sigma-1)/sigma)+(cf*ex(1)).^((sigma-1)/sigma)).^(sigma/(sigma-1)/(1-betta));
v0b = @(k,K) (((-ex(2)*cf+ ((alfa*(K/L(2))^(alfa-1))+ex(2)-1).*k-delta.*k )).^((sigma-1)/sigma)+(cf*ex(2)).^((sigma-1)/sigma)).^(sigma/(sigma-1)/(1-betta));

%%%%%%%%%%%%% Grid for k and K %%%%%%%%%%%%%%%%%%%%%%%%


k_grid=[0.1:0.1:5.1,5.4:0.3:50.1]; 
K_grid=[16:0.04:18.5];      




for j=1:size(K_grid,2)
V1g(:,j)= v1g(k_grid,K_grid(j))';
V1b(:,j)= v1b(k_grid,K_grid(j))';
V0g(:,j)= v0g(k_grid,K_grid(j))';
V0b(:,j)= v0b(k_grid,K_grid(j))';
end


% initial values
b0g=0;
b1g=1;
b2g=1;
b0b=0;
b1b=1;
b2b=1;


for iter_b=1:100
iter_b


H=@(K,ei) exp( (b0g+b1g*log(K)+b2g*(ei))+ (b0b+b1b*log(K)+b2b*((1-ei))));


Ha= @(K,ei) min(abs(K_grid-H(K,ei))); 


% e=1 employed



c= @(i,I,e,g) max(((alfa*(K_grid(I)/L(g))^(alfa-1))+ex(g)-1).*k_grid(i)+ ...
             h*(1-alfa)*(K_grid(I)/L(g))^(alfa)*e +(1-delta)*k_grid(i) ...
             - k_grid-ex(g),0) ;
       
for iter=1:1000 % VFI         
for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)        
        
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); % Ha
        
       
        
 
        V0gt(i,I)= max((((((c(i,I,0,1).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
         V1gt(i,I)= max((((((c(i,I,1,1).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
     
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0)));
         V0bt(i,I)= max((((((c(i,I,0,2).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
        V1bt(i,I)= max(((((((c(i,I,1,2)).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));

     
    end     
end

dev= max(max(abs( [V0gt-V0g,V1gt-V1g,V0bt-V0b,V1bt-V1b])));

if dev<0.0001
    break
else
    V0g=V0gt;
    V1g=V1gt;
    V0b=V0bt;
    V1b=V1bt;
end 
end          
        

for i=1:size(k_grid,2)
    for I=1:size(K_grid,2)
    
      
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),1))); 
        [V0gt(i,I),a(i,I,2,1)]= max((((((c(i,I,0,1).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(1,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
        [V1gt(i,I),a(i,I,1,1)]= max((((((c(i,I,1,1).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(2,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
        
       
        [dif,Ip]=min(abs(K_grid-H(K_grid(I),0))); % H(:,1): bad shock
        [V0bt(i,I),a(i,I,2,2)]= max((((((c(i,I,0,2).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(3,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
        [V1bt(i,I),a(i,I,1,2)]= max(((((((c(i,I,1,2)).^((sigma-1)/sigma))+cf.^((sigma-1)/sigma))).^(sigma/(sigma-1)))' + betta * ((([pize(4,:)]*([V0g(:,Ip),V1g(:,Ip),V0b(:,Ip),V1b(:,Ip)]'))'))));
        

        
    end 
end
%% Simulation


if iter_b==1  
ext(1)=1;  
for t=2:2000  
    draw=rand;
    ext(t)= 1+(rand>=piZ(ext(t-1),1));
end



extb=ext; 
extb(1:200)=0;


i_extg=find(ext==1);

i_extb=find(ext==2);

% initial distribution of assets and employment
indkk = find(k_grid==17); 
N_state(1:960,:,1)=ones(960,1)*[40,1]; 

% =2 unemployed
N_state(961:1000,:,1)=ones(40,1)*[40,2];

K_ind(1)=3;

for t=2:2000
for n=1:1000 
 
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),ext(t-1));

    N_state(n,2,t)= 2-(rand>=pize(1 + ext(t-1)*2 - N_state(n,2,t-1),ext(t)*2-1)/piZ(ext(t-1),ext(t))); 

end

[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));
end


else 
      
for t=2:2000
for n=1:1000
    
    N_state(n,1,t)=a(N_state(n,1,t-1),K_ind(t-1),N_state(n,2,t-1),ext(t-1));
    
   
end

[dev2, K_ind(t)]=min(abs(k_grid(round(mean(N_state(:,1,t))))-K_grid));


end

end

Yg=log(K_grid(K_ind(i_extg(20:end)))');
Xg=[ones(size(i_extg(20:end),2),1),log(K_grid(K_ind(i_extg(20:end)-1))'),ones(size(i_extg(20:end),2),1)*1.2] ;   
Bg=Xg\Yg
b0gp=Bg(1);
b1gp=Bg(2);
b2gp=Bg(3);

Yb=log(K_grid(K_ind(i_extb(20:end)))');
Xb=[ones(size(i_extb(20:end),2),1),log(K_grid(K_ind(i_extb(20:end)-1))'),ones(size(i_extb(20:end),2),1)*0.8]  ;  
Bb=Xb\Yb
b0bp=Bb(1);
b1bp=Bb(2);
b2bp=Bb(3);


dev_b = max(abs([b0g-b0gp b1g-b1gp b2g-b2gp b0b-b0bp b1b-b1bp b2b-b2bp]))

pause(1)
if dev_b<=0.000001 
    break
end


b0g=0.1*b0gp+0.9*b0g; 
b1g=0.1*b1gp+0.9*b1g;
b2g=0.1*b2gp+0.9*b2g;
b0b=0.1*b0bp+0.9*b0b;
b1b=0.1*b1bp+0.9*b1b;
b2b=0.1*b2bp+0.9*b2b;
end
%% 
SSTg = sum(((Yg - mean(Yg)).^2));
SSEg = sum(((Yg' - Bg'*Xg')).^2);
R2_g= 1- SSEg/SSTg

SSTb = sum(((Yb - mean(Yb)).^2));
SSEb = sum(((Yb' - Bb'*Xb')).^2);
R2_b= 1- SSEb/SSTb


%% some graphs
figure

subplot(1,2,1)
hold on
plot(a(:,:,2,1),'r-')
plot(a(:,:,1,1),'b-')

hold off

subplot(1,2,2)
hold on
plot(a(:,:,2,2),'r-')
plot(a(:,:,1,2),'b-')

hold off


subplot(1,2,1)
hold on
plot(a(:,:,1,1),'r-')
plot(a(:,:,1,2),'b-')

hold off

subplot(1,2,2)
hold on
plot(a(:,:,2,1),'r-')
plot(a(:,:,2,2),'b-')

hold off

figure
subplot(1,2,1)
hold on
plot(V0g,'r-')
plot(V0b,'b-')
hold off
subplot(1,2,2)
hold on
plot(V1g,'r-')
plot(V1b,'b-')

hold off

figure
subplot(1,2,1)
hold on
plot(V0g,'r-')
plot(V1g,'b-')
hold off
subplot(1,2,2)
hold on
plot(V0b,'r-')
plot(V1b,'b-')

hold off



figure
for t_ind=1:100
    
   hist(k_grid(reshape(N_state(:,1,t_ind),1,1000)),40)
   legend(num2str(t_ind))
   pause(1)
end
