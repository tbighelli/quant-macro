%% Question 2 - no tax
clear all
clc
eta=[1, 1.5, 2.5, 3];
% Random initial endowment
y0=.009 + (.009+.001)*rand(100,1);
for i=1:100
    if .0055<= y0(i,1) <=.0087
    y0(i,1)=.001;
    end
end

r= 0.05;
sigma=3;
k=4;
v=4;
beta=0.99;
tau=0.115;
T0=0;
T1=0;

%% for least-productive workers: eta_y =1 (z=1)
sol1=zeros(100,7);
for n=1:100
   
      i=y0(n,1);
      j=eta(1,1);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);
            %a = x(1);
            %h0 = x(2);
            %h1 = x(3);
            %lambda = x(4);   
           
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.001,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol1(n,1:7)=[u1,i,j,eps];
 
end

%% for second-least-productive workers: eta_y =1.5 (z=2)
sol2=zeros(100,7);
for n=1:100
   
      i=y0(n,1);
      j=eta(1,2);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);
            %a = x(1);
            %h0 = x(2);
            %h1 = x(3);
            %lambda = x(4);   
           
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.001,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol2(n,1:7)=[u1,i,j,eps];
 
end

%% for second-most-productive workers: eta_y =2.5 (z=3)
sol3=zeros(100,7);
for n=1:100
   
      i=y0(n,1);
      j=eta(1,3);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);
            %a = x(1);
            %h0 = x(2);
            %h1 = x(3);
            %lambda = x(4);   
           
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.001,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol3(n,1:7)=[u1,i,j,eps];
 
end

%% for most-productive workers: eta_y =3 (z=4)
sol4=zeros(100,7);
for n=1:100
   
      i=y0(n,1);
      j=eta(1,4);
      s = randi([1, 2], 1);
      e_y=[-.05;.05];
      eps=e_y(s,1);
            %a = x(1);
            %h0 = x(2);
            %h1 = x(3);
            %lambda = x(4);   
           
       f=@(a,h0,h1,lambda)[ ((1-tau)*j*h0 + i + T0 -a)^(-sigma)*(1-tau)*j-k*(h0)^(1/v);
            beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1-tau)*j-k*(h1)^(1/v);
          beta*(((1-tau)*j*h1)+(1+r)*a + T1)^(-sigma)*(1+r)-lambda-((1-tau)*j*h0 + i + T0 -a)^(-sigma);
           ((1-tau)*j*h0 + i + T0 -a) + (1/(1+r))*((1-tau)*(j+eps)*h1 + (1+r)*a + T1) - i - (1+r)*(j+eps)*h1];
       fun=@(u) f(u(1),u(2),u(3),u(4));
      xGuess= ([0.001,0.1,0.1, 1]);
      u1 = fsolve(fun, xGuess) ;
      sol4(n,1:7)=[u1,i,j,eps];
 
end

%% 

sol=[sol1;sol2;sol3;sol4]; %sol=[a, h0, h1, lambda, y0, eta_y, eps]

%c0=(1-tau)*j*h0+i+T0-a
%c1=(1-tau)*(j+eps)*h1+(1+r)*a+T1
%savings rate=a/(i+j*h0*(1-tau))
c0=(1-tau).*sol(:,6).*sol(:,2)+sol(:,5)+T0-sol(:,1);
c1=(1-tau).*(sol(:,6)+sol(:,7)).*sol(:,3)+(1+r).*sol(:,1)+T1;
sv=sol(:,1)./(sol(:,5)+sol(:,6).*sol(:,2).*(1-tau));
%sv1=sol1(:,1)./(sol1(:,5)+sol1(:,6).*sol1(:,2).*(1-tau)); 
c01=(1-tau).*sol1(:,6).*sol1(:,2)+sol1(:,5)+T0-sol1(:,1);
c02=(1-tau).*sol2(:,6).*sol2(:,2)+sol2(:,5)+T0-sol2(:,1);
c03=(1-tau).*sol3(:,6).*sol3(:,2)+sol3(:,5)+T0-sol3(:,1);
c04=(1-tau).*sol4(:,6).*sol4(:,2)+sol4(:,5)+T0-sol4(:,1);
sv1=sol1(:,1)./(sol1(:,5)+sol1(:,6).*sol1(:,2).*(1-tau));
sv2=sol2(:,1)./(sol2(:,5)+sol2(:,6).*sol2(:,2).*(1-tau));
sv3=sol3(:,1)./(sol3(:,5)+sol3(:,6).*sol3(:,2).*(1-tau));
sv4=sol4(:,1)./(sol4(:,5)+sol4(:,6).*sol4(:,2).*(1-tau));
ssol1=sort(sol1(:,1));
ssol2=sort(sol2(:,1));
ssol3=sort(sol3(:,1));
ssol4=sort(sol4(:,1));
c11=(1-tau).*(sol1(:,6)+sol1(:,7)).*sol1(:,3)+(1+r).*sol1(:,1)+T1;
c12=(1-tau).*(sol2(:,6)+sol2(:,7)).*sol2(:,3)+(1+r).*sol2(:,1)+T1;
c13=(1-tau).*(sol3(:,6)+sol3(:,7)).*sol3(:,3)+(1+r).*sol3(:,1)+T1;
c14=(1-tau).*(sol4(:,6)+sol4(:,7)).*sol4(:,3)+(1+r).*sol4(:,1)+T1;
figure 
hold on
plot(y0,sol1(:,1))
plot(y0,sol2(:,2))
plot(y0,sol3(:,3))
plot(y0,sol4(:,4))
legend('eta =1','eta=1.5','eta=2.5','eta=3')
hold off
figure 
hold on
plot(y0,c01)
plot(y0,c02)
plot(y0,c03)
plot(y0,c04)
legend('eta =1','eta=1.5','eta=2.5','eta=3')
hold off
figure 
hold on
plot(y0,c01)
plot(y0,c02)
plot(y0,c03)
plot(y0,c04)
legend('eta =1','eta=1.5','eta=2.5','eta=3')
hold off
figure 
hold on
plot(y0,sol1(:,2))
plot(y0,sol2(:,2))
plot(y0,sol3(:,2))
plot(y0,sol4(:,2))
hold off
figure
hold on
bbb1=sort(c01);
bbb2=sort(c02);
bbb3=sort(c03);
bbb4=sort(c04);
hold off
ssd=sort(sv1)
figure 
hold on
plot(y0,sol1(:,2))
plot(y0,sol2(:,2))
plot(y0,sol3(:,2))
plot(y0,sol4(:,2))
legend('eta =1','eta=1.5','eta=2.5','eta=3')
hold off
figure 
hold on
plot(y0,sol1(:,3))
plot(y0,sol2(:,3))
plot(y0,sol3(:,3))
plot(y0,sol4(:,3))
legend('eta =1','eta=1.5','eta=2.5','eta=3')
hold off