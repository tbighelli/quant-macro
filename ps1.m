% Tommaso Bighelli - Quantitative Macroeconomics PS I
%% EX  1
clear all
clc
x0=1;

t20=sym([21,1]);
syms y(x);
y(x)=x.^0.321;
for i=1:21
    t20(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-1)^(i-1);
    
end
t20(1)=y(x0);
T20=sum(t20);
t5=sym([6,1]);
for i=1:6
    t5(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-1)^(i-1);
    
end
t5(1)=y(x0);
T5=sum(t5);
t2=sym([3,1]);
for i=1:3
    t2(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-1)^(i-1);
    
end
t2(1)=y(x0);
T2=sum(t2);
t1=sym([2,1]);
for i=1:2
    t1(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-1)^(i-1);
    
end
t1(1)=y(x0);
T1=sum(t1);
figure(1)
hold on

axis([0 4 -2 7])
fplot(T20)
fplot(T5)
fplot(T2)
fplot(T1)
fplot(y(x))
xlabel('x');
ylabel('y(x)');
legend('T20(x)','T5(x)','T2(x)','T1(x)','y(x)');
%% EX 2
clear all
clc
x0=2;
t20=sym([21,1]);
syms y(x);
y(x)=(x+abs(x))/2;
for i=1:21
    t20(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-x0)^(i-1);
    
end
t20(1)=y(x0);
T20=sum(t20);

t5=sym([6,1]);
syms y(x);
y(x)=(x+abs(x))/2;
for i=1:5
    t5(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-x0)^(i-1);
    
end
t5(1)=y(x0);
T5=sum(t5);

t2=sym([3,1]);
syms y(x);
y(x)=(x+abs(x))/2;
for i=1:3
    t2(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-x0)^(i-1);
    
end
t2(1)=y(x0);
T2=sum(t2);

t1=sym([2,1]);
syms y(x);
y(x)=(x+abs(x))/2;
for i=1:2
    t1(i)=(feval(diff(y,i-1),x0)/factorial(i-1))*(x-x0)^(i-1);
    
end
t1(1)=y(x0);
T1=sum(t1);

figure(1)
hold on

axis([-2 6 -2 6])
fplot(T20)
fplot(T5)
fplot(T2)
fplot(T1)
fplot(y(x))
xlabel('x');
ylabel('y(x)');
legend('T20(x)','T5(x)','T2(x)','T1(x)','y(x)');



%%  EX 3 - runge function
clear all
clc

x=sym([10,1]);
x1=linspace(-1,1,4);
yy=1./(1+25*x1.^2);

syms y(x);
y=1./(1+25*x.^2);

p=polyfit(x1,yy,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
x2=linspace(-1,1,6);
yy2=1./(1+25*x2.^2);
%yyy=1./(1+25*x.^2);
p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
x3=linspace(-1,1,11);
yy3=1./(1+25*x3.^2);
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
    
   f3=sum(f3);

end
xq=-1:0.01:1;
for i=1:201
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);
 figure
 hold on
axis([-1 1 -0.2 2]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Evenly spaced interpolation nodes, cubic polynomials, monomials of order 5 and 10')
legend('Runge(x)','cubic','g5)','g10');

figure
hold on
axis([-1 1 0 2]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('Runge function approximation errors: evenly spaced and monomials');
legend('cubic error','g5 error','g10 error');



%% EX 3 - ramp function
clear all
clc

x=sym([10,1]);
x1=linspace(-1,1,4);
yy=(x1+abs(x1))/2;

syms y(x);
y=(x+abs(x))/2;

p=polyfit(x1,yy,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
x2=linspace(-1,1,6);
yy2=(x2+abs(x2))/2;

p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
x3=linspace(-1,1,11);
yy3=(x3+abs(x3))/2;
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=-1:0.01:1;
for i=1:201
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

 figure
 hold on
axis([-1 1 -0.5 2]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Evenly spaced interpolation nodes, cubic polynomials, monomials of order 5 and 10')
legend('ramp(x)','cubic','g5','g10');

figure
hold on
axis([-1 1 0 .25]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('ramp function approximation errors: evenly spaced and monomials');
legend('cubic error','g5 error','g10 error');

%% EX 3 - exponential

clear all
clc

x=sym([10,1]);
x1=linspace(-1,1,4);
yy=exp(1./x1);

syms y(x);
y=exp(1./x);

p=polyfit(x1,yy,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
x2=linspace(-1,1,6);
yy2=exp(1./x2);

p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
x3=linspace(-1,1,11);
yy3=exp(1./x3);
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=[-1:0.011:1];
for i=1:182
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

 figure
 hold on
axis([-5 8 -5 50]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Evenly spaced interpolation nodes, cubic polynomials, monomials of order 5 and 10')
legend('exp','cubic','g5','g10');
figure
hold on
axis([-1 1 0 30]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('exponential function approximation errors: evenly spaced and monomials');
legend('cubic error','g5 error','g10 error');


%%
%%  EX 3 -CHEB AND MONOMIALS runge function (cancellare)
clear all
clc

x=sym([10,1]);
x1=linspace(-1,1,4);
yy=1./(1+25*x1.^2);

syms y(x);
y=1./(1+25*x.^2);

p=polyfit(x1,yy,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
x2=linspace(-1,1,6);
yy2=1./(1+25*x2.^2);
%yyy=1./(1+25*x.^2);
p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
x3=linspace(-1,1,11);
yy3=1./(1+25*x3.^2);
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=-1:0.01:1;
for i=1:201
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

 figure
 hold on
axis([-1 1 -0.2 2]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Chebyshev interpolation nodes, cubic polynomials, monomials of order 5 and 10')

legend('Runge(x)','cubic','g5)','g10');
hold off
figure
hold on
axis([-1 1 0 2]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('runge function approximation errors: evenly spaced and monomials');
legend('cubic error','g5 error','g10 error');

%% EX 3  CHEBYSHEV INTERPOLATION NODES AND  mon
clear all
clc

for i=1:10
    x(i)=cos(((2*i-1)/20)*pi);
end
%x=sym([10,1]);
for i=1:10
    x1(i)=cos(((2*i-1)/20)*pi);
end

yy=1./(1+25*x1.^2);

syms y(x);
y=1./(1+25*x.^2);
%cubic
for i=1:4   
    xc(i)=cos(((2*i-1)/8)*pi);
end
yc=1./(1+25*xc.^2);
p=polyfit(xc,yc,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
for i=1:6
    x2(i)=cos(((2*i-1)/12)*pi);
end
yy2=1./(1+25*x2.^2);
%yyy=1./(1+25*x.^2);
p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
for i=1:11
    x3(i)=cos(((2*i-1)/22)*pi);
end
yy3=1./(1+25*x3.^2);
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=-1:0.01:1;
for i=1:201
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);
 figure
 hold on
axis([-1 1 -0.2 1]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Chebyshev interpolation nodes, cubic polynomials, polynomials of order 5 and 10')

legend('Runge(x)','cubic','g(5)','g(10)');

legend('Runge(x)','cubic','g5)','g10');
hold off
figure
hold on
axis([-1 1 0 1]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('runge function approximation errors: evenly spaced and monomials');
legend('cubic error','g5 error','g10 error');
%% EX 3 CHEBYSHEV INTERP NODES, AND mon, RAMP
clear all
clc

for i=1:10
    x(i)=cos(((2*i-1)/20)*pi);
end
%x=sym([10,1]);
for i=1:10
    x1(i)=cos(((2*i-1)/20)*pi);
end

yy=(x1+abs(x1))/2;

syms y(x);
y=(x+abs(x))/2;
%cubic
for i=1:4   
    xc(i)=cos(((2*i-1)/8)*pi);
end
yc=(xc+abs(xc))/2;
p=polyfit(xc,yc,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
for i=1:6
    x2(i)=cos(((2*i-1)/12)*pi);
end
yy2=(x2+abs(x2))/2;
%yyy=1./(1+25*x.^2);
p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
for i=1:11
    x3(i)=cos(((2*i-1)/22)*pi);
end
yy3=(x3+abs(x3))/2;
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=-1:0.01:1;
for i=1:201
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);
 figure
 hold on
axis([-1 1 -0.5 2]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Chebyshev interpolation nodes, cubic polynomials, monomials of order 5 and 10')

legend('Ramp(x)','cubic','g(5)','g(10)');

legend('Ramp(x)','cubic','g5)','g10');
hold off
figure
hold on
axis([-1 1 0 0.3]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('ramp function approximation errors: chebyshev nodes and monomials');
legend('cubic error','g5 error','g10 error');
%%
%% EX 3 CHEBYSHEV INTERP NODES, AND mon, exp
clear all
clc

for i=1:10
    x(i)=cos(((2*i-1)/20)*pi);
end
%x=sym([10,1]);
for i=1:10
    x1(i)=cos(((2*i-1)/20)*pi);
end

yy=exp(1./x1);

syms y(x);
y=exp(1./x);
%cubic
for i=1:4   
    xc(i)=cos(((2*i-1)/8)*pi);
end
yc=exp(1./xc);
p=polyfit(xc,yc,3);
syms f theta
f =  p(1)*theta^3 + p(2)*theta^2 + p(3)*theta + p(4);
% monomials 5
for i=1:6
    x2(i)=cos(((2*i-1)/12)*pi);
end
yy2=exp(1./x2);
%yyy=1./(1+25*x.^2);
p2=polyfit(x2,yy2,5);
syms f2 theta2
for i=1:6
    f2(i)=p2(i)*theta2^(6-i);
   f2=sum(f2);
end

% mon 10
for i=1:11
    x3(i)=cos(((2*i-1)/22)*pi);
end
yy3=exp(1./x3);
p3=polyfit(x3,yy3,10);
syms f3 theta3
for i=1:11
    f3(i)=p3(i)*theta3^(11-i);
   f3=sum(f3);
end
xq=[-1:0.011:1];
for i=1:182
err(i)=subs(f3,theta3,xq(i));
ery(i)=subs(y,x,xq(i));
er2(i)=subs(f2,theta2,xq(i));
er1(i)=subs(f,theta,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

 figure
 hold on
axis([-5 8 -5 50]);
fplot(y)
fplot(f)
fplot(f2)
fplot(f3)
title('Chebyshev interpolation nodes, cubic polynomials, monomials of order 5 and 10')
legend('exp','cubic','g5','g10');
figure
hold on
axis([-1 1 0 30]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('exponential function approximation errors: evenly chebyshev nodes and monomials');
legend('cubic error','g5 error','g10 error');

%% cheb nodes cheb poly runge
clear all
clc


m=4;
y3=zeros(m,1);
psi=zeros(1,m);
phi3=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi3(i,:)=psi;
    y3(i)=1/((1+25*x(i).^2));
end

theta3=phi3^(-1)*y3;

m=6;
y5=zeros(m,1);
psi=zeros(1,m);
phi5=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi5(i,:)=psi;
    y5(i)=1/((1+25*x(i).^2));
end
theta5=phi5^(-1)*y5;


m=11;
y10=zeros(m,1);
psi=zeros(1,m);
phi10=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi10(i,:)=psi;
    y10(i)=1/((1+25*x(i).^2));
end

theta10=phi10^(-1)*y10;

syms f(x) ff(x) ff5(x) ff10(x)
f(x)=1/((1+25*x.^2));
ff3(x)=[1 x chebyshevT([2,3],x)]*theta3;
ff5(x)=[1 x chebyshevT([2,3,4,5],x)]*theta5;
ff10(x)=[1 x chebyshevT([2,3,4,5,6,7,8,9,10],x)]*theta10;

xq=-1:0.01:1;
for i=1:201
err(i)=subs(ff10,x,xq(i));
ery(i)=subs(f,x,xq(i));
er2(i)=subs(ff5,x,xq(i));
er1(i)=subs(ff3,x,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

figure
hold on
axis([-1 1 -0.2 1]);
fplot(f,[-1 1]);
fplot(ff3,[-1 1]);
fplot(ff5,[-1 1]);
fplot(ff10,[-1 1]);
title('Chebyshev interpolation nodes, chebyshev polynomials  of order 3, 5 and 10')
legend('runge','ch3','ch5','ch10');

figure
hold on
axis([-1 1 0 1]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('runge function approximation errors: chebyshev nodes and chebyshev polys');
legend('ch3 error','ch5 error','ch10 error');
%%
%% cheb nodes cheb poly ramp
clear all
clc


m=4;
y3=zeros(m,1);
psi=zeros(1,m);
phi3=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2.*x(i).*psi(j)-psi(j-1);
    end
    phi3(i,:)=psi;
    y3(i)=(x(i)+abs(x(i)))./2;
end

theta3=phi3^(-1)*y3;

m=6;
y5=zeros(m,1);
psi=zeros(1,m);
phi5=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi5(i,:)=psi;
    y5(i)=(x(i)+abs(x(i)))./2;
end
theta5=phi5^(-1)*y5;


m=11;
y10=zeros(m,1);
psi=zeros(1,m);
phi10=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi10(i,:)=psi;
    y10(i)=(x(i)+abs(x(i)))./2;
end

theta10=phi10^(-1)*y10;

syms f(x) ff3(x) ff5(x) ff10(x)
%f(x)=(x(i)+abs(x(i)))./2;
f(x)=(x+abs(x))/2;
ff3(x)=[1 x chebyshevT([2,3],x)]*theta3;
ff5(x)=[1 x chebyshevT([2,3,4,5],x)]*theta5;
ff10(x)=[1 x chebyshevT([2,3,4,5,6,7,8,9,10],x)]*theta10;

xq=-1:0.01:1;
for i=1:201
err(i)=subs(ff10,x,xq(i));
ery(i)=subs(f,x,xq(i));
er2(i)=subs(ff5,x,xq(i));
er1(i)=subs(ff3,x,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

figure
hold on
axis([-1 1 -0.2 1]);
fplot(f,[-1 1]);
fplot(ff3,[-1 1]);
fplot(ff5,[-1 1]);
fplot(ff10,[-1 1]);
title('Chebyshev interpolation nodes, chebyshev polynomials  of order 3, 5 and 10')
legend('ramp','ch3','ch5','ch10');

figure
hold on
axis([-1 1 0 0.3]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('ramp function approximation errors: chebyshev nodes and chebyshev polys');
legend('ch3 error','ch5 error','ch10 error');
%%
%% cheb nodes cheb poly exp
clear all
clc


m=4;
y3=zeros(m,1);
psi=zeros(1,m);
phi3=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2.*x(i).*psi(j)-psi(j-1);
    end
    phi3(i,:)=psi;
    y3(i)=exp(1./x(i));
end

theta3=phi3^(-1)*y3;

m=6;
y5=zeros(m,1);
psi=zeros(1,m);
phi5=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi5(i,:)=psi;
    y5(i)=exp(1./x(i));
end
theta5=phi5^(-1)*y5;


m=11;
y10=zeros(m,1);
psi=zeros(1,m);
phi10=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=-cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi10(i,:)=psi;
    y10(i)=exp(1./x(i));
end

theta10=phi10^(-1)*y10;

syms f(x) ff3(x) ff5(x) ff10(x)

f(x)=exp(1./x);
ff3(x)=[1 x chebyshevT([2,3],x)]*theta3;
ff5(x)=[1 x chebyshevT([2,3,4,5],x)]*theta5;
ff10(x)=[1 x chebyshevT([2,3,4,5,6,7,8,9,10],x)]*theta10;


xq=[-1:0.011:1];
for i=1:182
err(i)=subs(ff10,x,xq(i));
ery(i)=subs(f,x,xq(i));
er2(i)=subs(ff5,x,xq(i));
er1(i)=subs(ff3,x,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);

figure
hold on
axis([-1 1 -0.2 100]);
fplot(f,[-1 1]);
fplot(ff3,[-1 1]);
fplot(ff5,[-1 1]);
fplot(ff10,[-1 1]);
title('Chebyshev interpolation nodes, chebyshev polynomials  of order 3, 5 and 10')
legend('exp','ch3','ch5','ch10');

figure
hold on
axis([-1 1 0 50]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);
title('exp function approximation errors: chebyshev nodes and chebyshev polys');
legend('ch3 error','ch5 error','ch10 error');
%% EX 4 PROBABILITY
clear all
clc
alpha=1;
rho1=1/0.2;
rho2=1/100;

m=4;
y3=zeros(m,1);
psi=zeros(1,m);
phi3=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi3(i,:)=psi;
    y3(i)=exp(-alpha*x(i))/(rho1+rho2*exp(-alpha*x(i)));
end

theta3=phi3^(-1)*y3;

m=6;
y5=zeros(m,1);
psi=zeros(1,m);
phi5=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi5(i,:)=psi;
    y5(i)=exp(-alpha*x(i))/(rho1+rho2*exp(-alpha*x(i)));
end
theta5=phi5^(-1)*y5;


m=11;
y10=zeros(m,1);
psi=zeros(1,m);
phi10=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi10(i,:)=psi;
    y10(i)=exp(-alpha*x(i))/(rho1+rho2*exp(-alpha*x(i)));
end

theta10=phi10^(-1)*y10;

syms f(x) ff(x) ff5(x) ff10(x)
f(x)=exp(-alpha*x)/(rho1+rho2*exp(-alpha*x));
ff3(x)=[1 x chebyshevT([2,3],x)]*theta3;
ff5(x)=[1 x chebyshevT([2,3,4,5],x)]*theta5;
ff10(x)=[1 x chebyshevT([2,3,4,5,6,7,8,9,10],x)]*theta10;

xq=-1:0.05:10;
for i=1:221
err(i)=subs(ff10,x,xq(i));
ery(i)=subs(f,x,xq(i));
er2(i)=subs(ff5,x,xq(i));
er1(i)=subs(ff3,x,xq(i));
end
er3=double(err);
ser=double(ery);
err3=abs(er3-ser);
er2=double(er2);
er1=double(er1);
err2=abs(er2-ser);
err1=abs(er1-ser);




brho1=1/0.25;


m=4;
by3=zeros(m,1);
psi=zeros(1,m);
phi3=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi3(i,:)=psi;
    by3(i)=exp(-alpha*x(i))/(brho1+rho2*exp(-alpha*x(i)));
end

btheta3=phi3^(-1)*by3;

m=6;
by5=zeros(m,1);
psi=zeros(1,m);
phi5=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi5(i,:)=psi;
    by5(i)=exp(-alpha*x(i))/(brho1+rho2*exp(-alpha*x(i)));
end
btheta5=phi5^(-1)*by5;


m=11;
by10=zeros(m,1);
psi=zeros(1,m);
phi10=zeros(m,m);
x =zeros(m,1);
psi(1)=1;


for i=1:m
    x(i)=0.5*10+0.5*10*cos(((2*i-1)*pi)/(2*m));
    psi(2)=x(i);
    for j=2:m-1
        psi(j+1)=2*x(i)*psi(j)-psi(j-1);
    end
    phi10(i,:)=psi;
    by10(i)=exp(-alpha*x(i))/(brho1+rho2*exp(-alpha*x(i)));
end

btheta10=phi10^(-1)*by10;

syms f(x) ff(x) ff5(x) ff10(x)
bf(x)=exp(-alpha*x)/(brho1+rho2*exp(-alpha*x));
bff3(x)=[1 x chebyshevT([2,3],x)]*theta3;
bff5(x)=[1 x chebyshevT([2,3,4,5],x)]*theta5;
bff10(x)=[1 x chebyshevT([2,3,4,5,6,7,8,9,10],x)]*theta10;

xq=-1:0.05:10;
for i=1:221
berr(i)=subs(bff10,x,xq(i));
bery(i)=subs(bf,x,xq(i));
ber2(i)=subs(bff5,x,xq(i));
ber1(i)=subs(bff3,x,xq(i));
end
ber3=double(berr);
bser=double(bery);
berr3=abs(ber3-bser);
ber2=double(ber2);
ber1=double(ber1);
berr2=abs(ber2-bser);
berr1=abs(ber1-bser);


figure
hold on
fplot(f,[0 10]);
fplot(ff3,[0 10]);
fplot(ff5,[0 10]);
fplot(ff10,[0 10]);
fplot(bf,[0 10]);
fplot(bff3,[0 10]);
fplot(bff5,[0 10]);
fplot(bff10,[0 10]);
title('probabiliy function approximation: chebyshev nodes and chebyshev polys');
legend('p(x)-1','ch 3-1','ch 5-1','ch 10-1','p(x)-2','ch 3-2','ch 5-2','ch 10-2');
figure
hold on
axis([-1 1 0 0.3]);
plot(xq,err1);
plot(xq,err2);

plot(xq,err3);

plot(xq,berr1);
plot(xq,berr2);

plot(xq,berr3);
title('probability function approximation errors');
legend('ch3 error','ch5 error','ch10 error','ch3-2 error ','ch5-2 error','ch10-2 error');

%% q2
clear all
clc

alpha=0.5;
s=0.25;
m=20;
H=4;
k=zeros(m,1);
h=zeros(m,1);
phi1=zeros(H,H);
phi2=zeros(H,H);
psi1=zeros(1,H);
psi2=zeros(1,H);
psi1(1) =1;
psi2(1) =1;
for i=1:m
    k(i)=cos(((2*i-1)*pi)/(2*m));
    k(i)=(k(i)+1)*5;
    h(i)=cos(((2*i-1)*pi)/(2*m));
    h(i)=(h(i)+1)*5;
    psi1(2)=k(i);
    psi2(2)=h(i);
    for j=2:H-1
        psi1(j+1)=2*k(i)*psi1(j)-psi1(j-1);
        psi2(j+1)=2*h(i)*psi2(j)-psi2(j-1);
    end
    phi1(i,:)=psi1;
    phi2(i,:)=psi2;
end


phi=kron(phi2,phi1);
K=repmat(k,m,1);
H=repmat(h,1,m);
H=H';
H=H(:);
x=[K H];
y=zeros(m*m,1);
for i=1:m*m
    y(i)=((1-alpha)*x(i,1)^((s-1)/s)+alpha*x(i,2)^((s-1)/s))^(s/(s-1));
end
theta=(phi'*phi)^(-1)*phi'*y; 

syms f(k,h) ff(k,h);
f(k,h)=((1-alpha)*k^((s-1)/s)+alpha*h^((s-1)/s))^(s/(s-1));
ff(k,h)=kron([1 h chebyshevT([2,3],h)],[1 k chebyshevT([2,3],k)])*theta;
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold off
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');

hold off

syms err(k,h);
err(k,h)=abs(f(k,h)-ff(k,h));
figure 
fsurf(err,[0,10]);
fcontour(err,[0,10]);
pc=prctile(y,[5 10 25 50 75 90 95]);
figure
fcontour(f,[0,10],'LevelList',pc);
figure
fcontour(ff,[0,10],'LevelList',pc);
%%
%% q2 b
clear all
clc

alpha=0.5;
s=0.9999;
m=20;
H=4;
k=zeros(m,1);
h=zeros(m,1);
phi1=zeros(H,H);
phi2=zeros(H,H);
psi1=zeros(1,H);
psi2=zeros(1,H);
psi1(1) =1;
psi2(1) =1;
for i=1:m
    k(i)=cos(((2*i-1)*pi)/(2*m));
    k(i)=(k(i)+1)*5;
    h(i)=cos(((2*i-1)*pi)/(2*m));
    h(i)=(h(i)+1)*5;
    psi1(2)=k(i);
    psi2(2)=h(i);
    for j=2:H-1
        psi1(j+1)=2*k(i)*psi1(j)-psi1(j-1);
        psi2(j+1)=2*h(i)*psi2(j)-psi2(j-1);
    end
    phi1(i,:)=psi1;
    phi2(i,:)=psi2;
end


phi=kron(phi2,phi1);
K=repmat(k,m,1);
H=repmat(h,1,m);
H=H';
H=H(:);
x=[K H];
y=zeros(m*m,1);
for i=1:m*m
    y(i)=((1-alpha)*x(i,1)^((s-1)/s)+alpha*x(i,2)^((s-1)/s))^(s/(s-1));
end
theta=(phi'*phi)^(-1)*phi'*y; 

syms f(k,h) ff(k,h);
f(k,h)=((1-alpha)*k^((s-1)/s)+alpha*h^((s-1)/s))^(s/(s-1));
ff(k,h)=kron([1 h chebyshevT([2,3],h)],[1 k chebyshevT([2,3],k)])*theta;
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold off
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');

hold off

syms err(k,h);
err(k,h)=abs(f(k,h)-ff(k,h));
figure 
fsurf(err,[0,10]);
fcontour(err,[0,10]);
pc=prctile(y,[5 10 25 50 75 90 95]);
figure
fcontour(f,[0,10],'LevelList',pc);
figure
fcontour(ff,[0,10],'LevelList',pc);
%% q2 c
clear all
clc

alpha=0.5;
s=5;
m=20;
H=4;
k=zeros(m,1);
h=zeros(m,1);
phi1=zeros(H,H);
phi2=zeros(H,H);
psi1=zeros(1,H);
psi2=zeros(1,H);
psi1(1) =1;
psi2(1) =1;
for i=1:m
    k(i)=cos(((2*i-1)*pi)/(2*m));
    k(i)=(k(i)+1)*5;
    h(i)=cos(((2*i-1)*pi)/(2*m));
    h(i)=(h(i)+1)*5;
    psi1(2)=k(i);
    psi2(2)=h(i);
    for j=2:H-1
        psi1(j+1)=2*k(i)*psi1(j)-psi1(j-1);
        psi2(j+1)=2*h(i)*psi2(j)-psi2(j-1);
    end
    phi1(i,:)=psi1;
    phi2(i,:)=psi2;
end


phi=kron(phi2,phi1);
K=repmat(k,m,1);
H=repmat(h,1,m);
H=H';
H=H(:);
x=[K H];
y=zeros(m*m,1);
for i=1:m*m
    y(i)=((1-alpha)*x(i,1)^((s-1)/s)+alpha*x(i,2)^((s-1)/s))^(s/(s-1));
end
theta=(phi'*phi)^(-1)*phi'*y; 

syms f(k,h) ff(k,h);
f(k,h)=((1-alpha)*k^((s-1)/s)+alpha*h^((s-1)/s))^(s/(s-1));
ff(k,h)=kron([1 h chebyshevT([2,3],h)],[1 k chebyshevT([2,3],k)])*theta;
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold off
figure
fsurf(f,[0,10],'ShowContours','on');
xlabel('h');
ylabel('k');
hold on
fsurf(ff,[0,10],'ShowContours','on');

hold off

syms err(k,h);
err(k,h)=abs(f(k,h)-ff(k,h));
figure 
fsurf(err,[0,10]);
fcontour(err,[0,10]);
pc=prctile(y,[5 10 25 50 75 90 95]);
figure
fcontour(f,[0,10],'LevelList',pc);
figure
fcontour(ff,[0,10],'LevelList',pc);