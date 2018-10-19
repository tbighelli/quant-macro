clear all;
close all;
clc;

theta = .67;
h = .31;
beta = .95;   
delta = .1; 
T = 99;

z = 1; 
k = (z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta);
c = (z*h)*(((1-theta)/(beta^(-1)-1+delta))^((1-theta)/theta)-delta*((1-theta)/(beta^(-1)-1+delta))^(1/theta));
y=((z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta))^(1-theta)*(z*h)^theta;
s=y-c;

z = 2; 
kssnew = (z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta);
cssnew = (z*h)*(((1-theta)/(beta^(-1)-1+delta))^((1-theta)/theta)-delta*((1-theta)/(beta^(-1)-1+delta))^(1/theta));
yssnew=((z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta))^(1-theta)*(z*h)^theta;
ssnew=yssnew-cssnew;
fun = @Euler;
x0 = k*ones(100,1);
x = fsolve(fun,x0);

t = [1:100];

% plots 
figure
plot(t,x0);
hold on 
plot(t,kssnew*ones(100,1));
plot(t, x);
ylim([0,3]);
hold off



consumption = ones(99,1);
yy=ones(99,1);
ss=ones(99,1);
for i = 1:99

   consumption(i,1) = (x(i,1)^(1-theta))*((z*h)^(theta))+(1-delta)*x(i,1) - x(i+1,1);
yy(i,1)=(x(i,1)^(1-theta))*((z*h)^(theta));
ss(i,1)=(x(i,1)^(1-theta))*((z*h)^(theta))-consumption(i,1);

end



t = [1:99];
figure
plot(t,c*ones(99,1));
hold on 
plot(t,cssnew*ones(99,1));
plot(t, consumption);
ylim([0,.8]);
hold off

t = [1:99];
figure
 plot(t,y*ones(99,1));
hold on 
plot(t,yssnew.*ones(99,1));
plot(t, yy);
ylim([0,1.5]);
hold off

fun = @Eul;
x10 = k*ones(100,1);
x1 = fsolve(fun,x10);

t = [1:100];

% plots 
figure
plot(t,x10);
hold on 
plot(t,kssnew*ones(100,1));
plot(t, x1);
ylim([0,2.5]);
hold off

consumption = ones(99,1);
yy1=ones(99,1);
ss1=ones(99,1);
for i = 1:99

   consumption1(i,1) = (x1(i,1)^(1-theta))*((z*h)^(theta))+(1-delta)*x1(i,1) - x1(i+1,1);
yy1(i,1)=(x1(i,1)^(1-theta))*((z*h)^(theta));
ss1(i,1)=(x1(i,1)^(1-theta))*((z*h)^(theta))-consumption1(i,1);

end
t = [1:99];
figure
plot(t,c*ones(99,1));
hold on 
plot(t,cssnew*ones(99,1));
plot(t, consumption1);
ylim([0,.8]);
hold off

t = [1:99];
figure
 plot(t,y*ones(99,1));
hold on 
plot(t,yssnew.*ones(99,1));
plot(t, yy1);
ylim([0,1.5]);
hold off