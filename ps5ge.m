clear
clc

t1 = clock;
sigma  = 2;                 
beta   = 1/1.06;            % 
prob   = [ .8 .2; .4 .6];
delta  = 0.97;            
A      = 1.00;            
alpha  = 0.25;           
theta  = 0.4;            
Kstart = 11.5119;            
g      = 0.2;          
r=0.04;
mpk=0.07;
kap=-1.56:0.025:16;
kapp=-1.56:0.025:16;
maxkap = 16;                        
inckap = 0.025;                   
nkap=703
nkap   = round((maxkap+1.56)/inckap+1);  
D = zeros(length(prob));
[ev,ed] = eig(prob);
[emax,inmax] = max(diag(ed));
D(inmax,inmax) = emax;
pinf = ev*D*inv(ev);
pempl = pinf(inmax,inmax);
N = 1;
liter   = 1;
maxiter = 5000;
toler   = 0.000001;
metric  = 10;
K = Kstart;
while  (metric > toler) & (liter <= maxiter);
   wage = (1-alpha) * A * K^(alpha)   * N^(-alpha);
   rent=r;
   util1=-10000*ones(nkap,nkap);  
   util2=-10000*ones(nkap,nkap);  
   for i=1:nkap;
             
         
         
         for j=1:nkap; 
               
               if kapp(j)>0
               cons1 = wage+theta + (rent + 1)*kap(i) - kapp(j); 
               if cons1 > 0;
                  util1(j,i)=((cons1)^(1-sigma)-1)/(1-sigma);
               end;
               else
                 util1(j,i)=-100;  
               end
               if kapp(j)>0
	       cons2 = wage-theta + (rent + 1)*kap(i) - kapp(j);
               if cons2 > 0;
                  util2(j,i)=((cons2)^(1-sigma)-1)/(1-sigma);
               end;
               else
                   util2(j,i)=-100;
               end
         end;
      
   end

   
   v       = zeros(nkap,2);
   decis   = zeros(nkap,2);
   test    = 10;
   [rs,cs] = size(util1);
   
   while test ~= 0;
       for i=1:cs;
           r1(:,i)=util1(:,i)+beta*(prob(1,1)*v(:,1)+ prob(1,2)*v(:,2));
           r2(:,i)=util2(:,i)+beta*(prob(2,1)*v(:,1)+ prob(2,2)*v(:,2));
       end;

       [tv1,tdecis1]=max(r1);
       [tv2,tdecis2]=max(r2);
       tdecis=[tdecis1' tdecis2'];
       tv=[tv1' tv2'];

       test=max(any(tdecis-decis));
       v=tv;
       decis=tdecis;

   end
   decis=(decis-1)*inckap;
   
   g2=sparse(cs,cs);
   g1=sparse(cs,cs);
   for i=1:cs
       g1(i,tdecis1(i))=1;
       g2(i,tdecis2(i))=1;
   end
   trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
   trans=trans';
   probst = (1/(2*nkap))*ones(2*nkap,1);
   test=1;
   while test > 10^(-8);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1;
   end
   
   kk=decis(:);
   meanK=probst'*kk;
   
   lambda=zeros(cs,2);
   lambda(:)=probst;
  
   [v1,d1]=eig(prob');
   [dmax,imax]=max(diag(d1));
   probst1=v1(:,imax);
   ss=sum(probst1);
   probst1=probst1/ss;
   probk=sum(lambda');     
   probk=probk';
 
   Kold = K;
   Knew = g*meanK + (1-g)*Kold;
   metric = abs((Kold-meanK)/Kold);
   K = Knew;
   disp([ liter metric meanK Kold ]);
   liter = liter+1;
   end


k = Kold;                  
n = 100;                
s0 = 1;                 
hist = zeros(n-1,2);
cons = zeros(n-1,1);
invest = zeros(n-1,1);
grid = [ (0:inckap:maxkap+1.56)' ];  
[chain,state] = markov(prob,n,s0);
for i = 1:n-1;
    hist(i,:) = [ k chain(i) ];
    I1 = round(k/inckap) ;
    I2 = round(k/inckap) + 1;
    if I1 == 0;
       I1=1;
       
    end;
    if I2 > nkap;
       I2 = nkap;
     
    end;
    weight = (grid(I2,1) - k)/inckap; 
    kprime = weight*(decis(I1,chain(i))) +  (1-weight)*(decis(I2,chain(i)));
    if chain(i) == 1;
        cons(i)=wage+theta + (rent + 1)*k - kprime;
       %cons(i) = wage + (rent + delta)*k - kprime;
    elseif chain(i) == 2;
        cons(i)=wage+theta + (rent + 1)*k - kprime;
       %cons(i) = wage*theta + (rent + delta)*k - kprime;
    else;
      
      chain
    end;
    k = kprime;
    invest(i) = kprime;
end;
figure
hold on




figure
income =  [ (rent*grid + wage-theta)  (rent*grid + wage+theta) ]  ; 
[ pinc, index ] = sort(income(:));
plambda = lambda(:);
plot(pinc,plambda(index));
%title('MODEL 2: INCOME DISTRIBUTION');
xlabel('income');
ylabel('% of agents');
print distmod2



figure
hold on
plot(probk);

xlabel('a');
ylabel('% of agents');


figure
hold on 
plot(tdecis)
refline(1,0);
xlabel('a');
ylabel('g(a,y');
% offerta aggr risparmio
for i=1:703
    EA(i)=probk(i)*grid(i);
    EA=sum(EA)
end
dd=K-EA;
figure
hold on
plot(v)
xlabel('a');
ylabel('v(a,y)');
hold off
