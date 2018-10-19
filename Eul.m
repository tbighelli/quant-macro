

function F = Eul(x1)
   
   
   theta = .67;
   h = .31;
   beta = .95;   
   delta = .1; 
   T = 100;

   z = 1; 
   k = (z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta);

    z = 2;
knewss = (z*h)*((1-theta)/(beta^(-1)-1+delta))^(1/theta);
   F(1) = (x1(1,1))^(1-theta)*(z*h)^(theta)+(1-delta)*x1(1,1) - x1(2,1) - beta*(k^(1-theta)*(z*h)^theta+(1-delta)*k-x1(1,1))*(1+(1-theta)*(x1(1,1)^(-theta))*((z*h)^(theta))-delta);          
  
   for t = 2:10
   
       F(t) = (x1(t,1))^(1-theta)*(z*h)^(theta)+(1-delta)*x1(t,1) - x1(t+1,1) - beta*(x1(t-1,1)^(1-theta)*(z*h)^theta+(1-delta)*x1(t-1,1)-x1(t,1))*(1+(1-theta)*(x1(t,1)^(-theta))*((z*h)^(theta))-delta);
   
   end
   F(11) = (x1(10-1,1))^(1-theta)*(z*h)^(theta)+(1-delta)*x1(10-1,1) - k - beta*(x1(10-2,1)^(1-theta)*(z*h)^theta+(1-delta)*x1(10-2,1)-x1(10-1,1))*(1+(1-theta)*(x1(10-1,1)^(-theta))*((z*h)^(theta))-delta);
   z=1;
    for t=11:99
        F(t) = (x1(t,1))^(1-theta)*(z*h)^(theta)+(1-delta)*x1(t,1) - x1(t+1,1) - beta*(x1(t-1,1)^(1-theta)*(z*h)^theta+(1-delta)*x1(t-1,1)-x1(t,1))*(1+(1-theta)*(x1(t,1)^(-theta))*((z*h)^(theta))-delta);
    end
F(100) = (x1(T-1,1))^(1-theta)*(z*h)^(theta)+(1-delta)*x1(T-1,1) - k - beta*(x1(T-2,1)^(1-theta)*(z*h)^theta+(1-delta)*x1(T-2,1)-x1(T-1,1))*(1+(1-theta)*(x1(T-1,1)^(-theta))*((z*h)^(theta))-delta);

end
