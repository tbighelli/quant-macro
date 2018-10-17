clear all
clc
aalpha = 0.321;           
bbeta  = 0.988;          
kappa = 5.24;
nu = 2.0;
ddelta = .031;  
GridCapital = 100;
GridLabor= 100;
            
laborSteadyState = 0.2890; 
capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;
vGridCapital = linspace(0.1*capitalSteadyState,2*capitalSteadyState,GridCapital);
vGridLabor=linspace(0.1*laborSteadyState,2*laborSteadyState,GridLabor);
Output           = zeros(GridCapital,GridLabor);
ValueFunction  = zeros(GridCapital,2);
ValueFunctionNew = zeros(GridCapital,2);
PolicyFunction   = zeros(GridCapital,2);
LaborFunction    = zeros(GridCapital,2);
expectedValueFunction = zeros(GridCapital,2);
for Capital=1:GridCapital
    for Labor=1:GridLabor
        Output(Capital,Labor) = (vGridCapital(Capital))^aalpha;
           
       
    end
end
tic 
maxDifference = 100.0;
tolerance = 0.0001;
iteration = 0;
while (maxDifference>tolerance)  
   
    expectedValueFunction = ValueFunction;
    
    if (mod(iteration,10)==0 || iteration ==1)
      
    
        
      
        gridCapitalNextPeriod = 1;
        
        for Capital = 1:GridCapital
                        
            valueHighSoFar = -100.0;
            capitalChoice  = vGridCapital(1);
            gridLabor=1;
            
            for CapitalNextPeriod = gridCapitalNextPeriod:GridCapital
                
                valueHighSoFar1=-100.0;
                
                for Labor = gridLabor:GridLabor
                    
                consumption = Output(Capital,Labor)-vGridCapital(CapitalNextPeriod)+vGridCapital(Capital)*(1-ddelta);
                valueProvisional =log(consumption)-kappa*vGridLabor(Labor)^(1+1/nu)/(1+1/nu)+...
                bbeta*ValueFunction(CapitalNextPeriod);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    gridLabor=Labor;
                    laborChoice=vGridLabor(Labor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(CapitalNextPeriod);
                    gridCapitalNextPeriod = CapitalNextPeriod;
                else
                    break; 
                end    
                  
            
            
            ValueFunctionNew(Capital) = valueHighSoFar;
            PolicyFunction(Capital) = capitalChoice;
            LaborFunction(Capital)=laborChoice;
            GridPolicyFunction(Capital)=gridCapitalNextPeriod;
        end
    end
    
    else
            
                
                 for Capital=1:GridCapital
                
                consumption = (vGridCapital(Capital))^aalpha...
                    (LaborFunction(Capital))^(1-aalpha)-PolicyFunction(Capital)+...
                    vGridCapital(Capital)*(1-ddelta); 
                
                gridCapitalNextPeriod1=GridPolicyFunction(Capital);
                
                ValueFunctionNew(Capital)=log(consumption)-...
                    kappa*mLaborFunction(nCapital)^(1+1/nu)/(1+1/nu)+ ...
                    bbeta*expectedValueFunction(gridCapitalNextPeriod1);
                 end
            
    end  
    
    maxDifference = max(max(abs(ValueFunctionNew-ValueFunction)));
    ValueFunction = ValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end


toc
figure(1)
plot(vGridCapital,ValueFunction)
xlim([vGridCapital(1) vGridCapital(GridCapital)])
xlabel('k');
ylabel('v(k)');


figure
plot(vGridCapital,PolicyFunction)
xlim([vGridCapital(1) vGridCapital(GridCapital)])
xlabel('k');
ylabel('g(k)');

figure
plot(vGridCapital,LaborFunction)
xlim([vGridCapital(1) vGridCapital(GridCapital)])
xlabel('k');
ylabel('h(k)')

