aalpha = 0.321;           
bbeta  = 0.988;          
ddelta = .013;        
kappa = 5.24;
nu = 2.0;
GridCapital = 1000;
GridLabor= 1000;
vProductivity = [1.01; 1/1.01;]';
Transition   = [1/3, 0;
                 0, 2/3];
laborSteadyState = 0.2890; 
capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;
vGridCapital = linspace(0.01*capitalSteadyState,2*capitalSteadyState,GridCapital);
vGridLabor=linspace(0.01*laborSteadyState,2*laborSteadyState,GridLabor);
GridProductivity = length(vProductivity);
Output           = zeros(GridCapital,GridLabor,GridProductivity);
ValueFunction  = zeros(GridCapital,GridProductivity);
ValueFunctionNew = zeros(GridCapital,GridProductivity);
PolicyFunction   = zeros(GridCapital,GridProductivity);
LaborFunction    = zeros(GridCapital,GridProductivity);
expectedValueFunction = zeros(GridCapital,GridProductivity);
for Capital=1:GridCapital
    for Labor=1:GridLabor
        for Productivity=1:GridProductivity
            Output(Capital,Labor,Productivity) = (vGridCapital(Capital))^aalpha*...
                vProductivity(Productivity)*(vGridLabor(Labor))^(1-aalpha);
        end
    end
end
tic 
maxDifference = 100.0;
tolerance = 0.0001;
iteration = 0;
while (maxDifference>tolerance) 
    
    expectedValueFunction = ValueFunction*Transition';
    
    if (mod(iteration,10)==0 || iteration ==1)
        
    for Productivity = 1:GridProductivity
        
        
        gridCapitalNextPeriod = 1;
        
        for Capital = 1:GridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            gridLabor=1;
            
            for CapitalNextPeriod = gridCapitalNextPeriod:GridCapital
                
                valueHighSoFar1=-1000.0;
                
                for Labor = gridLabor:GridLabor
                    
                consumption = Output(Capital,Labor,Productivity)-vGridCapital(CapitalNextPeriod)+vGridCapital(Capital)*(1-ddelta);
                valueProvisional =log(consumption)-kappa*vGridLabor(Labor)^(1+1/nu)/(1+1/nu)+...
                bbeta*expectedValueFunction(CapitalNextPeriod,Productivity);              
                
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
                  
            end
            
            ValueFunctionNew(Capital,Productivity) = valueHighSoFar;
            PolicyFunction(Capital,Productivity) = capitalChoice;
            LaborFunction(Capital,Productivity)=laborChoice;
            GridPolicyFunction(Capital,Productivity)=gridCapitalNextPeriod;
        end
    end
    
    else
             for Productivity = 1:GridProductivity
                
                 for Capital=1:GridCapital
                
                consumption = (vGridCapital(Capital))^aalpha*vProductivity(Productivity)*...
                    (LaborFunction(Capital,Productivity))^(1-aalpha)-PolicyFunction(Capital,Productivity)+...
                    vGridCapital(Capital)*(1-ddelta);
                
                gridCapitalNextPeriod1=GridPolicyFunction(Capital,Productivity);
                
                ValueFunctionNew(Capital,Productivity)=log(consumption)-...
                    kappa*LaborFunction(Capital,Productivity)^(1+1/nu)/(1+1/nu)+ ...
                    bbeta*expectedValueFunction(gridCapitalNextPeriod1,Productivity);
                 end
            end
    end  
    
    maxDifference = max(max(abs(ValueFunctionNew-ValueFunction)));
    ValueFunction = ValueFunctionNew;
    
    iteration = iteration+1;
    
           
end


figure(1)
plot(vGridCapital,ValueFunction)
xlim([vGridCapital(1) vGridCapital(GridCapital)])
xlabel('k');
ylabel('v(k)')

figure
plot(vGridCapital,PolicyFunction)
xlim([vGridCapital(1) vGridCapital(GridCapital)])
xlabel('k');
ylabel('g(k)')

y = PolicyFunction.^aalpha.*LaborFunction.^(1-aalpha);
k = PolicyFunction;
h = LaborFunction;