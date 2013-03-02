function [I,X,U] = abm(individualValues,network,sigmaFactors,populationSize,totalTime,networkAdj,initialAdopter,initialAware,initialUnaware,d,c,b,bb,k,P,sigma1,sigma2)

format long
TRUE = 1;
FALSE = 0;

numX = initialAdopter*populationSize; % number of adopter individuals
numI = initialAware*populationSize; % number of aware individuals
numU = initialUnaware*populationSize; % number of unaware individuals
X = zeros(1,totalTime);X(1) = initialAdopter; % fraction of adopters at each time step
I = zeros(1,totalTime);I(1) = initialAware; % fraction of aware at each time step
U = zeros(1,totalTime);U(1) = initialUnaware; % fraction of unaware at each time step

%avgNeighbourhoodSize = round(mean(sum(networkAdj,1)));
maxNeighbourhoodSize = max(sum(networkAdj,1));
Nz = maxNeighbourhoodSize;

% States %
unaware = 0;
aware = 1;
adopter = 2;

%%% Initializing matrix of individuals and their attributes %%%
attributes = zeros(9,populationSize);
%%row1% current state
%%row2% next turn state (i.e. state to transition to)
%%row3% c - advertising effect
%%row4% b - contact rate
%%row5% k - quickness of adoption rate
%%row6% P - personal perceived cost
%%row7% sigma1 - importance placed on being green
%%row8% sigma2 - importance placed on price vs. neighbourhood adoption
%%row9% whether an individual is an innovator (1) or not (0)
innovator = 1;

if numX > 0 % initializing current state for initial adopters
    for i = 1:numX
        attributes(1:2,i) = adopter;
        attributes(9,i) = innovator;
    end
end
if numI > 0 % initializing current state for initial aware
    for i = numX+1:numX+numI
        attributes(1:2,i) = aware;
    end
end

for i = 1:populationSize
    if individualValues == TRUE % ABMs 2 & 3
        attributes(3,i) = c(i); % setting advertising effects
        attributes(4,i) = b(i); % setting contact rates
        attributes(5,i) = k(i); % setting quickness of adoption rates
        attributes(6,i) = P(i); % setting cost
        attributes(7,i) = sigma1(i); % setting green factor
        attributes(8,i) = sigma2(i); % setting social influence factor
    else % ABM 1
        attributes(3,i) = c; % setting advertising effects
        attributes(4,i) = b; % setting contact rates
        attributes(5,i) = k; % setting quickness of adoption rates
        attributes(6,i) = P; % setting cost
        attributes(7,i) = sigma1; % setting green factor
        attributes(8,i) = sigma2; % setting social influence factor
    end
end

%%% Creating matrix to store individual's states over time for Netlogo %%%
states = zeros(totalTime+1,populationSize);
states(1,:) = attributes(1,:);


%%% Time Step %%%
for t = 1:totalTime
%%% Calculating number of aware/adopter individuals and average contact rate of aware/adopter individuals %%%
    if individualValues == TRUE && network == FALSE % ABM 2
        avgbAwareNeighbours = 0;
        avgbAdopterNeighbours = 0;
        numNeighbours = populationSize;
        numAwareNeighbours = numI;
        numAdopterNeighbours = numX;
        for i = 1:populationSize
            if attributes(1,i) == aware
                avgbAwareNeighbours = avgbAwareNeighbours + attributes(4,i);
            elseif attributes(1,i) == adopter
                avgbAdopterNeighbours = avgbAdopterNeighbours + attributes(4,i);
            end
        end
        if numAwareNeighbours ~= 0
            avgbAwareNeighbours = avgbAwareNeighbours/numAwareNeighbours;
        end
        if numAdopterNeighbours ~= 0
            avgbAdopterNeighbours = avgbAdopterNeighbours/numAdopterNeighbours;
        end
    elseif individualValues == FALSE && network == FALSE % ABM 1
        avgbAwareNeighbours = b;
        avgbAdopterNeighbours = bb;
        numNeighbours = populationSize;
        numAwareNeighbours = numI;
        numAdopterNeighbours = numX;
    end
%%% Action %%%
    for agent = 1:populationSize
        %%% Unaware %%%
        if attributes(1,agent) == unaware
            unawareActions()
        %%% Aware %%%
        elseif attributes(1,agent) == aware
            awareActions()
        %%% Adopter %%%
        elseif attributes(1,agent) == adopter
            adopterActions()
        end
    end
    attributes(1,:) = attributes(2,:); % transition to next turn
    states(t+1,:) = attributes(2,:); % storing states at time step t
    X(t+1) = numX/populationSize; % calculating fraction of adopters at end of time step t
    I(t+1) = numI/populationSize; % calculating fraction of aware at end of time step t
    U(t+1) = numU/populationSize; % calculating fraction of unaware at end of time step t
end


%%% Subfunctions %%%

%--------------------------------------------------------------------------

function unawareActions()

format long
%%% Calculating number of aware/adopter neighbours and average contact rate of aware/adopter neighbours %%%
if individualValues == TRUE && network == TRUE % ABM 3
    avgbAwareNeighbours = 0;
    avgbAdopterNeighbours = 0;
    numNeighbours = 0;
    numAwareNeighbours = 0;
    numAdopterNeighbours = 0;
    for neighbour = 1:populationSize
        if networkAdj(agent,neighbour) == 1
            numNeighbours = numNeighbours + 1;
            if attributes(1,neighbour) == aware
                avgbAwareNeighbours = avgbAwareNeighbours + attributes(4,neighbour);
                numAwareNeighbours = numAwareNeighbours + 1;
            elseif attributes(1,neighbour) == adopter
                avgbAdopterNeighbours = avgbAdopterNeighbours + attributes(4,neighbour);
                numAdopterNeighbours = numAdopterNeighbours + 1;
            end
        end
    end
    if numAwareNeighbours ~= 0
        avgbAwareNeighbours = avgbAwareNeighbours/numAwareNeighbours;
    end
    if numAdopterNeighbours ~= 0
        avgbAdopterNeighbours = avgbAdopterNeighbours/numAdopterNeighbours;
    end
end

A = (-6*((t/totalTime)-0.3)^2+0.6).*((t/totalTime)<=0.3) +...
    ((1/(sqrt(2*pi)))*(exp(-((20*((t/totalTime)-0.3)^2)/2))) -...
    (1/(sqrt(2*pi)))+0.6).*((t/totalTime)>0.3); % advertising function
%%% Transition %%%
if rand(1) < attributes(3,agent)*A +...
        avgbAwareNeighbours*(numAwareNeighbours/numNeighbours) +...
        avgbAdopterNeighbours*(numAdopterNeighbours/numNeighbours) % likelihood of becoming aware [from the DE model]
    attributes(2,agent) = aware;
    numU = numU - 1;
    numI = numI + 1;
end

end

%--------------------------------------------------------------------------

function awareActions()

format long
%%% Calculating number of adopter neighbours %%%
numAdopterNeighbours = 0;
if sigmaFactors == TRUE % ABM 4
    for neighbour = 1:populationSize
        if networkAdj(agent,neighbour) == 1
            if attributes(1,neighbour) == adopter
                numAdopterNeighbours = numAdopterNeighbours + 1;
            end
        end
    end
end
%%% Transition %%%
if rand(1) < exp(-d*attributes(6,agent))*(1-attributes(8,agent))+(attributes(8,agent)*(numAdopterNeighbours/Nz)) % likelihood of becoming a potential adopter [from the DE model]
    if rand(1) < attributes(5,agent) % likelihood of adopting after being a potential adopter [from the DE model]
        attributes(2,agent) = adopter;
        numI = numI - 1;
        numX = numX + 1;
    end
end

end

%--------------------------------------------------------------------------

function adopterActions()

format long
%%% Transition %%%
if rand(1) > exp(-d*attributes(6,agent)) && attributes(9,agent) ~= innovator % likelihood of finding the price unacceptable
    if rand(1) < attributes(5,agent)     % and disadopting
        attributes(2,agent) = aware;
        numX = numX - 1;
        numI = numI + 1;
    end
end

end

%--------------------------------------------------------------------------

end