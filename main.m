format long
%%% Variables %%%
populationSize = 1000;
totalTime = 300;
networkAdj = importdata('adj_pa_1000.txt'); % import adjacency matrix from netlogo output
initialAdopter = 0.01;
initialAware = 0.02;
initialUnaware = 0.97;


%%% ABM Toggles %%%
TRUE = 1;
FALSE = 0;
individualValues = TRUE;
network = TRUE;
sigmaFactors = TRUE;


%%% Parameters %%%
% ABM 2 with individual values
if individualValues == TRUE
    d = 0.00039; % price sensitivity
    c = betarnd(2,2,1,populationSize)/10; % advertising effect
    b = 0.0000316*betarnd(2,2,1,populationSize);bb = 0; % contact rate
    k = 0.11492*betarnd(2,2,1,populationSize); % how fast potential adopters move to adopt
    P = 7000*betarnd(2,2,1,populationSize); % personal price
    % ABM 4 with sigma factors
    if sigmaFactors == TRUE
        sigma1 = zeros(1,populationSize); % green factor
        unif = rand(1,populationSize);
        for i = 1:populationSize
            sigma1(i) = (-log(1-unif(i))/(0.4*exp(1/exp(0.4*unif(i)))))/9;
        end
        sigma2 = betarnd(2,2,1,populationSize); % social influence factor
    else
        sigma1 = zeros(1,populationSize);
        sigma2 = zeros(1,populationSize);
    end
else
    % ABM 1 based on DE
    d = 0.00039; % price sensitivity
    c = 0.05; % advertising effect
    b = 0.0000158; % aware contact rate
    bb = 0.0000158; % adopter contact rate
    k = 0.05746; % how fast potential adopters move to adopt
    P = 3500; % personal price
    sigma1 = 0;
    sigma2 = 0;
end
% ABM 3 with network
if network == FALSE
    networkAdj = 1;
end


%%% Run ABM %%%
[I,X,U] = abm(individualValues,network,sigmaFactors,populationSize,totalTime,networkAdj,initialAdopter,initialAware,initialUnaware,d,c,b,bb,k,P,sigma1,sigma2);

eq = [I(totalTime+1);X(totalTime+1);U(totalTime+1)];
nonCumulI = zeros(1,length(I));
nonCumulX = zeros(1,length(X));
nonCumulU = zeros(1,length(U));
for i = 2:length(I)
    nonCumulI(i) = I(i) - I(i-1);
    nonCumulX(i) = X(i) - X(i-1);
    nonCumulU(i) = U(i) - U(i-1);
end
time = [1:1:totalTime+1]';


%%% Plot %%%
hold on
box on
set(gca,'FontSize',16)
plot(time,I,'Color',[0,51/255,153/255],'LineStyle','--','LineWidth',4);
plot(time,X,'Color',[0,0,0],'LineStyle','-','LineWidth',4);
plot(time,U,'Color',[222/255,125/255,0],'LineStyle','-.','LineWidth',4);
%set(legend('Aware [I(t)]','Adopters [X(t)]','Unaware [U(t)]'),'Orientation','horizontal','Position',[0.508155067817327 0.0154382077564343 0.39626727046901 0.0288832363264762]);
title('Agent-based model')
ylabel('fraction of the population')
xlabel('t')
xlim([0 totalTime+1])
ylim([0 1])