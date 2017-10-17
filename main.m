%%
load('genParams')
load('gcParams')
load('partnerParams')
stepsPerYear = 12;
tspan = 1 / stepsPerYear;
popInitial = zeros(hivStatus , stiTypes , sites);
popInitial(1 , 1 , 1) = 100469; % N
popInitial(2 , 1 , 1) = 424 * 0.9; % I
popInitial(2 , 2 , 1 : 3) = 424 * 0.1 / 3;
popInitial(3 , 1 , 1) = 390; % K
popInitial(4 , 1 , 1) = 3972; % V
popInitial(5 , 1 , 1) = 10526 + 27366; % P
popInitial(: , 2 , 2 : 3) = popInitial(: , 1 , 2 : 3) + 50;
%%
dPop = ode45(gcHivTransMat + mixInfect , tspan , popInitial);
