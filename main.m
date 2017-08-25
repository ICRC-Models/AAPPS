
stepsPerYear = 12;
tspan = 1 / stepsPerYear;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 1) = 100469; % N
popInitial(2 , 1 , 1 , 1) = 424; % I
popInitial(3 , 1 , 1 , 1) = 390; % K
popInitial(4 , 1 , 1 , 1) = 3972; % V
popInitial(5 , 1 , 1 , 1) = 10526; % P
popInitial(5 , 1 , 1 , 2) = 27366; % P high risk
dPop = ode45(gcHivTransMat + mixInfect , tspan , popInitial);
