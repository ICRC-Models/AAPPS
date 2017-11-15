%%
close all;clear all;clc
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
stepsPerYear = 200;
tspan = 1 / stepsPerYear;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 3) =  100469 - 5000; % N (low risk)
popInitial(1 , 1 , 1 , 2) = 1000; % N (medium risk)
popInitial(1 , 1 , 1 , 1) = 4000; % N (high risk)
popInitial(2 , 1 , 1 , 3) = 424 * 0.9; % I (low risk)
popInitial(2 , 2 , 2 : 4 , 1 : 2) = 424 * 0.1 / 6; % I (high and medium risk)
popInitial(3 , 1 , 1 , 3) = 390; % K
popInitial(4 , 1 , 1 , 3) = 3972; % V
popInitial(5 , 1 , 1 , 3) = 10526 + 27366; % P
popInitial(: , 1 , 1 , :) = popInitial(: , 1 , 1 , :) + 1;
%%
partners = c;
pop = popInitial;
sympref('HeavisideAtOrigin' , 1);
%%
dPop = ode45(gcHivTransMat + mixInfect , tspan , popInitial);
