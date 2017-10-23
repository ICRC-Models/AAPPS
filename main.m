%%
close all;clear all;clc
load('genParams')
load('gcParams')
load('partnerParams')
stepsPerYear = 12;
tspan = 1 / stepsPerYear;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1) =  1000; % N
popInitial(1 , 1 , 2) = 100469 - 5000; % N , urethral
popInitial(1 , 1 , 3) = 4000; % N , pharyngeal
popInitial(2 , 1 , 1) = 424 * 0.9; % I
popInitial(2 , 2 , 2 : 3) = 424 * 0.1 / 3;
popInitial(3 , 1 , 1) = 390; % K
popInitial(4 , 1 , 1) = 3972; % V
popInitial(5 , 1 , 1) = 10526 + 27366; % P
popInitial(: , 2 , 2 : 3) = popInitial(: , 1 , 2 : 3) + 50;
%%
partners = c;
pop = popInitial;
sympref('HeavisideAtOrigin' , 1);
%%
dPop = ode45(gcHivTransMat + mixInfect , tspan , popInitial);
