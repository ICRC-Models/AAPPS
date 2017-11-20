%%
close all;clear all; clc
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
stepsPerYear = 200;
tspan = 2017 : 1 / stepsPerYear : 2050;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 3) =  100469 - 5000; % N (low risk)
popInitial(1 , 1 , 1 , 2) = 1000; % N (medium risk)
popInitial(1 , 1 , 1 , 1) = 4000; % N (high risk)
popInitial(2 , 1 , 1 , 3) = 424 * 0.5; % I (low risk)
popInitial(2 , 2 , 2 : 4 , 1 : 2) = 424 * 0.5 / 6; % I (high and medium risk)
popInitial(3 , 1 , 1 , 3) = 390; % K
popInitial(4 , 1 , 1 , 3) = 3972; % V
popInitial(5 , 1 , 1 , 3) = 10526 + 27366; % P
popInitial(: , 1 , 1 , :) = popInitial(: , 1 , 1 , :) + 1;
%%
partners = c;
sympref('HeavisideAtOrigin' , 1);
%%
disp('Running...')
[t , pop] = ode45(@(t , pop) mixInfect(t , pop , hivStatus , stiTypes , sites , ...
    risk , kBorn , kDie , pSite , gcTreat , k_gcPS , gcClear , kInt , ...
    perActInf , rTreat_v , k_toPrep, k_prepOut , kTest_i , ...
    kTreat_i , k_treatOut , kprep_u , kprep_u_ps , kTest_u , kTest_u_ps , ...
    kTreat_u , kTreat_u_ps , kTreat_k , kTreat_k_ps , kprep_p , kprep_p_ps , kTest_p , ...
    kTest_p_ps , kTreat_p , kTreat_p_ps , kprep_r , kprep_r_ps , kTest_r , ...
    kTest_r_ps , kTreat_r , kTreat_r_ps , partners , p_cond , acts) , tspan , popInitial);

disp('Finished solving')

% reshape pop vector
pop = reshape(pop , [size(pop , 1) , hivStatus , stiTypes , sites , risk]);




%% plot settings (completely optional)
colors = [241, 90, 90;
          240, 196, 25;
          78, 186, 111;
          45, 149, 191;
          149, 91, 165]/255;

set(groot, 'DefaultAxesColor', [10, 10, 10]/255);
set(groot, 'DefaultFigureColor', [10, 10, 10]/255);
set(groot, 'DefaultFigureInvertHardcopy', 'off');
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(groot, 'DefaultAxesColorOrder', colors);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultTextColor', [1, 1, 1]);
set(groot, 'DefaultAxesXColor', [1, 1, 1]);
set(groot, 'DefaultAxesYColor', [1, 1, 1]);
set(groot , 'DefaultAxesZColor' , [1 , 1 ,1]);
set(0,'defaultAxesFontSize',14)
ax = gca;
ax.XGrid = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.GridColor = [1, 1, 1];
ax.GridAlpha = 0.4;

%% plots

totalPop = sum(sum(sum(sum(pop , 2) , 3) , 4) , 5);
plot(t , totalPop)
title('Population Size'); xlabel('Year'); ylabel('Persons')

allGC = sum(sum(sum(pop(: , 1 : hivStatus , 2  , 1 : sites , 1 : 3) , 2) , 4) , 5);
figure()
plot(t , allGC ./ totalPop * 100)
title('GC Positive Individuals')
xlabel('Year'); ylabel('Prevalence (%)')

hivInf = sum(sum(sum(pop(: , 2 , : , : , :) , 3) , 4) , 5);
figure()
plot(t , hivInf ./ totalPop * 100)
title('HIV Infectious')
xlabel('Year'); ylabel('Prevalence (%)')



