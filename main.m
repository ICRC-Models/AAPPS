%%
close all;clear all; clc
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
stepsPerYear = 200;
startYear = 1990;
endYear = 2020;
tspan = startYear : 1 / stepsPerYear : endYear;
tVec = tspan;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 3) =  100469;% - 10000 + 10526 + 27366; % N (low risk)
popInitial(1 , 1 , 1 , 2) = 9000;
popInitial(1 , 1 , 1 , 1) = 1000;
% popInitial(1 , 2 , 2 : 4 , 1) = 10;
% popInitial(1 , 2 , 2 : 4 , 2) = 10;
% popInitial(1 , 2 , 2 : 4 , 3) = 10;
% popInitial(1 , 1 , 1 , 2) = 5000; % N (medium risk)
popInitial(1 , 2 , 2 , 1 : 2) = 20;
popInitial(2 , 1 , 1 , 1 : 2) = 2;
% popInitial(2 , 1 , 1 , 2) = 1; %424 * 0.9; % I (Medium risk)
% popInitial(2 , 1 , 1 , 3) = 1;
% popInitial(2 , 1 , 2 : 4 , 2) = 10  * 0.1 / 3; % I (medium risk)
% popInitial(3 , 1 , 1 , 1) = 390 * 0.1; % K
% popInitial(4 , 1 , 1 , 2) = 3972; % V
%popInitial(5 , 1 , 1 , 3) = 10526 + 27366; % P
popInitial(: , : , : , :) = popInitial(: , : , : , :) + 10^-9;
condUse = [0.44 , 0.25 , 0.23]; % condom usage by risk group (high, med, low)
riskVec = zeros(risk , 1);
for r = 1 : risk
    riskVec(r) = sum(sum(sum(popInitial(: , : , : , r)))) ./ sum(popInitial(:));
end
%%
partners = c;
sympref('HeavisideAtOrigin' , 1);
%% Scale up vectors

intStart = 2000; % start year for intervention
intPlat = 2010; % plateau year for intervention
targetVal = 0.3;% 0.05; %0.3; % target value in plateau year

% intial value before intervention starts
kTest_p_ps = zeros(size(tspan));
kTreat_p_ps = zeros(size(tspan));
kTest_u_ps = zeros(size(tspan));
kTreat_u_ps = zeros(size(tspan));
kTest_r_ps = zeros(size(tspan));
kTreat_r_ps = zeros(size(tspan));

k_gcPS_target = k_gcPS(1) / 10;
k_gcPS = zeros(size(tspan));


intStartInd = round((intStart - startYear) * stepsPerYear);
intPlatInd = round((intPlat - startYear) * stepsPerYear);
tScale = tspan(intStartInd : intPlatInd) - intStartInd;

% ramp up from intervention start year to plateau year
kTest_p_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;
kTreat_p_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;
kTest_u_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;
kTreat_u_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;
kTest_r_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;
kTreat_r_ps(intStartInd : intPlatInd) = targetVal / (intPlat - intStart) ... 
    * tScale;

k_gcPS(intStartInd : intPlatInd) = k_gcPS_target / (intPlat - intStart)...
    * tScale;


% rate following plateau year
kTest_p_ps(intPlatInd + 1 : end) = kTest_p_ps(intPlatInd);
kTreat_p_ps(intPlatInd + 1 : end) = kTreat_p_ps(intPlatInd);
kTest_u_ps(intPlatInd + 1 : end) = kTest_u_ps(intPlatInd);
kTreat_u_ps(intPlatInd + 1 : end) = kTreat_u_ps(intPlatInd);
kTest_r_ps(intPlatInd + 1 : end) = kTest_r_ps(intPlatInd);
kTreat_r_ps(intPlatInd + 1 : end) = kTreat_r_ps(intPlatInd);

k_gcPS(intPlatInd + 1 : end) = k_gcPS(intPlatInd);

%% ODE solver
% k_treatOut = 0.5 * k_treatOut; % test
disp('Running...')
[t , pop] = ode23s(@(t , pop) mixInfect(t , pop , hivStatus , stiTypes , sites , ...
    risk , kBorn , kDie , pSite , gcTreat , k_gcPS , gcClear , kInt , ...
    perActInf , rTreat_v , k_toPrep, k_prepOut , kTest_i , ...
    kTreat_i , k_treatOut , kprep_u , kprep_u_ps , kTest_u , kTest_u_ps , ...
    kTreat_u , kTreat_u_ps , kTreat_k , kTreat_k_ps , kprep_p , kprep_p_ps , kTest_p , ...
    kTest_p_ps , kTreat_p , kTreat_p_ps , kprep_r , kprep_r_ps , kTest_r , ...
    kTest_r_ps , kTreat_r , kTreat_r_ps , partners , p_cond , acts , riskVec , ...
    condUse , tVec) , tspan , popInitial);

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
figure()
totalPop = sum(sum(sum(sum(pop , 2) , 3) , 4) , 5);
plot(t , round(totalPop))
title('Population Size'); xlabel('Year'); ylabel('Persons')

allGC = sum(sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 4) , 5);
figure()
plot(t , allGC ./ totalPop * 100)
title('GC Positive Individuals')
xlabel('Year'); ylabel('Prevalence (%)')
axis([tVec(1) tVec(end) 0 100])


allGC_site = squeeze(bsxfun(@rdivide , ...
    sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 5) , ...
    totalPop)) * 100;
figure()
area(t , allGC_site)
title('GC Prevalence Site')
xlabel('Year'); ylabel('Prevalence (%)')
legend('Rectal' , 'Urethral' , 'Pharyngeal')
axis([tVec(1) tVec(end) 0 100])

% HIV
hivInf = sum(sum(sum(pop(: , 2 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;
hivTreated = sum(sum(sum(pop(: , 4 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;
hivTested = sum(sum(sum(pop(: , 3 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;
prepImm = sum(sum(sum(pop(: , 5 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;

figure()
plot(t , hivInf , t , hivTreated , t , hivTested , t , prepImm)
title('HIV Prevalence')
legend('Infected' , 'Treated' , 'Tested' , 'PrEP/Immune')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 20])

figure()
area(t , [hivInf , hivTreated , hivTested , prepImm])
title('HIV Prevalence')
legend('Infected' , 'Treated' , 'Tested' , 'PrEP/Immune')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 20])

% GC-HIV
gc_hivInf = sum(sum(pop(: , 2 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_hivTreated = sum(sum(pop(: , 4 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_hivTested = sum(sum(pop(: , 3 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_prepImm = sum(sum(pop(: , 5 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;

figure()
plot(t , round(gc_hivInf , 2) , t , round(gc_hivTreated , 2) , t , round(gc_hivTested , 2) , t , round(gc_prepImm , 2))
title('GC-HIV CoInfection')
legend('GC + HIV Infected' , 'GC + HIV Treated' , 'GC + HIV Tested' ,...
    'GC + PrEP/Immune')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 100])

%% GC and HIV susceptibles
figure()
gc_Sus = sum(sum(sum(sum(pop(: , 1 : hivStatus , 1 , 1 : 4 , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
subplot(2 , 1 , 1)
plot(t , round(gc_Sus , 2))
title('GC-Susceptible')
xlabel('Year'); ylabel('Proportion (%)')
axis([tVec(1) tVec(end) 0 100])

hiv_sus = sum(sum(sum(sum(pop(: , 1 , 1 : stiTypes , 1 : sites , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
subplot(2 , 1 , 2)
plot(t , round(hiv_sus , 2))
title('HIV-Susceptible')
xlabel('Year'); ylabel('Proportion (%)')
axis([tVec(1) tVec(end) 0 100])

%%
% HIV
totalPop_Risk = squeeze(sum(sum(sum(pop(: , : , : , : , :), 2), 3) , 4));

hivInf = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 2 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
hivTreated = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 4 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
hivTested = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 3 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
prepImm = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 5 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;

hivArray = {hivInf, hivTested, hivTreated, prepImm};
hivPlotTitle = {'Infected' , 'Tested' , 'Treated' , 'Immune'};

for i = 1 : size(hivArray,2)
    figure()
    plot(t , hivArray{i})
    title(['HIV ' , hivPlotTitle{i}])
    legend('High Risk' , 'Medium Risk' , 'Low Risk')
    xlabel('Year'); ylabel('Prevalence (%)')
end
