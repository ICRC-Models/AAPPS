%%
close all;clear all; clc
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
stepsPerYear = 50;
startYear = 1985;
endYear = 2020;
tspan = startYear : 1 / stepsPerYear : endYear;
tVec = tspan;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 3) =  100469 * 0.99;% - 10000 + 10526 + 27366; % N (low risk)
popInitial(1 , 1 , 1 , 2) = 9000 * 0.85; 
popInitial(1 , 1 , 1 , 1) = 1000 * 0.8;
popInitial(1 , 2 , 2 : 4 , 1 : 2) = 10;
popInitial(2 , 1 , 1 , 1) = 1000 * 0.2;
popInitial(2 , 1 , 1 , 2) = 9000 * 0.15;
popInitial(2 , 1 , 1 , 3) = 0.01 * 100460;
condUse = [0.44 , 0.25 , 0.23]; % condom usage by risk group (high, med, low)
riskVec = zeros(risk , 1);
for r = 1 : risk
    riskVec(r) = sum(sum(sum(popInitial(: , : , : , r)))) ./ sum(popInitial(:)); % risk structure remains constant thorughout simulation
end
%%
partners = c;
% sympref('HeavisideAtOrigin' , 1);
%% Scale up vectors

intStart = startYear; % start year for intervention
intPlat = startYear + 5; % plateau year for intervention

% partner services
psTreatMatTarget = 2 .* (1 - exp(-(p_ps .* kDiagTreat)));
% routine treatment scale-up
routineTreatMatTarget = 2 .* (1 - exp(-(p_routine .* kDiagTreat)));
routineTreatMat_init = 0.1 .* routineTreatMatTarget;

% intial value before intervention starts
intStartInd = round((intStart - startYear) * stepsPerYear) + 1; % index corresponding to intervention start year
intPlatInd = round((intPlat - startYear) * stepsPerYear) + 1; % index corresponding to intervention plateau year

hivScreenStart = 1995;
hivScreenPlat = 2000;

intStartInd_HivScreen = round((hivScreenStart - startYear) * stepsPerYear) + 1; % index corresponding to HIV screen start year
intPlatInd_HivScreen = round((hivScreenPlat - startYear) * stepsPerYear) + 1; % index corresponding to HIV screen plateau year

% ramp up from intervention start year to plateau year

% increment in GC and HIV screening from start year to plateau year
kHivScreen_init = 0;
kHivScreen = 0.5; % 50% HIV screen rate plateau (assumption) TEST 1/2
d_psTreatMat = psTreatMatTarget ./ (intPlatInd - intStartInd); % increment in GC and HIV screening through PS from start to plateau year 
d_routineTreatMat = (routineTreatMatTarget - routineTreatMat_init) ./ (intPlatInd - intStartInd); % increment in GC and HIV screening through routine screening from start to plateau year 
d_kHivScreen = (kHivScreen - kHivScreen_init) ./ (intPlatInd_HivScreen - intStartInd_HivScreen); % increment in HIV screening from start to plateau year 

% Scale factors for PS and routine screening
fScale = zeros(length(tspan) , 1);
fScale(intPlatInd : end) = intPlatInd - intStartInd; % scale factor between intervention start and plateau years
fScale(intStartInd : intPlatInd) = [0 : intPlatInd - intStartInd]; % scale factor for plateau year onward 

% Scale factors for HIV screening
fScale_HivScreen = zeros(length(tspan) , 1);
fScale_HivScreen(intPlatInd_HivScreen : end) = intPlatInd_HivScreen - intStartInd_HivScreen; % scale factor for plateau value
fScale_HivScreen(intStartInd_HivScreen : intPlatInd_HivScreen) = [0 : intPlatInd_HivScreen - intStartInd_HivScreen];

%%
% scale-up HIV serosorting from 1990 (hAssStart) to 2000 (hAssPlat)
hAssStart = startYear; % HIV assorting start year
hAssPlat = 2000; % HIV assorting plateau year
hAssStartInd = round((hAssStart - startYear) * stepsPerYear) + 1; % index corresponding to HIV assorting start year
hAssPlatInd = round((hAssPlat - startYear) * stepsPerYear) + 1; % index corresponding to HIV assorting plateau year
hAssortTarget = 0.8; % Target plateau value for HIV assortativity
hAssort_init = 0; % Initial HIV assortativity value
hScale = zeros(length(tspan) , 1);
d_hAssort = (hAssortTarget - hAssort_init) ./ (hAssPlatInd - hAssStartInd);
hScale(hAssStartInd : end) = hAssPlatInd - hAssStartInd;
hScale(hAssStartInd : hAssPlatInd) = [0 : hAssPlatInd - hAssStartInd];

rAssort = 0.5; % risk assortativity

kHivTreat = 1-exp(-0.9); % Hypothetical HIV treatment rate
% gcClear = gcClear;
%%
figure()
subplot(2 , 1 , 1)
plot(tVec , hScale .* d_hAssort)
title('HIV Assortativity')

% subplot(2 , 2 , 2)
% plot(tVec ,  fScale .* d_psTreatMat)
% title('Partner Services Scale-Up')
% 
% subplot(2 , 2 , 3)
% plot(tVec ,  fScale .* d_routineTreatMat)
% title('Routine Screen Scale-Up')

subplot(2 , 1 , 2)
plot(tVec , fScale_HivScreen .* d_kHivScreen)
title('HIV Screen Scale-Up')
%% ODE solver
% ODE options: Absolute tolerance: 10^-3 (ODE solver ignores differences in values less than 10^-3). 
% NonNegative option ensures solutions with negative values are not generated. 
options = odeset('AbsTol', 1e-4, 'NonNegative' , 1);
disp('Running...')
[t , pop] = ode45(@(t , pop) mixInfect(t , pop , hivStatus , stiTypes , sites , ...
    risk , kBorn , kDie , gcClear , d_routineTreatMat , routineTreatMat_init , ...
    p_symp , fScale ,fScale_HivScreen , d_psTreatMat , kDiagTreat , ...
    kHivScreen_init , d_kHivScreen , kHivTreat , partners , acts , riskVec ,...
    condUse , d_hAssort , hScale , hAssort_init , rAssort , tVec) , ...
    tspan , popInitial , options);

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
plot(t , totalPop)
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
plot(t , gc_hivInf , t , gc_hivTreated , t , gc_hivTested , t , gc_prepImm)
title('GC-HIV CoInfection')
legend('GC + HIV Infected' , 'GC + HIV Treated' , 'GC + HIV Tested' ,...
    'GC + PrEP/Immune')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 100])

%% GC and HIV susceptibles
figure()
gc_Sus = sum(sum(sum(sum(pop(: , 1 : hivStatus , 1 , 1 : 4 , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
subplot(2 , 1 , 1)
plot(t , gc_Sus)
title('GC-Susceptible')
xlabel('Year'); ylabel('Proportion (%)')
axis([tVec(1) tVec(end) 0 100])

hiv_sus = sum(sum(sum(sum(pop(: , 1 , 1 : stiTypes , 1 : sites , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
subplot(2 , 1 , 2)
plot(t , hiv_sus)
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
