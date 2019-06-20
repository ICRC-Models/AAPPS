%%
close all;clear all; clc

%%
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
stepsPerYear = 50;
startYear = 1980;
endYear = 2041;
tspan = startYear : 1 / stepsPerYear : endYear;
tVec = tspan;
popInitial = zeros(hivStatus , stiTypes , sites , risk);
popInitial(1 , 1 , 1 , 3) = 14100 * 1.8 - 300;% - 10000 + 10526 + 27366; % N (low risk)
popInitial(1 , 1 , 1 , 2) = 1800 * 1.8 - 300;
popInitial(1 , 1 , 1 , 1) = 1800 * 1.8 - 300;
popInitial(1 , 2 , 2 : 4 , 1) = 190 * 3;
popInitial(1 , 2 , 2 : 4 , 2) = 180 * 3;
popInitial(1 , 2 , 2 : 4 , 3) = 460 * 3;
% popInitial(2 , 2 , 2 : 4 , 1) = 500;
% popInitial(2 , 2 , 2 : 4 , 2) = 600;
popInitial(2 , 1 , 1 , 1) = 360 * 1.8;
popInitial(2 , 2 , 2:4 , 1) = 40 * 1.8;
popInitial(2 , 1 , 1 , 2) = 315 * 1.8;
popInitial(2 , 2 , 2:4 , 2) = 35 * 1.8;
% condUse = [0.23 , 0.25 , 0.44]; % condom usage by risk group (high, med, low)
% condUse = [0.23, 0.29, 0.4]; % condom usage by risk group (high, med, low) TEST!!!
riskVec = zeros(risk , 1);
popNew = popInitial .* 0.01;
popVec = zeros(length(tspan), hivStatus , stiTypes , sites , risk);

%%
partners = c;
kDie = 0.0018;
kBorn = 7 * kDie;

% sympref('HeavisideAtOrigin' , 1);
%% Interventions

% partner services
% p_symp2 = p_symp .* [1; 2; 1];
% psTreatMatTarget = (1 - exp(-((p_ps + .7) .* p_symp )));

psTreatMatTarget = (p_ps + .3) .* p_symp .* 0;

% routine treatment scale-up
routine1TreatMat_init =  ones(3,5) .* 0.001; 
routine1TreatMatTarget = (1 - exp(-(p_routine .* .5)));

routine2TreatMat_init =  routine1TreatMatTarget; 
routine2TreatMatTarget = (1 - exp(-(p_routine .* 1.5)));

%% Scale up vectors for GC interventions and screening 

% Scale factors for PS and routine screening
intStart = 2020; % start year for intervention
intPlat = intStart + 5; % plateau year for intervention

% intial value before intervention starts
intStartInd = round((intStart - startYear) * stepsPerYear) + 1; % index corresponding to intervention start year
intPlatInd = round((intPlat - startYear) * stepsPerYear) + 1; % index corresponding to intervention plateau year
fScale = zeros(length(tspan) , 1);

d_psTreatMat = psTreatMatTarget ./ (intPlatInd - intStartInd); % increment in GC and HIV screening through PS from start to plateau year 

% ramp up routine screening for GC from 2000 year to end year
rout1Start = 2010;
rout1Plat = intStart;
rout1StartInd = round((rout1Start - startYear) * stepsPerYear) + 1; % index corresponding to intervention start year
rout1PlatInd = round((rout1Plat - startYear) * stepsPerYear) + 1;

d1_routineTreatMat = (routine1TreatMatTarget - routine1TreatMat_init) ./ (rout1PlatInd - rout1StartInd); % increment in GC and HIV screening through routine screening from start to plateau year 
rout1Scale = fScale;
rout1Scale(rout1PlatInd : end) = rout1PlatInd - rout1StartInd; % scale factor for plateau year onward
rout1Scale(rout1StartInd : rout1PlatInd) = [0 : rout1PlatInd - rout1StartInd];

rout2Start = intStart;
rout2Plat = 2040;
rout2StartInd = round((rout2Start - startYear) * stepsPerYear) + 1; % index corresponding to intervention start year
rout2PlatInd = round((rout2Plat - startYear) * stepsPerYear) + 1;

d2_routineTreatMat = (routine2TreatMatTarget - routine2TreatMat_init) ./ (rout2PlatInd - rout2StartInd); % increment in GC and HIV screening through routine screening from start to plateau year 2
rout2Scale = fScale;
rout2Scale(rout2PlatInd : end) = rout2PlatInd - rout2StartInd; % scale factor for plateau year onward
rout2Scale(rout2StartInd : rout2PlatInd) = [0 : rout2PlatInd - rout2StartInd];

%% Scale factor for HIV screening and treatment 
hivScreenStart = startYear + 5;
hivScreenPlat = 2020;

intStartInd_HivScreen = round((hivScreenStart - (startYear)) * stepsPerYear) + 1; % index corresponding to HIV screen start year
intPlatInd_HivScreen = round((hivScreenPlat - (startYear)) * stepsPerYear) + 1; % index corresponding to HIV screen plateau year

fScale(intPlatInd : end) = intPlatInd - intStartInd; % scale factor for plateau year onward
fScale(intStartInd : intPlatInd) = [0 : intPlatInd - intStartInd]; % scale factor between intervention start and plateau years  

% increment in GC and HIV screening from start year to plateau year
kHivScreen_init = 0.01;
kHivScreen = 0.8; %5; % 50% HIV screen rate plateau (assumption) TEST 1/2
d_kHivScreen = (kHivScreen - kHivScreen_init) ./ (intPlatInd_HivScreen - intStartInd_HivScreen); % increment in HIV screening from start to plateau year 

% Scale factors for HIV screening
fScale_HivScreen = zeros(length(tspan) , 1);
fScale_HivScreen(intPlatInd_HivScreen : end) = intPlatInd_HivScreen - intStartInd_HivScreen; % scale factor for plateau value
fScale_HivScreen(intStartInd_HivScreen : intPlatInd_HivScreen) = [0 : intPlatInd_HivScreen - intStartInd_HivScreen];

% HIV death rate 
muHiv_init = 1 - exp(-0.03);
muHiv_plat = 1 - exp(-0.001);
d_muHiv = (muHiv_plat - muHiv_init) ./ (intPlatInd_HivScreen - intStartInd_HivScreen);

% HIV treatment 
hTreatStart = startYear + 5;
hTreatPlat = 2020;
hTreatStartInd = round((hTreatStart - startYear) * stepsPerYear) + 1 ;
hTreatPlatInd = round((hTreatPlat - startYear) * stepsPerYear) + 1; 
hTreatTarget = 0.8; %1-exp(-0.45); % Hypothetical HIV treatment rate
hTreat_init = 0.001;
d_hTreat = (hTreatTarget - hTreat_init) ./ (hTreatPlatInd - hTreatStartInd);
hTscale = zeros(length(tspan) , 1);
hTscale(hTreatStartInd : end) = hTreatPlatInd - hTreatStartInd;
hTscale(hTreatStartInd : hTreatPlatInd) = [0 : hTreatPlatInd - hTreatStartInd];

%kHivTreat = 1-exp(-0.08); % Hypothetical HIV treatment rate

%%
% scale-up HIV serosorting (source: Khosropour 2016)
hAssStart = startYear ; % HIV assorting start year
hAssPlat = 2010; % HIV assorting plateau year
hAssStartInd = round((hAssStart - startYear) * stepsPerYear) + 1; % index corresponding to HIV assorting start year
hAssPlatInd = round((hAssPlat - startYear) * stepsPerYear) + 1; % index corresponding to HIV assorting plateau year
hAssortTarget = 0.5; % Target plateau value for HIV assortativity
hAssort_init = 0.2; % Initial HIV assortativity value
hScale = zeros(length(tspan) , 1);
d_hAssort = (hAssortTarget - hAssort_init) ./ (hAssPlatInd - hAssStartInd);
hScale(hAssStartInd : end) = hAssPlatInd - hAssStartInd;
hScale(hAssStartInd : hAssPlatInd) = [0 : hAssPlatInd - hAssStartInd];

% original rAssort = 0.5; % risk assortativity
%scale up risk assortativity
rAssStart = startYear ; % Risk assorting start year
rAssPlat = 2020; % HIV assorting plateau year
rAssStartInd = round((rAssStart - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting start year
rAssPlatInd = round((rAssPlat - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting plateau year
rAssortTarget = 0.5; % Target plateau value for risk assortativity
rAssort_init = 0.8; % Initial risk assortativity value
rScale = zeros(length(tspan) , 1);
d_rAssort = -(rAssort_init - rAssortTarget) ./ (rAssPlatInd - rAssStartInd);
rScale(rAssStartInd : end) = rAssPlatInd - rAssStartInd;
rScale(rAssStartInd : rAssPlatInd) = [0 : rAssPlatInd - rAssStartInd];

%scale down condom use
cAssStart = startYear ; % Risk assorting start year
cAssPlat = 2010; % HIV assorting plateau year
cAssStartInd = round((cAssStart - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting start year
cAssPlatInd = round((cAssPlat - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting plateau year
cAssortTarget = [.6, .5, .5]; % Target plateau values for condom use
cAssort_init = [.2, .25, .3]; % Initial condom use values
cScale = zeros(length(tspan) , 3);
d_cAssort = (cAssortTarget - cAssort_init) ./ (cAssPlatInd - cAssStartInd);
cScale(:, :) = (cAssPlatInd - cAssStartInd);
cScale(cAssStartInd:cAssPlatInd, : ) = [0 : cAssPlatInd - cAssStartInd; 0 : cAssPlatInd - cAssStartInd; 0 : cAssPlatInd - cAssStartInd]';

cotestStart = 1990 ; % cotesting start year 
cotestPlat = 2020; % cotesting plateau year
cotestStartInd = round((cotestStart - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting start year
cotestPlatInd = round((cotestPlat - startYear) * stepsPerYear) + 1; % index corresponding to risk assorting plateau year
cotestTarget = 0.6;
cotestInit = 0.01; % Initial risk assortativity value
cotestScale = zeros(length(tspan) , 1);
d_cotest = (cotestTarget - cotestInit) ./ (cotestPlatInd - cotestStartInd);
cotestScale(cotestPlatInd : end) = cotestPlatInd - cotestStartInd;
cotestScale(cotestStartInd : cotestPlatInd) = [0 : cotestPlatInd - cotestStartInd];

% gcClear = gcClear;

% years = endYear - startYear;
% s = 1 : (1 / stepsPerYear) : years + 1;
% newHiv = zeros(length(s) - 1, stiTypes, sites, risk);

%%
figure()
subplot(2 , 1 , 1)
plot(tVec , hScale .* d_hAssort + hAssort_init)
xlim([startYear endYear])
title('HIV Assortativity')

% subplot(2 , 2 , 2)
% plot(tVec ,  fScale .* d_psTreatMat)
% title('Partner Services Scale-Up')
% % 
% d_routineTreatMat1 = d_routineTreatMat(:, 1) + routineTreatMat_init(:, 1);
% d_routineTreatMat2 = d_routineTreatMat(:, 2) + routineTreatMat_init(:, 2);
% d_routineTreatMat3 = d_routineTreatMat(:, 3) + routineTreatMat_init(:, 3);
% d_routineTreatMat4 = d_routineTreatMat(:, 4) + routineTreatMat_init(:, 4);
% 
% subplot(3 , 2 , 3)
% plot(tVec ,  fScale .* d_routineTreatMat(:, 1))
% title('Routine Screen Scale-Up')

subplot(2 , 1 , 2)
plot(tVec , fScale_HivScreen .* d_kHivScreen + kHivScreen_init)
xlim([startYear endYear])
title('HIV Screen Scale-Up')
% 

%% ODE solver
newHiv = zeros(length(tspan)-1,stiTypes,sites,risk);
newSti = zeros(length(tspan)-1,hivStatus,sites,risk);
newGcPs = zeros(length(tspan)-1, hivStatus, stiTypes, sites, risk);

popVec(1,:,:,:,:) = popInitial;
popIn = popInitial;

disp('Running...')
tic

for time = 1:length(tspan)-1
    tspanStep = [tspan(time) , tspan(time + 1)]; % evaluate diff eqs over one time interval
    
    [t , pop, newHiv(time,:,:,:), newSti(time,:,:,:), newGcPs(time, :, :, :, :)] = ode4xtra(@(t , pop) mixInfect_Hiv_rout(t , pop , hivStatus , stiTypes , sites ,  ...
    risk , kDie , kBorn , gcClear , rout1Scale, d1_routineTreatMat , routine1TreatMat_init , ...
    rout2Scale, d2_routineTreatMat , routine2TreatMat_init, ...
    p_symp , fScale ,fScale_HivScreen , d_psTreatMat , kDiagTreat , ...
    kHivScreen_init , d_kHivScreen , d_muHiv, muHiv_init,  hTscale, d_hTreat, hTreat_init, partners , acts , riskVec ,...
    cScale, d_cAssort, cAssort_init , d_hAssort , hScale , hAssort_init, d_cotest, cotestInit, cotestStartInd, ...
    d_rAssort , rScale, rAssort_init, tVec) , tspanStep , popIn);

    popIn = reshape(pop(end,:) , [hivStatus , stiTypes , sites , risk]);
    popVec(time+1,:,:,:,:) = popIn;

end
    
disp('Finished solving')
toc
disp(' ')

% rename to avoid changes to plotting code
pop = popVec;
t = tspan';

%% output additional mixInfect internal variables (ie. incidence)
% 
% disp('Getting internal vars...')
% tic
% popCheck = pop;
% newHiv = zeros(size(pop , 1)-1,stiTypes,sites,risk);
% newSti = zeros(size(pop , 1)-1,stiTypes,sites,risk);
% for time = 1 : length(t) - 1 
%     popStep(:,:,:,:) = pop(time,:,:,:,:);
%     [dpop , newHivStep, newStiStep] = mixInfect(t(time) , popStep , hivStatus , stiTypes , sites ,  ...
%         risk , kDie , kBorn , gcClear , d_routineTreatMat , routineTreatMat_init , ...
%         p_symp , fScale ,fScale_HivScreen , d_psTreatMat , kDiagTreat , ...
%         kHivScreen_init , d_kHivScreen ,  hTscale, d_hTreat, hTreat_init, partners , acts , riskVec ,...
%         cScale, d_cAssort, cAssort_init , d_hAssort , hScale , hAssort_init, ...
%         d_rAssort , rScale, rAssort_init, tVec);
%     popCheck(time + 1,:,:,:,:) = reshape(dpop , [hivStatus , stiTypes , sites , risk]);
%     newHiv(time,:,:,:) = newHivStep;
%     newSti(time,:,:,:) = newStiStep;
% end
% disp('Finished')
% toc
% disp(' ')

%size(newHiv)
%newHiv = reshape(newHiv, [size(newHiv, 1), stiTypes , sites , risk]);

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
set(groot, 'DefaultAxesZColor' , [1 , 1 ,1]);
set(0, 'defaultAxesFontSize', 14)
ax = gca;
ax.XGrid = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.GridColor = [1, 1, 1];
ax.GridAlpha = 0.4;

%% Population over time plots
figure()
totalPop = sum(sum(sum(sum(pop , 2) , 3) , 4) , 5);
%totalHiv = sum(sum(sum(newHiv, 2), 3), 4);
plot(t , totalPop)
xlim([startYear endYear])
title('Population Size'); xlabel('Year'); ylabel('Persons')

figure()
totalPop_hr = sum(sum(sum(sum(pop(: , : , :  , :  , 1 ), 2), 3), 4), 5) ;
totalPop_mr = sum(sum(sum(sum(pop(: , : , :  , :  , 2 ), 2), 3), 4), 5) ;
totalPop_lr = sum(sum(sum(sum(pop(: , : , :  , :  , 3 ), 2), 3), 4), 5) ;
plot(tVec , totalPop_lr, tVec, totalPop_mr, tVec, totalPop_hr)
xlim([startYear endYear])
legend({'Low risk' , 'Medium risk' , 'High risk'},'Location', 'northwest')
legend ('boxoff')
title('Population Size'); xlabel('Year'); ylabel('Persons')

% Overall pop 
% filename = 'Result_basecase_RoutPrep3.xlsx';
% sheet = 'Overall_pop';
% xlswrite (filename, [t(1: stepsPerYear * 10 :end), totalPop(1: stepsPerYear * 10: end)], sheet)
% 
% % Pop by risk group
% filename = 'Result_basecase_RoutPrep3.xlsx';
% test = {'1' , '2', '3'};
% riskPop{1} = totalPop_hr;
% riskPop{2} = totalPop_mr;
% riskPop{3} = totalPop_lr;
% for r = 1 : risk
% sheet = ['PopByRisk', test{r}];
% xlswrite (filename, [t(1: stepsPerYear * 10 :end), riskPop{r}(1: stepsPerYear * 10: end)], sheet)
% end

%% GC prevalence plots
allGC = sum(sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 4) , 5) ./ totalPop * 100;
% figure()
% plot(tVec , allGC )
% title('GC Positive Individuals')
% xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 50])
%% 
allGC_site = squeeze(bsxfun(@rdivide , ...
    sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 5) , ...
    totalPop)) * 100;
figure()
plot(t, allGC, t , allGC_site)
xlim([startYear endYear])
title('GC Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')
legend({'Overall', 'Rectal' , 'Urethral' , 'Pharyngeal'}, 'Location', 'northwest')
% axis([tVec(1) tVec(end) 0 50])

% % 
% filename = 'Result_basecase_RoutPrep3.xlsx';
%     sheet = 'GC_prev';
%     xlswrite (filename, [t(1: stepsPerYear :end),  allGC(1: stepsPerYear : end), allGC_site(1: stepsPerYear : end, :)], sheet, 'A2');
%     col_header={'Year', 'Overall', 'Rectal' , 'Urethral' , 'Pharyngeal'};
%     xlswrite(filename, col_header, sheet, 'A1') ;

%% GC counts 
allGC_num = sum(sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 4) , 5);
GCsite_num = squeeze(sum(sum(pop(: , 1 : hivStatus , 2  , 2 : sites , 1 : 3) , 2) , 5));
figure()
plot(t, allGC_num, t, GCsite_num)
title('Number of GC cases')
xlabel('Year'); ylabel('Number of cases')
legend({'Overall', 'Rectal' , 'Urethral' , 'Pharyngeal'}, 'Location', 'northwest')
xlim([startYear endYear])
% 
% filename = 'Result_basecase_RoutPrep3.xlsx';
%     sheet = 'GC_Counts';
%     xlswrite(filename, [t(1: stepsPerYear :end), allGC_num(1: stepsPerYear : end),GCsite_num(1: stepsPerYear : end, :)], sheet, 'A2');
%     col_header={'Year', 'Overall', 'Rectal' , 'Urethral' , 'Pharyngeal'};
%     xlswrite(filename, col_header, sheet, 'A1') ;

%% HIV prevalence 
hivInf = sum(sum(sum(pop(: , 2 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;
hivTreated = sum(sum(sum(pop(: , 4 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;
hivTested = sum(sum(sum(pop(: , 3 , : , : , :) , 3) , 4) , 5) ./ totalPop * 100;

%
figure()
plot(t , hivInf , t , hivTreated , t , hivTested)
xlim([startYear endYear])
title('HIV Prevalence')
legend({'Infected' , 'Treated' , 'Tested'} , 'Location', 'northwest')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 20])
%
figure()
area(t , [hivInf , hivTreated , hivTested ])
xlim([startYear endYear])
title('HIV Prevalence')
legend('Infected' , 'Treated' , 'Tested' , 'Location', 'northwest')
xlabel('Year'); ylabel('Prevalence (%)')
% axis([tVec(1) tVec(end) 0 20])

%% Overall HIV prevalence
figure()
hivAll = sum(sum(sum(sum(pop(: , 2 : 4 , : , : , :), 2) , 3) , 4) , 5) ./ totalPop * 100;
plot(t , hivAll)
xlim([startYear endYear])
title('Overall HIV Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')

% filename = 'Result_basecase_RoutPrep3.xlsx';
%     sheet = 'HIV_prev';
%      xlswrite(filename, [t(1: stepsPerYear :end),  hivAll(1: stepsPerYear : end), ...
%         hivInf(1: stepsPerYear : end),  hivTested(1: stepsPerYear : end), hivTreated(1: stepsPerYear : end)], sheet, 'A2');
%         col_header = {'Year', 'Overall', 'Infected' , 'Tested' , 'Treated'};
%         xlswrite(filename, col_header, sheet, 'A1');

%% Prep use by risk group 
prepImm = sum(sum(sum(pop(: , 5 , : , : , :) , 3) , 4) , 5) ./ (sum(sum(sum(pop(: , 1 , : , : , :) , 3) , 4) , 5) + sum(sum(sum(pop(: , 5 , : , : , :) , 3) , 4) , 5)) * 100;
prepImmRisk = squeeze(sum(sum(sum(pop(: , 5 , : , : , 1: 3), 2), 3), 4)) ./ ...
    (squeeze(sum(sum(sum(pop(: , 1 , : , : , 1: 3), 2), 3), 4)) + squeeze(sum(sum(sum(pop(: , 5 , : , : , 1: 3), 2), 3), 4))) * 100;

figure()
plot(t, prepImm, t , prepImmRisk )
xlim([startYear endYear])
title('PrEP use by risk group')
legend('Overall', 'High' , 'Medium' , 'Low' , 'Location', 'northwest')
xlabel('Year'); ylabel('Prevalence (%)')
%  
% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite(filename, [t(1: stepsPerYear :end), prepImm(1: stepsPerYear : end), prepImmRisk(1: stepsPerYear : end, :) ], 'PrEP', 'A2');
%         col_header = {'Year', 'Overall', 'High risk' , 'Medium risk' , 'Low risk'};
%         xlswrite(filename, col_header, 'PrEP', 'A1');

%% HIV prevalence among GC+ 
figure() 
gc_hivpos = squeeze(sum(sum(sum(pop(: , 2  , 2 , 2: sites, 1 : 3) , 3) , 4), 5)) ...
           ./ squeeze(sum(sum(sum(sum(pop(: , : , 2 , 2: sites , :), 2), 3) , 4) , 5)) * 100;
gc_hivtest = squeeze(sum(sum(sum(pop(: , 3  , 2 , 2: sites , 1 : 3) , 3) , 4), 5)) ...
           ./ squeeze(sum(sum(sum(sum(pop(: , : , 2 , 2: sites , :), 2), 3) , 4) , 5)) * 100;
gc_hivtreat = squeeze(sum(sum(sum(pop(: , 4 , 2 , 2: sites , 1 : 3) , 3) , 4), 5)) ...
           ./ squeeze(sum(sum(sum(sum(pop(: , : , 2 , 2: sites , :), 2), 3) , 4) , 5)) * 100;
gc_hivall = squeeze(sum(sum(sum(sum(pop(: , 2:4 , 2 , 2: sites , 1 : 3), 2), 3), 4), 5)) ...
           ./ squeeze(sum(sum(sum(sum(pop(: , : , 2 , 2: sites , :), 2), 3) , 4) , 5)) * 100;
plot(t, gc_hivpos, t, gc_hivtest, t, gc_hivtreat, t, gc_hivall)
legend('Infected' , 'Tested' , 'Treated', 'All HIV')
xlim([startYear endYear])
title('HIV prevalence among GC+')
xlabel('Year'); ylabel('Prevalence (%)')
% 


%% GC-HIV coinfection
gc_hivInf = sum(sum(pop(: , 2 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_hivTreated = sum(sum(pop(: , 4 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_hivTested = sum(sum(pop(: , 3 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_prepImm = sum(sum(pop(: , 5 , 2 , 2 : sites , :) , 4) , 5) ./ totalPop * 100;
gc_hivAll = sum(sum(sum(pop(: , 2 : 4 , 2 , 2 : sites , :) , 2),  4) , 5) ./ totalPop * 100;

figure()
plot(t , gc_hivInf , t , gc_hivTreated , t , gc_hivTested , t , gc_hivAll)
title('GC-HIV CoInfection')
legend('GC + HIV Infected' , 'GC + HIV Treated' , 'GC + HIV Tested' ,...
    'GC + All HIV', 'Location', 'northwest')
xlim([startYear endYear])
xlabel('Year'); ylabel('Prevalence (%)')
axis([tVec(1) tVec(end) 0 20])
% % 
% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite (filename, [t(1: stepsPerYear :end), gc_hivAll(1: stepsPerYear : end), gc_hivInf(1: stepsPerYear : end), ...
%         gc_hivTested(1: stepsPerYear: end), gc_hivTreated(1: stepsPerYear : end)], 'gc_hiv_coInf', 'A2');
%         col_header = {'Year', 'GC + All HIV', 'GC + HIV Infected', 'GC + HIV Tested', 'GC + HIV Treated'};
%         xlswrite(filename, col_header, 'gc_hiv_coInf', 'A1');
% %         
%% 
% gc_hivAll = sum(sum(sum(pop(: , 2 : 4 , 2 , 2 : sites , :) , 2),  4) , 5) ./ totalPop * 100;
% 
% figure()
% plot(t , gc_hivAll)
% title('GC-HIV CoInfection')
% %legend('GC + HIV Infected' , 'GC + HIV Treated' , 'GC + HIV Tested' ,...
% %    'GC + All HIV', 'Location', 'northeast')
% xlabel('Year'); ylabel('Prevalence (%)')
%axis([tVec(1) tVec(end) 0 20])

%% GC Prevalence by HIV status
gc_hivInf = sum(sum(pop(: , 2 , 2 , 2 : sites , :) , 4) , 5) ./...
    sum(sum(sum(pop(: , 2 , : , : , :) , 3) , 4) , 5) * 100;
gc_hivTreated = sum(sum(pop(: , 4 , 2 , 2 : sites , :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 4 , : , : , :) , 3) , 4) , 5) * 100;
gc_hivTested = sum(sum(pop(: , 3 , 2 , 2 : sites , :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 3 , : , : , :) , 3) , 4) , 5) * 100;
% gc_prepImm = sum(sum(pop(: , 5 , 2 , 2 : sites , :) , 4) , 5) ./ ...
%     sum(sum(sum(pop(: , 5 , : , : , :) , 3) , 4) , 5) * 100;
gc_hivNeg = sum(sum(pop(: , 1 , 2 , 2 : sites , :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 1 , : , : , :) , 3) , 4) , 5) * 100;

figure()
subplot(2 , 1, 1)
plot(t , gc_hivInf , t, gc_hivTested, t , gc_hivTreated , t , gc_hivNeg)
title('GC prevalence by HIV Status')
legend('HIV Infected' ,  'HIV Tested' , 'HIV Treated' ,...
    'HIV Negative')
xlim([startYear endYear])
xlabel('Year'); ylabel('Prevalence (%)')

Rgc_hivInf = sum(sum(pop(: , 2 , 2 , 2  , :) , 4) , 5) ./...
    sum(sum(sum(pop(: , 2 , : , : , :) , 3) , 4) , 5) * 100;
Rgc_hivTreated = sum(sum(pop(: , 4 , 2 , 2 , :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 4 , : , : , :) , 3) , 4) , 5) * 100;
Rgc_hivTested = sum(sum(pop(: , 3 , 2 , 2, :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 3 , : , : , :) , 3) , 4) , 5) * 100;
Rgc_hivNeg = sum(sum(pop(: , 1 , 2 , 2 , :) , 4) , 5) ./ ...
    sum(sum(sum(pop(: , 1 , : , : , :) , 3) , 4) , 5) * 100;

subplot(2 , 1, 2)
plot(t , Rgc_hivInf , t, Rgc_hivTested, t , Rgc_hivTreated , t , Rgc_hivNeg)
title('Rectal GC prevalence by HIV Status')
% legend('HIV Infected' ,  'HIV Tested' , 'HIV Treated' ,...
%     'HIV Negative')
xlim([startYear endYear])
xlabel('Year'); ylabel('Prevalence (%)')


% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite (filename, [t(1: stepsPerYear :end), gc_hivInf(1: stepsPerYear : end), gc_hivTested(1: stepsPerYear : end), ...
%        gc_hivTreated(1: stepsPerYear : end),  gc_hivNeg(1: stepsPerYear : end)], 'GC-by-HIVstat', 'A2');
%        col_header = {'Year', 'HIV Infected' , 'HIV Treated' , 'HIV Tested',  'HIV Negative'};
%        xlswrite(filename, col_header , 'GC-by-HIVstat', 'A1'); 
% %        
%% GC and HIV susceptibles
% figure()
% gc_Sus = sum(sum(sum(sum(pop(: , 1 : hivStatus , 1 , 1 : 4 , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
% subplot(2 , 1 , 1)
% plot(t , gc_Sus)
% title('GC-Susceptible')
% xlabel('Year'); ylabel('Proportion (%)')
% axis([tVec(1) tVec(end) 0 100])
% 
% hiv_sus = sum(sum(sum(sum(pop(: , 1 , 1 : stiTypes , 1 : sites , 1 : risk) , 2) , 3), 4) , 5) ./ totalPop * 100;
% subplot(2 , 1 , 2)
% plot(t , hiv_sus)
% title('HIV-Susceptible')
% xlabel('Year'); ylabel('Proportion (%)')
% axis([tVec(1) tVec(end) 0 100])

%%
% % HIV
% totalPop_Risk = squeeze(sum(sum(sum(pop(: , : , : , : , :), 2), 3) , 4));
% 
% hivInf = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 2 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
% hivTreated = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 4 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
% hivTested = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 3 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
% prepImm = bsxfun(@rdivide , squeeze(sum(sum(pop(: , 5 , : , : , :) , 3) , 4)) , totalPop_Risk) * 100;
% 
% hivArray = {hivInf, hivTested, hivTreated};
% hivPlotTitle = {'Infected' , 'Tested' , 'Treated'};
% 
% for i = 1 : size(hivArray,2)
%     figure()
%     plot(t , hivArray{i})
%     title(['HIV ' , hivPlotTitle{i}])
%     legend('High Risk' , 'Medium Risk' , 'Low Risk')
%     xlim([startYear endYear])
%     xlabel('Year'); ylabel('Prevalence (%)')
% end
%% Distribution of HIV by risk group
riskDistHiv = squeeze(sum(sum(sum(pop(: , 2 : 4 , : , : , :), 2), 3), 4)) ./ ...
    squeeze(sum(sum(sum(sum(pop(: , 2 : 4 , : , : , :), 2), 3), 4),5));
figure()
area(t , riskDistHiv)
xlim([startYear endYear])
legend({'High Risk' , 'Medium Risk' , 'Low Risk'}, 'Location', 'southeast')
xlabel('Year'); ylabel('Proportion')
title('Distribution of HIV by risk group')


riskDistHiv = squeeze(sum(sum(sum(pop(: , 2  , : , : , :), 2), 3), 4)) ./ ...
    squeeze(sum(sum(sum(sum(pop(: , 2 , : , : , :), 2), 3), 4),5));
figure()
area(t , riskDistHiv)
xlim([startYear endYear])
legend({'High Risk' , 'Medium Risk' , 'Low Risk'}, 'Location', 'southeast')
xlabel('Year'); ylabel('Proportion')
title('Distribution of undiagnosed HIV by Risk Group')

%% HIV Prevalence in Low Risk
% lowRiskHiv = squeeze(sum(sum(pop(: , 2 : 4 , : , : , 3), 3), 4)) ./ ...
%     squeeze(sum(sum(sum(sum(pop(: , :, : , : , 3), 2), 3), 4))) * 100;
% figure()
% area(t , lowRiskHiv)
% legend('Infected' , 'Tested' , 'Treated')
% xlabel('Year'); ylabel('Percent')
% title('HIV Prevalence in Low Risk')
%%
%GC
% totalPop_Risk = squeeze(sum(sum(sum(pop(: , : , : , : , :), 2), 3) , 4));
% gcInf = bsxfun(@rdivide , squeeze(sum(sum(pop(: , : , 2 , 2 : 4 , :) , 2) , 4)) , totalPop_Risk) * 100;
% figure()
% plot(t , gcInf)
% title('GC Prevalence by Risk Group')
% legend('High Risk' , 'Medium Risk' , 'Low Risk')
% xlabel('Year'); ylabel('Prevalence (%)')

%% HIV incidence by risk
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values

figure()

for i=1:risk
    hivInf_noSti = zeros(3,length(newHiv)/stepsPerYear);
    hivInf_sti = zeros(3,length(newHiv)/stepsPerYear);
    hivSus_noSti = zeros(3,length(newHiv)/stepsPerYear);
    hivSus_sti = zeros(3,length(newHiv)/stepsPerYear);
    
    hivInf_noSti(i,:) = annlz(squeeze(newHiv(:,1,1,i))); % new HIV infections with no STI
    hivInf_sti(i,:) = annlz(squeeze(sum(newHiv(:,2,2:4,i),3))); % new HIV infections with STI based on risk group
    
    hivSus_noSti(i,:) = annlz(pop(1:end-1,1,1,1,i))./stepsPerYear; % number susceptible to HIV over year without STI
    hivSus_sti(i,:) = annlz(sum(pop(1:end-1,1,2,2:4,i),4))./stepsPerYear; % number susceptible to HIV over year with STI
    
    hivInc_noSti = (hivInf_noSti(i,:)./hivSus_noSti(i,:)).*100000;
    hivInc_sti = (hivInf_sti(i,:)./hivSus_sti(i,:)).*100000;

%     subplot(1,3,1)
%     plot(t(1:stepsPerYear:end-1),hivInc_noSti)
%     xlabel('Year'); ylabel('Incidence Per 100,000'); 
%     axis([startYear , endYear , 0 , 5000]);   
%     title('HIV Incidence among those w/o an STI')
% %     legend('High Risk' , 'Medium Risk' , 'Low Risk')
%     hold all;
%     subplot(1,3,2)
%     plot(t(1:stepsPerYear:end-1),hivInc_sti)
%     title('HIV Incidence among those with an STI')
% %     legend('High Risk' , 'Medium Risk' , 'Low Risk')
%     xlabel('Year'); ylabel('Incidence Per 100,000');
%     axis([startYear , endYear , 0 , 5000]);
%     hold all;
%     subplot(1,3,3)
    hivInc_tot_risk = ((hivInf_noSti(i,:)+hivInf_sti(i,:))./(hivSus_noSti(i,:)+hivSus_sti(i,:))).*100000;
    plot(t(1:stepsPerYear:end-1),hivInc_tot_risk)
    xlabel('Year'); ylabel('Incidence Per 100,000'); 
%     axis([startYear , endYear , 0 , 5000]);   
    title('HIV Incidence by risk group')
    legend('High Risk' , 'Medium Risk' , 'Low Risk')
    hold all;
end

%% HIV incidence

hivInf_noSti = annlz(squeeze(sum(newHiv(:,1,1,:),4))); % new HIV infections with no STI
hivInf_sti = annlz(squeeze(sum(sum(newHiv(:,2,2:4,:),3),4))); % new HIV infections with STI based on risk group

hivSus_noSti = annlz(squeeze(sum(pop(1:end-1,1,1,1,:),5)))./stepsPerYear; % number susceptible to HIV over year without STI
hivSus_sti = annlz(sum(sum(pop(1:end-1,1,2,2:4,:),4),5))./stepsPerYear; % number susceptible to HIV over year with STI

hivInc_noSti = (hivInf_noSti./hivSus_noSti).*100000;
hivInc_sti = (hivInf_sti./hivSus_sti).*100000;
%
figure()
% subplot(1,3,1)
% plot(t(1:stepsPerYear:end-1),hivInc_noSti)
% xlabel('Year'); ylabel('Incidence Per 100,000'); 
% axis([startYear , endYear , 0 , 2000]);   
% title('HIV Incidence among those w/o an STI')
% subplot(1,3,2)
% plot(t(1:stepsPerYear:end-1),hivInc_sti)
% title('HIV Incidence among those with an STI')
% xlabel('Year'); ylabel('Incidence Per 100,000');
% axis([startYear , endYear , 0 , 2000]);
% subplot(1,3,3)


hivInc_tot = ((hivInf_noSti+hivInf_sti)./(hivSus_noSti+hivSus_sti)).* 100000;
plot(t(1:stepsPerYear:end-1),hivInc_tot)
xlabel('Year'); ylabel('Incidence Per 100,000'); 
%axis([startYear , endYear , 0 , 2000]);   
title('HIV Incidence, any risk')
% 
% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite (filename, [t(1: stepsPerYear :end-1), hivInc_tot'], 'HIV_inc', 'A2');
%        col_header = {'Year', 'Overall incidence' };
%        xlswrite(filename, col_header , 'HIV_inc', 'A1'); 
%        
       %% HIV incidence counts 
%        figure()
% plot(t(1:stepsPerYear:end-1), (hivInf_noSti+hivInf_sti))
% title('Annual new HIV cases')
% xlabel('Year'); ylabel('Number of new cases');
% 
%        
%% Overall GC incidence
siteNames = ["none","rectal","urethral","pharyngeal"];
gcInf_r = annlz(squeeze(sum(sum(newSti(:,:,2,:),2),4))); 
gcInf_u = annlz(squeeze(sum(sum(newSti(:,:,3,:),2),4)));
gcInf_p = annlz(squeeze(sum(sum(newSti(:,:,4,:),2),4)));
gcSus = annlz(sum(sum(pop(1:end-1,:,1,1,:),2),5))./stepsPerYear;
gcInc_r = (gcInf_r ./ gcSus).* 100000;
gcInc_u = (gcInf_u ./ gcSus).* 100000;
gcInc_p = (gcInf_p ./ gcSus).* 100000;

figure()
subplot(3, 1, 1)
plot(t(1:stepsPerYear:end-1),gcInc_r)
xlabel('Year'); ylabel('Incidence Per 100,000');
% axis([startYear , endYear , 0 , 50000]);
title('Rectal GC Incidence, any risk')

subplot(3, 1, 2)
plot(t(1:stepsPerYear:end-1),gcInc_u)
xlabel('Year'); ylabel('Incidence Per 100,000');
title('Urethral GC Incidence, any risk')

subplot(3, 1, 3)
plot(t(1:stepsPerYear:end-1),gcInc_p)
xlabel('Year'); ylabel('Incidence Per 100,000');
title('Pharyngeal GC Incidence, any risk')
% 
% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite (filename, [t(1: stepsPerYear :end-1), gcInc_r', gcInc_u', gcInc_p'], 'GC_inc', 'A2');
%        col_header = {'Year', 'Rectal GC' , 'Urethral GC', 'Pharyngeal GC'};
%        xlswrite(filename, col_header , 'GC_inc', 'A1'); 
% %        
% %% GC incidence
% figure()
% siteNames = ["none","rectal","urethral","pharyngeal"];
% for i=2:sites
% 
%         gcInf_site = annlz(squeeze(sum(newSti(:,1,i,:), 4))); 
%         gcSus_site = zeros(1,length(newSti)/stepsPerYear);
%         %for k=1:sites
%         %    if k==i
%         %        continue;
%         %    else
%         %        gcSus = gcSus + annlz(sum(pop(1:end-1,:,1,k,j),2))./stepsPerYear; % number susceptible to GC over time
%         %    end
%         %end
%         gcSus_site = annlz(sum(sum(pop(1:end-1,:,1,1,:),2),5)) ./stepsPerYear;
%         gcInc_site = (gcInf_site./gcSus_site).*100000;
% 
%         subplot(1,3,i-1);
%         plot(t(1:stepsPerYear:end-1),gcInc_site)
%         xlabel('Year'); ylabel('Incidence Per 100,000');
%         axis([startYear , endYear , 0 , 10000]);
%         title(siteNames(i)+' Incidence')
%         %legend('High Risk' , 'Medium Risk' , 'Low Risk')
%         hold all;
% 
%         j = 'C2';
%         if i == 3 
%             j = 'D2';
%         end
%         if i == 4
%             j = 'E2';
%         end
% % filename = 'Result_psOnly.xlsx';
% %      xlswrite (filename, [gcInc_site'], 'GC_inc', j);
% %      col_header = {'Year' , 'Overall Inc', 'Rectal GC' , 'Urethral GC', 'Pharyngeal GC'};
% %      xlswrite(filename, col_header , 'GC_inc', 'A1');   
% end
% 
%        
% %% GC incidence by HIV stat
gcInf_N = annlz(squeeze(sum(sum(newSti(:,1,2:4,:),3),4))); 
gcSus_N = annlz(sum(sum(pop(1:end-1,1,1,1,:),2),5))./stepsPerYear;
gcInc_N = (gcInf_N ./ gcSus_N).*100000;

gcInf_I = annlz(squeeze(sum(sum(newSti(:,2,2:4,:),3),4))); 
gcSus_I = annlz(sum(sum(pop(1:end-1,2,1,1,:),2),5))./stepsPerYear;
gcInc_I = (gcInf_I ./ gcSus_I).*100000;

gcInf_S = annlz(squeeze(sum(sum(newSti(:,3,2:4,:),3),4))); 
gcSus_S = annlz(sum(sum(pop(1:end-1,3,1,1,:),2),5))./stepsPerYear;
gcInc_S = (gcInf_S ./ gcSus_S).*100000;

gcInf_T = annlz(squeeze(sum(sum(newSti(:,4,2:4,:),3),4))); 
gcSus_T = annlz(sum(sum(pop(1:end-1,4,1,1,:),2),5))./stepsPerYear;
gcInc_T = (gcInf_T ./ gcSus_T).*100000;

gcInf_P = annlz(squeeze(sum(sum(newSti(:,5,2:4,:),3),4))); 
gcSus_P = annlz(sum(sum(pop(1:end-1,5,1,1,:),2),5))./stepsPerYear;
gcInc_P = (gcInf_P ./ gcSus_P).*100000;

figure()
plot(t(1:stepsPerYear:end-1), gcInc_N, t(1:stepsPerYear:end-1), gcInc_I, t(1:stepsPerYear:end-1), gcInc_S, t(1:stepsPerYear:end-1), gcInc_T, t(1:stepsPerYear:end-1), gcInc_P)
legend({'Negative', 'Infected', 'Tested', 'Treated', 'PREP'},'Location', 'northwest');
xlabel('Year'); ylabel('Incidence Per 100,000');
% axis([startYear , endYear , 0 , 100000]);
title('GC Incidence among HIV positive')

figure()
plot(t(1:stepsPerYear:end-1), gcInc_N, t(1:stepsPerYear:end-1), gcInc_P)
legend({'Not on PrEP', 'on PrEP'},'Location', 'northwest');
xlabel('Year'); ylabel('Incidence Per 100,000');
% axis([startYear , endYear , 0 , 100000]);
title('GC Incidence among HIV negative MSM')

% filename = 'Result_basecase_RoutPrep3.xlsx';
%      xlswrite (filename, [t(1: stepsPerYear :end-1), gcInc_N', gcInc_I', gcInc_S', gcInc_T', gcInc_P'], 'GC_inc_by_HIV', 'A2');
%        col_header = {'Year', 'Negative' , 'Infected' , 'Tested', 'Treated', 'PrEP'};
%        xlswrite(filename, col_header , 'GC_inc_by_HIV', 'A1'); 

%% GC treatment and HIV testing due to PS 
% 
% gcPs = annlz(squeeze(sum(newGcPs(:, 1, 2, 2, :), 5)));
% 
%         filename = 'Result_psAndHIV.xlsx'
% xlswrite (filename, [t(1: stepsPerYear :end-1), gcPs(1:  end)'], 'GC_PsHIV', 'A2');
%         col_header = {'Year', 'Rectal' };
%         xlswrite(filename, col_header , 'GC_PsHIV', 'A1');
% 
% gcHivPs = annlz(squeeze(sum(sum(newGcPs(:, 2, 2, 2:4, :), 4), 5)));
% xlswrite (filename, [t(1: stepsPerYear :end-1), gcHivPs(1:  end)'], 'GC_PsHIV', 'D2');
%         col_header = {'Year', 'HIV'};
%         xlswrite(filename, col_header , 'GC_PsHIV', 'D1');