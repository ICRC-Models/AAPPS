function [negSumLogL] = mainCalibrate(paramSet)

%% Load pre-set parameters
load('genParams')
load('gcParams')
load('gcHivParams')
load('partnerParams')
load('calibTargets')
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
% popInitial(2 , 2 , 2:4 , 1) = 40 * 1.8;
popInitial(2 , 1 , 1 , 2) = 315 * 1.8;
popInitial(2 , 2 , 2:4 , 2) = 35 * 1.8;
% condUse = [0.23 , 0.25 , 0.44]; % condom usage by risk group (high, med, low)
% condUse = [0.23, 0.29, 0.4]; % condom usage by risk group (high, med, low) TEST!!!
riskVec = zeros(risk , 1);
popNew = popInitial .* 0.01;
popVec = zeros(length(tspan), hivStatus , stiTypes , sites , risk);

% partners = c;
kDie = 0.0018;
kBorn = 7 * kDie;

% sympref('HeavisideAtOrigin' , 1);

%% Load calibrated parameters
pIdx = load([pwd , '\' , 'pIdx_21June19.dat']); % load indices for parameters being calibrated
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx) % Load info into paramsSub cell array for the parameters being calibrated
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end
% Partners
if any(1 == pIdx)
    idx = find(1 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rh = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rl = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersAnal(1:hivStatus , 1:risk) = [paramSet(rh) , paramSet(rm) , paramSet(rl)];
end
if any(2 == pIdx)
    idx = find(2 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rh = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rl = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersOral(1:hivStatus , 1:risk) = [paramSet(rh) , paramSet(rm) , paramSet(rl)];
end
partners(:,:,1) = partnersAnal;
partners(:,:,2) = partnersOral;

% Acts
if any(3 == pIdx)
    idx = find(3 == pIdx);
    acts = paramSet(paramsSub{idx}.inds(:));
end

% HIV anal transmission
if any(4 == pIdx)
    idx = find(4 == pIdx);
    perActHiv = paramSet(paramsSub{idx}.inds(:));
end

% GC transmission probabilities by site and mode of transmission
% 'x' indicates that transmission does not occur by this route
if any(5 == pIdx)
    idx = find(5 == pIdx);
    recUre = paramSet(paramsSub{idx}.inds(1));
    recPha = paramSet(paramsSub{idx}.inds(2));
    ureRec = paramSet(paramsSub{idx}.inds(3));
end
perAct_Anal = [0 ,    recUre,  recPha;... % rectal -> (x, urethral , pharyngeal)
               ureRec,  0,      0 ; ... % urethral -> (rectal , x , x)
               0,     0,      0] .* .35; % pharyngeal -> (rectal , x , x)
if any(6 == pIdx)
    idx = find(6 == pIdx);
    urePha = paramSet(paramsSub{idx}.inds(1));
    phaRec = paramSet(paramsSub{idx}.inds(2));
    phaUre = paramSet(paramsSub{idx}.inds(3));
end     
perAct_Oral = [0 ,   0,    0 ; ... % rectal -> (x , x , x)
               0 ,   0,    urePha ;... % urethral -> (x , x , pharyngeal)
               phaRec, phaUre, 0 ] .* .3 ; % pharyngeal -> (rectal , urethral , x)

% Condom use 
if any(7 == pIdx)
    idx = find(7 == pIdx);
    cAssortTarget = paramSet(paramsSub{idx}.inds(:))';
end  
if any(8 == pIdx)
    idx = find(8 == pIdx);
    cAssort_init = paramSet(paramsSub{idx}.inds(:))';
end  

%% Interventions

% partner services
% p_symp2 = p_symp .* [1; 2; 1];
% psTreatMatTarget = (1 - exp(-((p_ps + .7) .* p_symp )));

psTreatMatTarget = (p_ps + .3) .* p_symp .* 0;

% routine treatment scale-up
routine1TreatMat_init =  ones(3,4) .* 0.001; 
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
% cAssortTarget = [.6, .5, .5]; % Target plateau values for condom use
% cAssort_init = [.2, .25, .3]; % Initial condom use values
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

%% ODE solver
error = 0;

newHiv = zeros(length(tspan)-1,stiTypes,sites,risk);
newSti = zeros(length(tspan)-1,hivStatus,sites,risk);
newGcPs = zeros(length(tspan)-1, hivStatus, stiTypes, sites, risk);

popVec(1,:,:,:,:) = popInitial;
popIn = popInitial;

% disp('Running...')
% tic

for time = 1:length(tspan)-1
    tspanStep = [tspan(time) , tspan(time + 1)]; % evaluate diff eqs over one time interval
    
    [t , pop, newHiv(time,:,:,:), newSti(time,:,:,:), newGcPs(time, :, :, :, :)] = ode4xtra(@(t , pop) mixInfect_Hiv_rout_calibrate(t , pop , ...
        hivStatus , stiTypes , sites , risk , ...
        kDie , kBorn , gcClear , rout1Scale, d1_routineTreatMat , routine1TreatMat_init , ...
        rout2Scale, d2_routineTreatMat , routine2TreatMat_init, ...
        p_symp , fScale ,fScale_HivScreen , d_psTreatMat , kDiagTreat , ...
        kHivScreen_init , d_kHivScreen , d_muHiv, muHiv_init,  hTscale, d_hTreat, hTreat_init, partners , acts , riskVec ,...
        cScale, d_cAssort, cAssort_init , d_hAssort , hScale , hAssort_init, d_cotest, cotestInit, cotestStartInd, ...
        d_rAssort , rScale, rAssort_init, perAct_Anal, perAct_Oral, perActHiv, tVec) , tspanStep , popIn);

    popIn = reshape(pop(end,:) , [hivStatus , stiTypes , sites , risk]);
    popVec(time+1,:,:,:,:) = popIn;
    
    if any(popVec(:,:,:,:,:) <  0)
        disp('Error!')
        error = 1;
        break
    end

end
% disp('Finished solving')
% toc
% disp(' ')

if error
    negSumLogL = 1000000000;
else
    negSumLogL = likeFun(popVec , newHiv , ...
        hivPrev_obs , hivInc_obs , gcRecPrev_obs , gcUrePrev_obs , gcPhaPrev_obs , ...
        hivStatus , stiTypes , sites , risk , stepsPerYear , startYear);
end
