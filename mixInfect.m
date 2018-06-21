function dPop = mixInfect(t , pop , hivStatus , stiTypes , sites , ...
    risk , popNew , kDie , kBorn , gcClear , d_routineTreatMat , routineTreatMat_init , ...
    p_symp , fScale ,fScale_HivScreen , d_psTreatMat , kDiagTreat , ...
    kHivScreen_init , d_kHivScreen , kHivTreat , partners , acts , riskVec ,...
    condUse , d_hAssort , hScale , hAssort_init , rAssort , tVec)
%%
sumall = @(x) sum(x(:));
pop(5 , : , : , :) = 0;
%% Scale-up treatments
% kTest_p_ps = interp1(tVec , kTest_p_ps , t);
% kTreat_p_ps = interp1(tVec , kTreat_p_ps , t);
% kTest_u_ps = interp1(tVec , kTest_u_ps , t);
% kTreat_u_ps = interp1(tVec , kTreat_u_ps , t);
% kTest_r_ps = interp1(tVec , kTest_r_ps , t);
% kTreat_r_ps = interp1(tVec , kTreat_r_ps , t);
% k_gcPS(1 : 3) = psProp .* interp1(tVec , k_gcPS , t); % partner services treatment rate (proportion of GC treated by PS) scale-up

fScale = interp1(tVec , fScale , t);
psTreatMat = fScale .* d_psTreatMat;
routineTreatMat = fScale .* d_routineTreatMat + routineTreatMat_init;
fScale_HivScreen = interp1(tVec , fScale_HivScreen , t);
kHivScreen = fScale_HivScreen .* d_kHivScreen + kHivScreen_init;
hScale = interp1(tVec , hScale , t);
hAssort = hScale .* d_hAssort + hAssort_init;
%%
pop = reshape(pop , [hivStatus , stiTypes , sites , risk]);

% GC transmission probabilties by site and mode of transmission
% 'x' indicates that transmission does not occur by this route

perAct_Anal = [0 , 0.84 , 0 ;... % rectal -> (x, urethral , pharyngeal)
    0.243 , 0  0.0865 ; ... % urethral -> (rectal , x , pharyngeal)
    0 , 0.62 , 0]; % pharyngeal -> (x , urethral , x)

perAct_Oral = [0 , 0 , 0 ; ... % rectal -> (x , x , x)
    0 , 0 , 0.62 ;... % urethral -> (x , x , pharyngeal)
    0 , 0.62 , 0] .* 0.01; % pharyngeal -> (x , urethral , x)

perActHiv = 0.82 * 10 ^ -2; % anal transmission

acts = acts ./ 2; % half anal, half oral for testing

% transmission probabilities assume independence of potential transmission events
perPartner_Anal = min(1 - (1 - perAct_Anal) .^ acts , 0.99999); % GC anal transmission probability
perPartner_Oral = min(1 - (1 - perAct_Oral) .^ acts , 0.99999); % GC oral transmission probability
perPartnerHiv = 1 - (1 - perActHiv) .^ acts; % HIV transmission probability

% assortativity matrices by HIV status and risk
hAssortMat = eye(3);% (1) HIV-/unknown status, (2) HIV+/Tested, (3) PrEP
rAssortMat = eye(risk); % (1) High, (2) Medium, (3) Low

% probability of individuals in HIV state h mixing with individuals in HIV state hh
mixMatHiv = zeros(3 , 3);
% probability of individuals in HIV state h and risk group r mixing with individuals in risk group rr
mixMatRisk = zeros(3 , risk , risk); 
% probability of individuals in HIV state h and risk group r mixing with individuals in HIV state hh and risk group rr
mixMat = zeros(3 , 3 , risk , risk);

byHiv = squeeze(sum(sum(pop , 2) , 3)); % sum across sti types and sites --> dimensions: [hivStatus(5) x risk(3)]
% for mixing purposes, HIV-negatives and HIV-infectious who do not know
% their own status are combined into one group
hivNeg = sum(sum(byHiv(1 : 2 , :))); % negative, infectious (own status unknown) --> [1 x risk(3)]
hivPos = sum(sum(byHiv(3 : 4 , :))); % tested, treated (HIV-positive with known status) --> [1 x risk(3)]
hivImm = sum(byHiv(5 , :)); %PrEP/HIV Immune --> [1 x risk(3)]
total = sumall(byHiv); % total population

% Mixing matrix by HIV status
% random mixing component
mixMatHiv(: , 3) = hivImm / total .* (1 - hAssort); % proportion of partnerships from random mixing with PrEP
mixMatHiv(: , 2) = hivPos / total .* (1 - hAssort); % proportion of partnerships from random mixing with HIV positive
mixMatHiv(: , 1) = hivNeg / total .* (1 - hAssort); % proportion of partnerships from random mixing with HIV-/status unknown

% HIV status mixing matrix resulting from combination of assortative and random mixing
% by HIV status
mixMatHiv = mixMatHiv + hAssort .* hAssortMat; % random + assortative mixing by HIV status

%  
byHiv_Risk = squeeze(sum(sum(pop , 2) , 3)); % [hivStatus x sites]

% Risk distribution by HIV status
hivNeg_Risk = zeros(1 , risk); 
hivPos_Risk = zeros(1 , risk); 
hivImm_Risk = zeros(1 , risk); 

if hivNeg > 0
    hivNeg_Risk = sum(byHiv_Risk(1 : 2 , :) , 1) ./ hivNeg; % proportion of negatives (1) and infectious (2) in each risk group 
end

if hivPos > 0
    hivPos_Risk = sum(byHiv_Risk(3 : 4 , :) , 1) ./ hivPos; % proportion of tested (3) and treated (4) in each risk group
end

if hivImm > 0
    hivImm_Risk = sum(byHiv_Risk(5 , :) , 1) ./ hivImm; % proportion of PrEP takers (5) in each risk group
end

% Mixing matrix by risk group for individuals in each HIV group
% For random mixing, the proportion of partners from risk group "r" is
% equal to the proportion of people in the population belonging to risk
% group "r"
for r = 1 : risk
    mixMatRisk(1 , : , r) = hivNeg_Risk(r) .* (1 - rAssort); % proportion of contacts from random mixing by risk group with HIV negative and HIV infectious
    mixMatRisk(2 , : , r) = hivPos_Risk(r) .* (1 - rAssort); % proportion of contacts from random mixing by risk group with HIV positive
    mixMatRisk(3 , : , r) = hivImm_Risk(r) .* (1 - rAssort); % proportion of contacts from random mixing by risk group with HIV immune
end

for h = 1 : 3
    mixMatRisk(h , : , :) = squeeze(mixMatRisk(h , : , :)) + rAssort .* rAssortMat; % adding assortative mixing by risk
end

% Mixing matrix by HIV status and risk group
% mixMat(h , hh , rr , r)
% probability of a person from HIV group h and risk group r mixing with a
% person from HIV group hh and risk group rr
for h = 1 : 3
    for hh = 1 : 3
        mixMat(h , hh , : , :) = mixMatHiv(h , hh) .* mixMatRisk(hh , : , :);
    end
end

% Adjust partners
% Ensures number of partners that group a has coming from group b is equal 
% to the number of partners group b has coming from group a
% dim(partners) = [hiv state x risk x act type]
partners_hivImm = squeeze(mean(partners(5 , : , :) , 1)); % [risk x act type]
partners_hivPos = squeeze(mean(partners(3 : 4 , : , :) , 1)); % [risk x act type]
partners_hivNeg = squeeze(mean(partners(1 : 2 , : , :) , 1)); % [risk x act type]

partners_anal = [partners_hivNeg(: , 1) , partners_hivPos(: , 1) , partners_hivImm(: , 1)];
partners_oral = [partners_hivNeg(: , 2) , partners_hivPos(: , 2) , partners_hivImm(: , 2)];

% increase partners in high risk groups
partners_anal(1 , :) = 1.5 .* partners_anal(1, :);
partners_oral(1 , :) = 1.5 .* partners_oral(1, :);

fracPop = [hivNeg_Risk ; hivPos_Risk ; hivImm_Risk];
% fracAct = zeros(size(partners));
% actTypes = 2;
% for f = 1 : actTypes % (1) anal , (2) oral
%     fracAct(: , : , f) = partners(: , : , f) ./ sum(partners , 3);
% end
hivNegPosImm = 3;
adjustFac_anal = zeros(hivNegPosImm , hivNegPosImm , risk , risk);
adjustFac_oral = adjustFac_anal;
% ratio of HIV+ with infection at site : HIV- with infection at site
% obtain adjustment factors for partnerships
for r = 1 : risk
    for rr = 1 : risk
        for h = 1 : hivNegPosImm
            for hh = 1 : hivNegPosImm
                frac = 0;
                if fracPop(hh , rr) > 10 ^ -6
                    frac = fracPop(h , r) / fracPop(hh , rr);
                end
                if sum(partners_anal(rr , hh) * mixMat(hh , h , r , rr)) > 10 ^ -6
                    if sum(partners_anal(r , h) * mixMat(h , hh , rr , r)) > 10 ^ -6
                        adjustFac_anal(h , hh , rr , r) = ...
                            sum(partners_anal(r , h) * mixMat(h , hh , rr , r)) ...
                            / sum(partners_anal(rr , hh) * mixMat(hh , h , r , rr)) * frac;
                    end
                    if sum(partners_oral(r , h) * mixMat(h , hh , rr , r)) > 10 ^ -6
                        adjustFac_oral(h , hh , rr , r) = ...
                            sum(partners_oral(r , h) * mixMat(h , hh , rr , r)) ...
                            / sum(partners_oral(rr , hh) * mixMat(hh , h , r , rr)) * frac;
                    end
                end
            end
        end
    end
end
%%
analPartnersAdj = zeros(3 , 3 , risk , risk);
oralPartnersAdj = analPartnersAdj;
theta = 0.5;

for r = 1 : risk
    for rr = 1 : risk
        for h = 1 : hivNegPosImm
            for hh = 1 : hivNegPosImm
                frac = 0;
                if fracPop(hh , rr) > 10 ^ -6
                    frac = fracPop(h , r) / fracPop(hh , rr);
                end
                % Adjust anal partners
                if adjustFac_anal(h , hh , rr , r) ~= 0
                    analPartnersAdj(h , hh , rr , r) = partners_anal(r , h) ...
                        .* frac ^ theta ...
                        * adjustFac_anal(h , hh , rr , r) .^ -(1 - theta);
               
                    analPartnersAdj(hh , h , r , rr) = partners_anal(rr , hh) ...
                        .* frac ^ -(1 - theta)...
                        * adjustFac_anal(h , hh , rr , r) .^ theta;
                end
                % Adjust oral partners
                if adjustFac_oral(h , hh , rr , r) ~= 0
                    oralPartnersAdj(h , hh , rr , r) = partners_oral(r , h) ...
                        .* frac ^ theta...
                        * adjustFac_oral(h , hh , rr , r) .^ -(1 - theta);

                    oralPartnersAdj(hh , h , r , rr) = partners_oral(rr , hh) ...
                        .* frac ^ -(1 - theta)...
                        * adjustFac_oral(h , hh , rr , r) .^ theta;
                end
            end
        end
    end
end

%%
% STIs (other than HIV) per year per partnership transmission. Per year per
% event rate of infection (site i -> site ii)
% per year inf calculated as -log(1-per_Partner_per_Year_Transmission)
perYear_AnalInf = zeros(hivNegPosImm , stiTypes , 3 , risk , sites , sites);
perYear_OralInf = perYear_AnalInf;
popSubs = [sum(byHiv(1 : 2 , :)) ; sum(byHiv(3 : 4 , :)) ; byHiv(5 , :)]; % pop subtotals by hiv status. hivNeg includes HIV+ who think they are negative
hivNegPop = sum(pop(1 : 2 , : , : , :) , 1); % (1) HIV-negative and HIV-infected/status unknown group
hivPosActual = sum(pop(2 , : , : , :) , 1); % HIV-infected/status unknown (a subset of (1))
hivPosPop = sum(pop(3 : 4 , : , : , :) , 1); % (2) HIV-tested and HIV-treated group
hivImmPop = pop(5 , : , : , :); % (3) PrEP group 
pops = {hivNegPop , hivPosPop , hivImmPop};
% Anal transmission probability for non-HIV STIs
perPartner_Anal = [zeros(1 , 3) ; perPartner_Anal];
perPartner_Anal = [zeros(4 , 1) , perPartner_Anal];
% Oral transmission probability for non-HIV STIs
perPartner_Oral = [zeros(1 , 3) ; perPartner_Oral];
perPartner_Oral = [zeros(4 , 1) , perPartner_Oral];

for h = 1 : hivNegPosImm
    for ty = 1 : stiTypes
        for r = 1 : risk
            for s = 1 : sites
                for ss = 1 : sites
                    popSubtotal = popSubs(h , r);
                    popGroup = squeeze(pops{h});
                    if popSubtotal > 0
                        % if h ~= 1 (not hiv negative group) 
                        contactProb = 0;
                        if popGroup(ty , ss , r) > 10 ^ -6
                            contactProb = popGroup(ty , ss , r) ./ popSubtotal; % proportion of pop with sti 
                        end
                        contactProbHiv = 0; % initialize/for HIV positive or HIV-immune group (HIV status known)
                        % per partner transmission probabilities
                        joint = perPartner_Anal(s , ss) * perPartnerHiv; % sti and hiv
                        p_hiv = perPartnerHiv; % hiv
                        
                        if h == 1 && hivPosActual(1 , ty , s , r) > 10 ^ -6
                            contactProbHiv = hivPosActual(1 , ty , s , r) ./ popSubtotal; % proportion of "hiv-neg" that are really hiv-pos
                        end
                        
                        
                        perYear_AnalInf(h , ty , 1 , r , s , ss) = ... % probability of getting STI
                            - log(1 - (ty > 1) * perPartner_Anal(s , ss)) ... % evaluate probability of getting STI for STI indices, i.e. if t > 1
                            .* contactProb;
                        
                        if h < 3 % getting HIV from infectious, no transmission from HIV-positive on ART/tested
                            perYear_AnalInf(h , ty , 2 , r , s , ss) = ... % probability of getting HIV
                                - log(1 - (p_hiv - joint * (ty > 1))) ...
                                .* contactProbHiv;  
                            perYear_AnalInf(h , ty , 3 , r , s , ss) = ... probability of getting HIV and STI
                                - log(1 - joint * (ty > 1)) ...
                                .* contactProbHiv;   
                            
                        end
                        
                        % oral (non-HIV STIs)
                        perYear_OralInf(h , ty , 1 , r , s , ss) = ...
                            - log(1 - (perPartner_Oral(s , ss)) * (ty > 1)) ...
                            .* contactProb;
                    end
                end
            end
        end
    end
end


%%
% Force of infection
% Rate of infection at a given site based on probability of infection by
% individuals with infections at any modelled sites.
infs = 3; % Other STI, HIV , HIV + Other STI
lambda_Anal = zeros(hivNegPosImm , risk , infs , stiTypes , sites);
lambda_Oral = lambda_Anal;
for hivStat = 1 : hivNegPosImm
    for hh = 1 : hivStatus
        hivStatPartner = 1;
        if hh > 1 && hh < 5
            hivStatPartner = 2;
        elseif hh == 5
            hivStatPartner = 3;
        end
        for r = 1 : risk
            for rr = 1 : risk
                for ty = 1 : stiTypes % (2) GC
                    for i = 1 : infs % Other STI, HIV , HIV + Other STI        
                        for s = 1 : sites
                            for ss = 1 : sites
                                lambda_Anal(hivStat , r , i , ty , s) = lambda_Anal(hivStat , r , i , ty , s) ... 
                                    + analPartnersAdj(hivStat , hivStatPartner , rr , r) ...
                                    * mixMat(hivStat , hivStatPartner , rr , r)...
                                    * perYear_AnalInf(hivStatPartner , ty , i , rr , s , ss) ...
                                    .* min((1 - condUse(r)) , (1 - condUse(rr)));
                                
                                lambda_Oral(hivStat , r , i , ty , s) = lambda_Oral(hivStat , r , i , ty , s) ...
                                    + oralPartnersAdj(hivStat , hivStatPartner , rr , r) ...
                                    * mixMat(hivStat , hivStatPartner , rr , r)...
                                    * perYear_OralInf(hivStatPartner , ty , i , rr , s , ss);
                            end
                        end
                    end
                end
            end
        end
    end
end
lambda = lambda_Anal + lambda_Oral;
if any(imag(lambda(:)))
    disp('stop')
end
%% Calculate derivatives
dPop = zeros(size(pop));
%% Transitions to infection(s)

for sTo = 1 : sites
    for r = 1 : risk 
        % HIV-negative, STI positive
        infected = lambda(1 , r , 2 , 1 , sTo) ...
            .* pop(1 , 2 : stiTypes , sTo , r);
        dPop(1 , 2 : stiTypes , sTo , r) = ...
            dPop(1 , 2 : stiTypes , sTo , r) - infected;
        dPop(2 , 2 : stiTypes , sTo , r) = ...
            dPop(2 , 2 : stiTypes , sTo , r) + infected;
        
        for tTo = 1 : stiTypes
            %STI-negative, HIV-positive / HIV-immune
            for h = 2 : hivStatus
                hivStat = 2; % infectious, tested, treated
                if h > 4 % PrEP
                    hivStat = 3;
                end
                infected = lambda(hivStat , r , 1 , tTo , sTo) ...
                    .* pop(h , 1 , 1 , r);
                dPop(h , 1 , 1 , r) = ...
                    dPop(h , 1 , 1 , r) - infected;
                dPop(h , tTo , sTo , r) = ...
                    dPop(h , tTo , sTo , r) + infected;
            end
        end
        
        for tTo = 1 : stiTypes
            %  All-negative
            for i = 1 : 3
                hTo = 1;
                if i > 1
                    hTo = 2;
                end
                infected = lambda(1 , r , i , tTo , sTo) ...
                    .* pop(1 , 1 , 1 , r);
                dPop(1 , 1 , 1 , r) = ...
                    dPop(1 , 1 , 1 , r) - infected;
                dPop(hTo , tTo , sTo , r) = ...
                    dPop(hTo , tTo , sTo , r) + infected;
            end
        end
    end
end

%% GC transitions
% % dimensions: [N, I, K, V] x [R, U, P]

% GC transition matrix
% natural recovery
gcTransMat_rec = gcClear(: , 1 : 4); % use [N, I, K, V] x [R, U, P] dimensions

% Treatment
% Uses proportion of all GC-infected by site that are treated by any 
% program/route and proportion treated by specific programs/routes to find
% proportion of all GC-infected by site that are treated by a specific 
% program/route

% symptomatic
sympTreatMat = p_symp .* kDiagTreat;

% preallocate matrices
gc_n_rec = zeros(size(pop));
gc_routTreat = zeros(size(pop));
gc_sympTreat = zeros(size(pop));
gc_psTreat = zeros(size(pop));
gc_n_routTreat = gc_routTreat;
gc_n_sympTreat = gc_sympTreat;
gc_n_psTreat = gc_psTreat;

% change in GC infected by site of infection and HIV status
for r = 1 : risk
    gc_n_rec(1 : 4 , 2 , 2 : 4 , r) = ... % natural clearance of GC
        (gcTransMat_rec .* squeeze(pop(1 : 4 , 2 , 2 : 4 , r))')';
    %HIV-negative
    gc_n_routTreat(1 , 2 , 2 : 4 , r) =  ... % GC treatment through routine screening
        routineTreatMat(: , 1) .* squeeze(pop(1 , 2 , 2 : 4 , r));
    gc_n_sympTreat(1 , 2 , 2 : 4 , r) =... % GC treatment of symptomatics
        sympTreatMat(: , 1) .* squeeze(pop(1 , 2 , 2 : 4 , r));
    gc_n_psTreat(1 , 2 , 2 : 4 , r) =... % GC treatment through partner services
        psTreatMat(: , 1) .* squeeze(pop(1 , 2 , 2 : 4 , r));
    
    % HIV-positive
    gc_routTreat(2 , 2 , 2 : 4 , r) =  ... % GC treatment and HIV testing through routine screening
        routineTreatMat(: , 2) .* squeeze(pop(2 , 2 , 2 : 4 , r));
    gc_sympTreat(2 , 2 , 2 : 4 , r) =... % GC treatment and HIV testing of symptomatics
        sympTreatMat(: , 2) .* squeeze(pop(2 , 2 , 2 : 4 , r));
    gc_psTreat(2 , 2 , 2 : 4 , r) =... % GC treatment and HIV testing through partner services
        psTreatMat(: , 2) .* squeeze(pop(2 , 2 , 2 : 4 , r));
    
end

% Natural clearance
dPop(1 : 4 , 2 , 2 : 4 , :) = dPop(1 : 4 , 2 , 2 : 4 , :) ...
    - gc_n_rec(1 : 4 , 2 , 2 : 4 , :);
dPop(1 : 4 , 1 , 1 , :) = dPop(1 : 4 , 1 , 1 , :)...
    + sum(gc_n_rec(1 : 4 , 2 , 2 : 4 , :) , 3);

% HIV transitions
% HIV negative and GC positive -> HIV negative and GC negative
dPop = dPop - gc_n_routTreat - gc_n_sympTreat - gc_n_psTreat; % HIV negative and GC positive ->
dPop(1 , 1 , 1 , :) = dPop(1 , 1 , 1 , :) + ... % -> HIV negative and GC negative
    + sum(gc_n_routTreat(1 , 2 , 2 : 4 , :) , 3) ...
    + sum(gc_n_sympTreat(1 , 2 , 2 : 4 , :) , 3)...
    + sum(gc_n_psTreat(1 , 2 , 2 : 4 , :) , 3);

% HIV infectious and GC positive -> HIV tested and GC negative
dPop = dPop - gc_routTreat - gc_sympTreat - gc_psTreat;% HIV infectious and GC positive ->
dPop(3 , 1 , 1 , :) = dPop(3 , 1 , 1 , :) ... % -> HIV tested and GC negative
    + sum(gc_routTreat(2 , 2 , 2 : 4 , :) , 3)...
    + sum(gc_sympTreat(2 , 2 , 2 : 4 , :) , 3) + ...
    + sum(gc_psTreat(2 , 2 , 2 : 4 , :) , 3);

%% HIV 
% Screening (Partner Services and Routine)
psScreen(2 , 1 , 1 , :) =... %HIV testing through partner services
    mean(psTreatMat(: , 2)) .* pop(2 , 1 , 1 , :) * 0; %no PS screening for HIV for now

routScreen(2 , 1 , 1 , :) =  ... % HIV testing through routine screening
    kHivScreen .* pop(2 , 1 , 1 , :);

% infectious -> screened
dPop(2 , 1 , 1 , :) = dPop(2 , 1 , 1 , :) - psScreen(2 , 1 , 1 , :) - routScreen(2 , 1 , 1 ,:); % infectious -> 
dPop(3 , 1 , 1 , :) = dPop(3 , 1 , 1 , :) + psScreen (2 , 1 , 1 , :) + routScreen(2 , 1 , 1 , :); % -> screened

% Treatment
% Screened -> treated
% 
hivTreat(3 , 1 , 1 , :) = kHivTreat * pop(3 , 1 , 1 , :);
gcHivTreat(3 , 2 , : , :) = kHivTreat * pop(3 , 2 , : , :);
dPop(3 , 1 , 1 , :) = dPop(3 , 1 , 1 , :) - hivTreat(3 , 1 , 1 , :); % GC-, HIV screened -> 
dPop(4 , 1 , 1 , :) = dPop(4 , 1 , 1 , :) + hivTreat(3 , 1 , 1 , :);  % -> GC-, HIV treated
dPop(3 , 2 , : , :) = dPop(3 , 2 , : , :) - gcHivTreat(3 , 2 , : , :); % GC+, HIV screened ->
dPop(4 , 1 , : , :) = dPop(4 , 1 , : , :) + gcHivTreat(3 , 2 , : , :); % -> GC-, HIV treated

%% Births and Deaths
dPop_bd = -kDie .* pop; % uniformly remove a fraction of the population
% replace fraction with a scaled-down version of the seed population
%dPop_bd = dPop_bd + popNew; 
for r = 1 : 3
    dPop_bd(1 , 1 , 1 , r) = dPop_bd(1 , 1 , 1 , r) ...
        + sum(kBorn .* pop(:)) .* riskVec(r);
end
dPop = dPop + dPop_bd;

% convert dPop back to column vector for ODE solver
dPop = dPop(:);
