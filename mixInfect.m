function dPop = mixInfect(t , pop , hivStatus , stiTypes , sites , ...
    risk , kBorn , kDie , pSite , gcTreat , k_gcPS , gcClear , kInt , ...
    perActInf , rTreat_v , k_toPrep, k_prepOut , kTest_i , ...
    kTreat_i , k_treatOut , kprep_u , kprep_u_ps , kTest_u , kTest_u_ps , ...
    kTreat_u , kTreat_u_ps , kTreat_k , kTreat_k_ps , kprep_p , kprep_p_ps , kTest_p , ...
    kTest_p_ps , kTreat_p , kTreat_p_ps , kprep_r , kprep_r_ps , kTest_r , ...
    kTest_r_ps , kTreat_r , kTreat_r_ps , partners , p_cond , acts , riskVec , ...
    condUse , tVec)
%%
sumall = @(x) sum(x(:));

%% Scale-up vectors
kTest_p_ps = interp1(tVec , kTest_p_ps , t);
kTreat_p_ps = interp1(tVec , kTreat_p_ps , t);
kTest_u_ps = interp1(tVec , kTest_u_ps , t);
kTreat_u_ps = interp1(tVec , kTreat_u_ps , t);
kTest_r_ps = interp1(tVec , kTest_r_ps , t);
kTreat_r_ps = interp1(tVec , kTreat_r_ps , t);
k_gcPS(1 : 3) = interp1(tVec , k_gcPS , t);

%%
pop = reshape(pop , [hivStatus , stiTypes , sites , risk]);

% perActInf = [0 , 0.84 , 0; 0.243 , 0 , 0.0865 ; 0 , 0.62 , 0];
perAct_Anal = [0 , 0.84 , 0 ; 0.243 , 0  0.0865 ; 0 , 0.62 , 0];
perAct_Oral = [0 , 0 , 0 ; 0 , 0.62 , 0.62 ; 0 , 0.62 , 0.62];
perActHiv = 0.82 * 10 ^ -2;

acts = acts / 2; % half anal, half oral for testing

perPartner_Anal = min(1 - (1 - perAct_Anal) .^ acts , 0.99); % GC anal transmission probability
perPartner_Oral = min(1 - (1 - perAct_Oral) .^ acts , 0.99); % GC oral transmission probability
perPartnerHiv = min(1 - (1 - perActHiv) .^ acts , 0.99);

hAssortMat = eye(3);% (1) HIV-/unknown status, (2) HIV+/Tested, (3) PrEP
rAssortMat = eye(risk);

hAssort = 0.8;%0.8; % extent of HIV serosorting
rAssort = 0.5;%0.5; % assorting by risk

mixMatHiv = zeros(3 , 3); % positive, negative, PrEP
mixMatRisk = zeros(3 , risk , risk);
mixMat = zeros(3 , 3 , risk , risk);

byHiv = squeeze(sum(sum(pop , 2) , 3)); % sum across sti types and sites
hivNeg = sum(sum(byHiv(1 : 2 , :))); % negative, infectious (own status unknown)
hivPos = sum(sum(byHiv(3 : 4 , :))); % tested, treated
hivImm = sum(sum(byHiv(5 , :))); %PrEP/HIV Immune
total = sumall(byHiv); % total population

% Mixing matrix by HIV status
mixMatHiv(: , 3) = hivImm / total .* (1 - hAssort); % random mixing with PrEP
mixMatHiv(: , 2) = hivPos / total .* (1 - hAssort); % random mixing with HIV positive
mixMatHiv(: , 1) = hivNeg / total .* (1 - hAssort); % random mixing with HIV-/status unknown
mixMatHiv = mixMatHiv + hAssort .* hAssortMat; % random + assortative mixing by HIV status

%
byHiv_Risk = squeeze(sum(sum(pop , 2) , 3)); % [hivStatus x sites]
hivNeg_Risk = sum(byHiv_Risk(1 : 2 , :) , 1) ./ hivNeg; % [1(HIV-susceptible) x sites] -> [1 x sites]
hivPos_Risk = sum(byHiv_Risk(3 : 4 , :) , 1) ./ hivPos; % sum[(HIV-positive states) 2 x sites] -> [1 x sites]
hivImm_Risk = sum(byHiv_Risk(5 , :) , 1) ./ hivImm; % sum[(HIV-immune) 1 x sites] -> [1 x sites]

% Mixing matrix by risk group
for r = 1 : risk
    mixMatRisk(1 , : , r) = hivNeg_Risk(r) .* (1 - rAssort); % pure random mixing by risk group for HIV negative
    mixMatRisk(2 , : , r) = hivPos_Risk(r) .* (1 - rAssort); % pure random mixing by risk group for HIV positive
    mixMatRisk(3 , : , r) = hivImm_Risk(r) .* (1 - rAssort); % pure random mixing by risk group for HIV immune
end

for h = 1 : 3
    mixMatRisk(h , : , :) = squeeze(mixMatRisk(h , : , :)) + rAssort .* rAssortMat; % assortative mixing by risk
end

% Mixing matrix by HIV status and risk group (using infection site as
% proxy)
for h = 1 : 3
    for hh = 1 : 3
        mixMat(h , hh , : , :) = mixMatHiv(h , hh) .* mixMatRisk(hh , : , :);
    end
end

% Adjust partners to account for discrepancies in reporting
% Balance according to serosorting report?
% dim(partners) = [hiv state x risk x act type]
partners_hivImm = squeeze(mean(partners(5 , : , :) , 1)); % [risk x act type]
partners_hivPos = squeeze(mean(partners(3 : 4 , : , :) , 1)); % [risk x act type]
partners_hivNeg = squeeze(mean(partners(1 : 2 , : , :) , 1)); % [risk x act type]

partners_anal = [partners_hivNeg(: , 1) , partners_hivPos(: , 1) , partners_hivImm(: , 1)];
partners_oral = [partners_hivNeg(: , 2) , partners_hivPos(: , 2) , partners_hivImm(: , 2)];
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
for r = 1 : risk
    for rr = 1 : risk
        for h = 1 : hivNegPosImm
            for hh = 1 : hivNegPosImm
                frac = 0;
                if fracPop(hh , rr) > 10 ^ -6
                    frac = fracPop(h , r) / fracPop(hh , rr);
                end
                if(sum(partners_anal(rr , hh)) * mixMat(hh , h , rr , r)) > 10 ^ -6
                    adjustFac_anal(h , hh , r , rr) = ...
                        sum(partners_anal(r , h) * mixMat(h , hh , r ,rr)) ...
                        / (sum(partners_anal(rr , hh)) * mixMat(hh , h , rr , r)) * frac;
                    adjustFac_oral(h , hh , r , rr) = ...
                        sum(partners_oral(r , h) * mixMat(h , hh , r ,rr)) ...
                        / (sum(partners_oral(rr , hh)) * mixMat(hh , h , rr , r)) * frac;
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
                % Adjust anal partners
                if adjustFac_anal(h , hh , r , rr) ~= 0
                    analPartnersAdj(h , hh , r , rr) = partners_anal(r , h) ...
                        * adjustFac_anal(h , hh , r , rr) .^ -(1 - theta);
               
                    analPartnersAdj(hh , h , rr , r) = partners_anal(rr , hh) ...
                        * adjustFac_anal(h , hh , r , rr) .^ theta;
                end
                % Adjust oral partners
                if adjustFac_oral(h , hh , r , rr) ~= 0
                    oralPartnersAdj(h , hh , r , rr) = partners_oral(r , h) ...
                        * adjustFac_oral(h , hh , r , rr) .^ -(1 - theta);

                    oralPartnersAdj(hh , h , rr , r) = partners_oral(rr , hh) ...
                        * adjustFac_oral(h , hh , r , rr) .^ theta;
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
popSubs = [hivNeg ; hivPos ; hivImm]; % pop subtotals by hiv status. hivNeg includes HIV+ who think they are negative
hivNegPop = sum(pop(1 : 2 , : , : , :) , 1);
hivPosActual = sum(pop(2 , : , : , :) , 1);
hivPosPop = sum(pop(3 : 4 , : , : , :) , 1);
hivImmPop = pop(5 , : , : , :);
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
                    popSubtotal = popSubs(h);
                    popGroup = squeeze(pops{h});
                    if popSubtotal > 10 ^ -6 
                        % if h ~= 1 (not hiv negative group) 
                        contactProb = 0;
                        if popGroup(ty , s , r) > 10 ^ -6
                            contactProb = popGroup(ty , s , r) ./ popSubtotal; % proportion of pop with sti 
                        end
                        contactProbHiv = 0; % for HIV positive or HIV-immune group (HIV status known)
                        % per partner transmission probabilities
                        joint = perPartner_Anal(ss , s) * perPartnerHiv; % sti and hiv
                        p_hiv = perPartnerHiv; % hiv
                                                     
                        if h == 1 && hivPosActual(1 , ty , s , r) > 10 ^ -6
                            contactProbHiv = hivPosActual(1 , ty , s , r) ./ popSubtotal; % proportion of "hiv-neg" that are really hiv-pos
                        end
                        
                        % ADD protection from condoms here
                        perYear_AnalInf(h , ty , 1 , r , s , ss) = ... % probability of getting STI
                            - log(1 - (ty > 1) * (perPartner_Anal(ss , s) - joint)) ... % evaluate probability of getting STI for STI indices, i.e. if t > 1
                            .* (1 - condUse(r)) .* contactProb;
                        
                        if h < 3 % getting HIV from non-immune
                            perYear_AnalInf(h , ty , 2 , r , s , ss) = ... % probability of getting HIV
                                - log(1 - (p_hiv - joint * (ty > 1))) ...                        
                                .* (1 - condUse(r)) .* contactProbHiv;  % ADD protection from condoms here
                            perYear_AnalInf(h , ty , 3 , r , s , ss) = ... probability of getting HIV and STI
                                - log(1 - joint * (ty > 1)) ...
                                .* (1 - condUse(r)) .* contactProbHiv;   % ADD protection from condoms here

                        end
                        
                        % oral (non-HIV STIs)
                        perYear_OralInf(h , ty , 1 , r , s , ss) = ...
                            - log(1 - (perPartner_Oral(ss , s)) * (ty > 1)) ...
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
                                    + analPartnersAdj(hivStat , hivStatPartner , r , rr) ...
                                    * mixMat(hivStat , hivStatPartner , r , rr)...
                                    * perYear_AnalInf(hivStatPartner , ty , i , rr , s , ss);
                                
                                lambda_Oral(hivStat , r , i , ty , s) = lambda_Oral(hivStat , r , i , ty , s) ...
                                    + oralPartnersAdj(hivStat , hivStatPartner , r , rr) ...
                                    * mixMat(hivStat , hivStatPartner , r , rr)...
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
% lambda(2 : 3 , : , 2 : 3 , : , :) = 0; % HIV-positive/immune cannot be reinfected   
% dims: (hivNegPosImm (3) , risk , infs (Other STI, HIV , HIV + Other STI) , stiTypes , sites)
% lambda = min(lambda , 0.999);
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
        
        for tTo = 2 : stiTypes
            %STI-negative, HIV-positive / HIV-immune
            for h = 2 : hivStatus
                hivStat = 2;
                if h > 3
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
% dimensions: [S, U, P, R] x [S, U, P, R]

% GC transition matrix
% Infection
% gcTransMat_Inf = [-(lambda_u + lambda_p + lambda_r) , ...
%     0 , 0 , 0 ; ...
%     lambda_u , 0 , 0 , 0 ; ...
%     lambda_p , 0 , 0 , 0; ...
%     lambda_r , 0 , 0 , 0];

% recovery
gcTransMat_rec = [0 , gcClear(1), gcClear(2) , gcClear(3) ; ...
    0 , -gcClear(1) , 0 , 0 ; ...
    0 , 0 , -gcClear(2) , 0; ...
    0 , 0 , 0 , -gcClear(3)];

% treatment
gcTransMat_trt = [0 , gcTreat(1) , gcTreat(2) , gcTreat(3); ...
    0 , -gcTreat(1) , 0 , 0 ; ...
    0 , 0 , -gcTreat(2) 0; ...
    0 , 0 , 0 , -gcTreat(3)];

% partner services
gcTransMat_ps = [0 , k_gcPS(1) , k_gcPS(2) , k_gcPS(3); ...
    0 , -k_gcPS(1) , 0 , 0 ; ...
    0 , 0 , -k_gcPS(2) , 0; ...
    0 , 0 , 0 , -k_gcPS(3)];


dPop_gc_rec = zeros(size(pop));
dPop_gc_trt = zeros(size(pop));
dPop_gc_ps = zeros(size(pop));
for h = 1 : hivStatus
    for r = 1 : risk
        dPop_gc_rec(h , 2 , : , r) = dPop_gc_rec(h , 2 , : , r) ...
            + reshape(gcTransMat_rec * squeeze(pop(h , 2 , : , r)) , [1 1 4]);
        dPop_gc_trt(h , 2 , : , r) = dPop_gc_trt(h , 2 , : , r) ...
            + reshape(gcTransMat_trt * squeeze(pop(h , 2 , : , r)) , [1 1 4]);
        dPop_gc_ps(h , 2 , : , r) = dPop_gc_ps(h , 2 , : , r) ...
            + reshape(gcTransMat_ps * squeeze(pop(h , 2 , : , r)) , [1 1 4]);
    end
end

dPop = dPop + dPop_gc_rec + dPop_gc_trt + dPop_gc_ps;
dPop(: , 1 , 1 , :) = dPop(: , 1 , 1 , :) + dPop(: , 2 , 1 , :); % put GC recovered back into susceptible pool
dPop(: , 2 , 1 , :) = 0;
% HIV state transition matrix
% 
% Matrices for HIV-state transitions for GC-Susceptible, GC-Urethral, ...
% GC-Pharyngeal , GC-Recovered.
% Dimensions: [5 x 5] , [N , I , K , V , P] x [N , I , K , V , P]
% N - negative(1), I - infectious(2) , K - tested(3) , V - HIV Treated & Virally Suppressed(4) , ...
% P - On PrEP, HIV immune(5)
% 
% GC-Susceptible HIV transition matrix
gcHivTransMat(: , : , 1) = ...
    [-(k_toPrep) , 0 , 0 , 0 , k_prepOut;
    0 , -(kTest_i) , 0 , 0 , 0;
    0 , kTest_i , -kTreat_k , k_treatOut , 0;
    0 , 0 , kTreat_k , -k_treatOut , 0;
    k_toPrep , 0 , 0 , 0 , -k_prepOut];
% 
% GC-Rectal HIV transition matrix
gcHivTransMat(: , : , 2) = ...
    [-(kprep_r + kprep_r_ps) , 0 , 0 , 0 , k_prepOut;
    0 , -(kTest_r + kTest_r_ps) , 0 , 0 , 0;
    0 , kTest_r + kTest_r_ps , -(kTreat_r + kTreat_r_ps) , rTreat_v , 0;
    0 , 0 , kTreat_r + kTreat_r_ps , -rTreat_v , 0;
    kprep_r + kprep_r_ps , 0 , 0 , 0 , -k_prepOut];

% 
% GC-Urethral HIV transition matrix
gcHivTransMat(: , : , 3) = ...
    [-(kprep_u + kprep_u_ps) , 0 , 0 , 0 , k_prepOut;
    0 , -(kTest_u + kTest_u_ps) , 0 , 0 , 0;
    0 , kTest_u + kTest_u_ps , -(kTreat_u + kTreat_u_ps) , rTreat_v , 0;
    0 , 0 , kTreat_u + kTreat_u_ps , -rTreat_v , 0;
    kprep_u + kprep_u_ps , 0 ,0 , 0 , -k_prepOut];
% 
% GC-Pharyngeal HIV transition matrix
gcHivTransMat(: , : , 4) = ...
    [-(kprep_p + kprep_p_ps) , 0 , 0 , 0 , k_prepOut;
    0 , -(kTest_p + kTest_p_ps) , 0 , 0 , 0;
    0 , kTest_p + kTest_p_ps , -(kTreat_p + kTreat_p_ps) , rTreat_v , 0;
    0 , 0 , kTreat_p + kTreat_p_ps , - rTreat_v , 0;
    kprep_p + kprep_p_ps , 0 , 0 , 0 , -k_prepOut];


dPop_Hiv = zeros(size(pop));

for s = 1 : sites
    for r = 1 : risk
        % update HIV status for HIV and GC coinfected
        dPop_Hiv(: , 2 , s , r) = dPop_Hiv(: , 2 , s , r) + ...
            gcHivTransMat(: , : , s) ...
            * squeeze(pop(: , 2 , s , r));
        % update HIV status for HIV infected
        dPop_Hiv(: , 1 , s , r) = dPop_Hiv(: , 1 , s , r) + ...
            gcHivTransMat(: , : , s) ...
            * squeeze(pop(: , 1 , s , r));
    end
end

dPop = dPop + dPop_Hiv;

%% Births and Deaths
dPop_bd = -kDie .* pop;
for r = 1 : 3
    dPop_bd(1 , 1 , 1 , r) = dPop_bd(1 , 1 , 1 , r) ...
        + sum(kBorn .* pop(:)) .* riskVec(r);
end
dPop = dPop + dPop_bd;

% convert dPop back to column vector for ODE solver
dPop = dPop(:);
