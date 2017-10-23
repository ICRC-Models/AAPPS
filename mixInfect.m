function dPop = mixInfect(pop , perActInf , acts , partners)
%%
sumall = @(x) sum(x(:));

% perActInf = [0 , 0.84 , 0; 0.243 , 0 , 0.0865 ; 0 , 0.62 , 0]; % [U , P , R]
perAct_Anal = [0.84 , 0 , 0 ; 0 , 0.243  0 ; 0 , 0 , 0.62];
perAct_Oral = [0 , 0 , 0 ; 0 , 0.62 , 0 ; 0 , 0 , 0.62];
perActHiv = 0.82 * 10 ^ -2;

perPartner_Anal = 1 - (1 - perAct_Anal) .^ acts;
perPartner_Oral = 1 - (1 - perAct_Oral) .^ acts;
perPartnerHiv = 1 - (1 - perActHiv) .^ acts;

hAssortMat = eye(2);
rAssortMat = eye(sites);

hAssort = 0.8; % extent of HIV serosorting
rAssort = 0.7; % assorting by risk, tied to pharyngeal infection or lack thereof

mixMatHiv = zeros(2 , 2); % positive or negative
mixMatRisk = zeros(2 , risk , risk);
mixMat = zeros(2 , 2 , risk , risk);

byHiv = sum(sum(pop , 2) , 3); % sum across sti types and sites
hivNeg = sum(byHiv(1 : 2)); % negative, infectious (own status unknown)
hivPos = sum(byHiv(3 : 5)); % tested, treated, PrEP/HIV Immune
total = sumall(byHiv); % total population

% Mixing matrix by HIV status
mixMatHiv(: , 2) = hivPos / total .* (1 - hAssort); % random mixing
mixMatHiv(: , 1) = hivNeg / total .* (1 - hAssort); % random mixing
mixMatHiv = mixMatHiv + hAssort .* hAssortMat; % assortative mixing by HIV status

% Site of infection being used as a proxy for risk
byHiv_Risk = sum(pop , 2); % [hivStatus x sites]
hivNeg_Risk = sum(byHiv_Risk(1 : 2 , :) , 1) ./ hivNeg; % [1(HIV-susceptible) x sites] -> [1 x sites]
hivPos_Risk = sum(byHiv_Risk(3 : 5 , :) , 1) ./ hivPos; % sum[(HIV-positive states) 3 x sites] -> [1 x sites]

% Mixing matrix by risk group
for h = 1 : 2
    for r = 1 : sites
        mixMatRisk(1 , : , r) = hivNeg_Risk(r) .* (1 - rAssort); % pure random mixing by risk group for HIV negative
        mixMatRisk(2 , : , r) = hivPos_Risk(r) .* (1 - rAssort); % pure random mixing by risk group for HIV positive
    end
end
for h = 1 : 2
    mixMatRisk(h , : , :) = squeeze(mixMatRisk(h , : , :)) + rAssort .* rAssortMat; % assortative mixing by risk
end

% Mixing matrix by HIV status and risk group (using infection site as
% proxy)
for h = 1 : 2
    for hh = 1 : 2
        mixMat(h , hh , : , :) = mixMatHiv(h , hh) .* mixMatRisk(hh , : , :);
    end
end

% Adjust partners to account for discrepancies in reporting
% Balance according to serosorting report?
% dim(partners) = [hiv state x infection site x act type]
partners_hivPos = squeeze(sum(partners(3 : 5 , : , :) , 1)); % [infection site x act type]
partners_hivNeg = squeeze(sum(partners(1 : 2 , : , :) , 1)); % [infection site x act type]

partners_anal = [partners_hivNeg(: , 1) , partners_hivPos(: , 1)];
partners_oral = [partners_hivNeg(: , 2) , partners_hivPos(: , 2)];
fracPop = [hivNeg_Risk ; hivPos_Risk];
% fracAct = zeros(size(partners));
% actTypes = 2;
% for f = 1 : actTypes % (1) anal , (2) oral
%     fracAct(: , : , f) = partners(: , : , f) ./ sum(partners , 3);
% end
hivPosNeg = 2;
adjustFac_anal = zeros(hivPosNeg , hivPosNeg , risk , risk);
adjustFac_oral = adjustFac_anal;
% ratio of HIV+ with infection at site : HIV- with infection at site
for r = 1 : risk
    for rr = 1 : risk
        for h = 1 : hivPosNeg
            for hh = 1 : hivPosNeg
                frac = fracPop(h , r) / fracPop(hh , rr);
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
%%
analPartnersAdj = zeros(2 , 2 , sites , sites);
oralPartnersAdj = analPartnersAdj;
theta = 0.5;

for r = 1 : sites
    for rr = 1 : sites
        for h = 1 : hivPosNeg
            for hh = 1 : hivPosNeg
                % Adjust anal partners
                analPartnersAdj(h , hh , r , rr) = partners_anal(r , h) ...
                    * adjustFac_anal(h , hh , r , rr) .^ -(1 - theta);
                
                analPartnersAdj(hh , h , rr , r) = partners_anal(rr , hh) ...
                    * adjustFac_anal(h , hh , r , rr) .^ theta;
                
                % Adjust oral partners
                oralPartnersAdj(h , hh , r , rr) = partners_oral(r , h) ...
                    * adjustFac_oral(h , hh , r , rr) .^ -(1 - theta);
                
                oralPartnersAdj(hh , h , r , r) = partners_oral(rr , hh) ...
                    * adjustFac_oral(h , hh , r , rr) .^ theta;
            end
        end
    end
end
%%
% STIs (other than HIV) per year per partnership transmission. Per year per
% event rate of infection
perYear_AnalInf = zeros(5 , stiTypes - 1 , 3 , sites , sites);
perYear_OralInf = perYear_AnalInf;
for h = 1 : 5
    for s = 2 : stiTypes
        for i = 1 : sites
            for ii = 1 : sites
                popSubtotal = sum(sum(pop(h , s , :)));
                if popSubtotal > 0
                    % anal
                    joint = perPartner_Anal(ii , i) * perPartnerHiv * ...
                        heaviside(h - 3);
                    perYear_AnalInf(h , s - 1 , 1 , ii , i) = ...
                        - log(1 - (perPartner_Anal(ii , i) - joint)) ...
                        .* pop(h , s , i) ./ popSubtotal;
                    perYear_AnalInf(h , s - 1 , 2 , ii , i) = ...
                        - log(1 - perPartnerHiv * heaviside(h - 3) - joint) ...
                        .* pop(h , s , i) ./ popSubtotal;
                    perYear_AnalInf(h , s - 1 , 3 , ii , i) = ...
                        - log(1 - joint) ...
                        .* pop(h , s , i) ./ popSubtotal;
                    
                    % oral
%                    joint = perPartner_Oral(ii , i) * perPartnerHiv * ...
%                         heaviside(h - 3);
                    perYear_OralInf(h , s - 1 , 1 , ii , i) = ...
                        - log(1 - (perPartner_Oral(ii , i))) ...
                        .* pop(h , s , i) ./ popSubtotal;
%                     perYear_OralInf(h , s - 1 , 2 , ii , i) = ...
%                         - log(1 - perPartnerHiv * heaviside(h - 3) - joint) ...
%                         .* pop(h , s , i) ./ popSubtotal;
%                     perYear_OralInf(h , s - 1 , 3 , ii , i) = ...
%                         - log(1 - joint) ...
%                         .* pop(h , s , i) ./ popSubtotal;  
                end
            end
        end
    end
end


%%
% Force of infection
infs = 3; % Other STI, HIV , HIV + Other STI
lambda_Anal = zeros(2 , sites , infs , sites , stiTypes);
lambda_Oral = lambda_Anal;
for hivStat = 1 : hivPosNeg
    for hh = 1 : hivStatus
        hivStatPartner = 1;
        if hh > 2
            hivStatPartner = 2;
        end
        for s = 1 : sites
            for ss = 1 : sites
                for t = 1 : stiTypes - 1
                    for i = 1 : infs
                    lambda_Anal(hivStat , s , i , ss , t) = min(lambda_Anal(hivStat , s , i , ss , t) ...
                        + analPartnersAdj(hivStat , hivStatPartner , s , ss) ...
                        * mixMat(hivStat , hivStatPartner , s , ss)...
                        * perYear_AnalInf(hh , t , i , ss , s) , 0.99);
                    
                    lambda_Oral(hivStat , s , i , ss , t) = min(lambda_Oral(hivStat , s , i , ss , t) ...
                        + oralPartnersAdj(hivStat , hivStatPartner , s , ss) ...
                        * mixMat(hivStat , hivStatPartner , s , ss)...
                        * perYear_OralInf(hh , t , i , ss , s) , 0.99);
                    end
                end
            end
        end
    end
end             

%% HIV state transition matrix

% Matrices for HIV-state transitions for GC-Susceptible, GC-Urethral, ...
% GC-Pharyngeal , GC-Recovered.
% Dimensions: [5 x 5] , [N , I , K , V , P] x [N , I , K , V , P]
% N - negative, I - infectious , K - tested , V - HIV Treated & Virally Suppressed , ...
% P - On PrEP, HIV immune

% GC-Susceptible HIV transition matrix
gcHivTransMat(: , : , 1) = ...
    [-(lambdaHIV + k_toPrep) , 0 , 0 , 0 , k_prepOut;
    lambdaHIV , -(kTest_i + kTreat_i) , 0 , 0 , 0;
    0 , kTest_i , -kTreat_k , k_treatOut , 0;
    0 , kTreat_i , kTreat_k , k_treatOut , 0;
    k_toPrep , 0 , 0 , 0 , -k_prepOut];

% GC-Urethral HIV transition matrix
gcHivTransMat(: , : , 2) = ...
    [-(lambdaHIV + kprep_u + kprep_u_ps) , 0 , 0 , 0 , k_prepOut;
    lambdaHIV , -(kTest_u + kTest_u_ps + kTreat_u + kTreat_u_ps) , 0 , 0 , 0;
    0 , kTest_u + kTest_u_ps , -(kTreat_k + kTreat_k_ps) , rTreat_v , 0;
    0 , kTreat_u + kTreat_u_ps , kTreat_k + kTreat_k_ps , -rTreat_v , 0;
    kprep_u + kprep_u_ps , 0 ,0 , 0 , -k_prepOut];

% GC-Pharyngeal HIV transition matrix
gcHivTransMat(: , : , 3) = ...
    [-(lambdaHIV + kprep_p + kprep_p_ps) , 0 , 0 , 0 , k_prepOut;
    lambdaHIV , -(kTest_p + kTest_p_ps + kTreat_p + kTreat_p_ps) , 0 , 0 , 0;
    0 , kTest_p + kTest_p_ps , -(kTreat_p + kTreat_p_ps) , rTreat_v , 0;
    0 , kTreat_p + kTreat_p_ps , 0 , - rTreat_v , 0;
    kprep_p + kprep_p_ps , 0 , 0 , 0 , -k_prepOut];

% GC-Rectal HIV transition matrix
gcHivTransMat(: , : , 4) = ...
    [-(lambdaHIV + kprep_r + kprep_r_ps) , 0 , 0 , 0 , k_prepOut;
    lambdaHIV , -(kTest_r + kTest_r_ps + kTreat_r + kTreat_r_ps) , 0 , 0 , 0;
    0 , kTest_r + kTest_r_ps , -(kTreat_r + kTreat_r_ps) , rTreat_v , 0;
    0 , kTreat_r + kTreat_r_ps , 0 , -rTreat_v , 0;
    kprep_r + kprep_r_ps , 0 , 0 , 0 , -k_prepOut];

