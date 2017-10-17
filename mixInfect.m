function dPop = mixInfect(pop , perActInf , acts , partners)
%%
sumall = @(x) sum(x(:));

perActInf = [0 , 0.84 , 0; 0.243 , 0 , 0.0865 ; 0 , 0.62 , 0]; % [U , P , R]
perActHiv = 0.82 * 10 ^ -2;

perPartnerInf = min(1 - (1 - perActInf) .^ acts , 0.99);
perPartnerHiv = min(1 - (1 - perActHiv) .^ acts , 0.99);

hAssortMat = eye(2);
rAssortMat = eye(sites);

hAssort = 0.8; % extent of HIV serosorting
rAssort = 0.7; % assorting by risk, tied to pharyngeal infection or lack thereof

mixMatHiv = zeros(2 , 2); % positive or negative
mixMatRisk = zeros(2 , sites , sites);
mixMat = zeros(2 , 2 , sites , sites);

byHiv = sum(sum(pop , 2) , 3); % sum across sti types and sites
hivNeg = byHiv(1) + byHiv(5); % negative, on PrEP, HIV immune
hivPos = sum(byHiv(2 : 4)); % infectious, tested, treated
total = sumall(byHiv); % total population

% Mixing matrix by HIV status
mixMatHiv(: , 1) = hivPos / total .* (1 - hAssort); % random mixing
mixMatHiv(: , 2) = hivNeg / total .* (1 - hAssort); % random mixing
mixMatHiv = mixMatHiv + hAssort .* hAssortMat; % assortative mixing by HIV status

% Site of infection being used as a proxy for risk
byHiv_Risk = sum(pop , 2); % [hivStatus x sites]
hivNeg_Risk = (byHiv_Risk(1 , :) + byHiv_Risk(5 , :)) ./ hivNeg; % [1(HIV-susceptible) x sites] -> [1 x sites]
hivPos_Risk = sum(byHiv_Risk(2 : 4 , :) , 1) ./ hivPos; % sum[(HIV-positive states) 3 x sites] -> [1 x sites]

% Mixing matrix by risk group
for h = 1 : 2
    for r = 1 : sites
        mixMatRisk(h , : , r) = hivNeg_Risk(r) .* (1 - rAssort); % pure random mixing by site for non-pharyngeal / pharyngeal infected HIV negative / HIV positive
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
partners_hivPos = squeeze(sum(partners(2 : 4 , : , :) , 1)); % [infection site x act type]
partners_hivNeg = squeeze(partners(1 , : , :) + partners(5 , : , :)); % [infection site x act type]

partners_anal = [partners_hivNeg(: , 1), partners_hivPos(: ,1)];
fracPop = [hivNeg_Risk ; hivPos_Risk];
% fracArt = zeros(size(partners));
% actTypes = 2;
% for f = 1 : actTypes % (1) anal , (2) oral
%     fracAct(: , : , f) = partners(: , : , f) ./ sum(partners , 3);
% end
hivPosNeg = 2;
adjustFac = zeros(hivPosNeg , hivPosNeg , sites , sites);
% ratio of HIV+ with infection at site : HIV- with infection at site
for s = 1 : sites
    for ss = 1 : sites
        for h = 1 : hivPosNeg
            for hh = 1 : hivPosNeg
                frac = fracPop(h , s) / fracPop(hh , ss);
                adjustFac(h , hh , s , ss) = ...
                    sum(partners_anal(s , h) * mixMat(h , hh , s ,ss)) ...
                    / (sum(partners_anal(ss , hh)) * mixMat(hh , h , ss , s)) * frac;
            end
        end
    end
end
%%
partnersAdj = zeros(2 , 2 , sites , sites);
theta = 0.5;

for s = 1 : sites
    for ss = 1 : sites
        for h = 1 : hivPosNeg
            for hh = 1 : hivPosNeg
                partnersAdj(h , hh , s , ss) = partners_anal(s , h) ...
                    * adjustFac(h , hh , s , ss) .^ -(1 - theta);

                partnersAdj(hh , h , ss , s) = partners_anal(ss , hh) ...
                    * adjustFac(h , hh , s , ss) .^ theta;
            end
        end
    end
end
%%
% STIs (other than HIV) per year per partnership transmission
perYearInf = zeros(5 , stiTypes - 1 , 3 , sites , sites);
for h = 1 : 5
    for s = 2 : stiTypes
        for i = 1 : sites
            for ii = 1 : sites
                popSubtotal = sum(sum(pop(h , s , :)));
                if popSubtotal > 0
                    joint = perPartnerInf(ii , i) * perPartnerHiv * ...
                        heaviside(h - 3);
                    perYearInf(h , s - 1 , 1 , ii , i) = ...
                        - log(1 - (perPartnerInf(ii , i) - joint)) ...
                        .* pop(h , s , i) ./ popSubtotal;
                    perYearInf(h , s - 1 , 2 , ii , i) = ...
                        - log(1 - perPartnerHiv * heaviside(h - 3) - joint) ...
                        .* pop(h , s , i) ./ popSubtotal;
                    perYearInf(h , s - 1 , 3 , ii , i) = ...
                        - log(1 - perPartnerInf(ii , i) * perPartnerHiv...
                        * heaviside(h - 3)) ...
                        .* pop(h , s , i) ./ popSubtotal;
                end
            end
        end
    end
end


%%
% STIs (other than HIV) lambda
infs = 3; % Other STI, HIV , HIV + Other STI
lambda = zeros(2 , sites , infs , sites , stiTypes);
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
                    lambda(hivStat , s , i , ss , t) = lambda(hivStat , s , i , ss , t) ...
                        + partnersAdj(hivStat , hivStatPartner , s , ss) ...
                        * mixMat(hivStat , hivStatPartner , s , ss)...
                        * perYearInf(hh , t , i , ss , s);
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

