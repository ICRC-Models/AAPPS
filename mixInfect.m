function dPop = mixInfect(pop , perActInf , acts , partners)
sumall = @(x) sum(x(:));

perPartnerInf = 1 - (1 - perActInf) .^ acts;

hAssortMat = eye(2);
rAssortMat = eye(3);

hAssort = 0.8; % extent of HIV serosorting
rAssort = 0.7; % assorting by risk, tied to pharyngeal infection or lack thereof

mixMatHiv = zeros(2 , 2); % positive or negative
mixMatRisk = zeros(2 , risk , risk);
mixMat = zeros(2 , 2 , risk , risk);

byHiv = sum(sum(pop , 2) , 3);
hivNeg = byHiv(1) + byHiv(5);
hivPos = sum(byHiv(2 : 4));
total = sumall(byHiv);

% Mixing matrix by HIV status
mixMatHiv(1 , :) = hivPos / total .* (1 - hAssort); % random mixing
mixMatHiv(2 , :) = hivNeg / total .* (1 - hAssort); % random mixing
mixMatHiv = mixMatHiv + hAssort .* hAssortMat; % assortative mixing by HIV status

byHiv_Risk = sum(pop , 2); % [hivStatus x sites]
hivNeg_Risk = byHiv_Risk(1 , :) + byHiv_Risk(5 , :); % [1(HIV-susceptible) x sites] -> [1 x sites]
hivPos_Risk = sum(byHiv_Risk(2 : 4 , :) , 1); % sum[(HIV-positive states) 3 x sites] -> [1 x sites]

% Mixing matrix by risk group
mixMatRisk(1 , 1 : 2 , :) = hivNeg_Risk; % pure random mixing by site for non-pharyngeal infected HIV negative
mixMatRisk(1 , 3 , :) = hivNeg_Risk .* (1 - rAssort); % random mixing by risk among HIV negative with pharyngeal infection
mixMatRisk(2 , 1 : 2 , :) = hivPos_Risk; % pure random mixing by site for non-pharyngeal infected HIV positive
mixMatRisk(2 , 3 , :) = hivPos_Risk .* (1 - rAssort); % random mixing by risk among HIV positive with pharyngeal infection
mixMatRisk = mixMatRisk + rAssort .* rAssortMat; % assortative mixing by risk

% Mixing matrix by HIV status and risk group
for h = 1 : 2
    for hh = 1 : 2
        mixMat(hh , h , : , :) = mixMatHiv(hh , h) .* mixMatRisk(h , : , :);
    end
end

% Adjust partners to account for discrepancies in reporting
% Balance according to serosorting report?
adjustFac = zeros(risk , risk);
for s = 1 : sites
    for ss = 1 : sites
        adjustFac(s ,ss) = ...
            sum(partners(1 , 1 , s)) * rho(1 , 2 , s ,ss) * hivNeg ...
            / (sum(partners(2 , 2 : hivStatus , ss)) * rho(2 , 1 , ss , s) * hivPos);
    end
end

partnersAdj = zeros(2 , 2 , sites , sites);
theta = 0.5;
for site = 1  : sites
    partnersAdj(1 , : , site , :) = partners(1 , :) .* adjustFac(site , :) .^ -(1 - theta);
    partnersAdj(2 , : , : , site) = partners(2 , :) .* adjustFac(site , :) .^ theta;
end


perYearInf = zeros(risk , risk);
for h = 1 : hivStatus
  for hh = 1 : hivStatus
    for s = 1 : stiTypes
      for i = 1 : sites
        popSubtotal = sum(pop(h , s , :));
        perYearInf(hh , h , s , i) = perYearInf(hh , h , s , i)...
          -log(1 - perPartnerInf(h , s , i)) .* pop(h , s , i) ./ popSubtotal;
        end
    end
  end
end

% Calculate lambda for each STI
lambda = zeros(2 , stiTypes , sites);
  for site = 1 : sites
    for h = 1 : 2
      for hh = 1 : 2
        for s = 1 : stiTypes
          lambda(h , site , s) = lambda(h , site , s) + partnersAdj(h , hh , : , site)...
            * mixMat(h , hh , : , site) * perYearInf(h , hh ,  site , s);
        end
      end
    end
  end
