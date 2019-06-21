function negSumLogL = likeFun(popVec , newHiv , ...
    hivPrev_obs , hivInc_obs , gcRecPrev_obs , gcUrePrev_obs , gcPhaPrev_obs , ...
    hivStatus , stiTypes , sites , risk , stepsPerYear , startYear)

%% HIV prevalence
hivYearVec = unique(hivPrev_obs(:,1));
hivPrev = zeros(1 , length(hivYearVec));
for t = 1 : length(hivYearVec)
    time = (hivYearVec(t) - startYear) * stepsPerYear;
    totalPop = sum(sum(sum(sum(popVec(time , : , : , : , :) , 2) , 3) , 4) , 5);
    hivPop = sum(sum(sum(sum(popVec(time , 2 : 4 , : , : , :), 2) , 3) , 4) , 5); 
    hivPrev(1,t) = (hivPop / totalPop) * 100;
end
pPos = [hivPrev(:)];
N = [hivPrev_obs(:,3)];
nPos = [hivPrev_obs(:,2)];

%% HIV incidence
% need to look up how to construct likelihood function for incidence
% assuming a normal distribution

%% GC prevalence
gcYearVec = unique(gcRecPrev_obs(:,1));

% Rectal GC
% gcRecPrev = zeros(1 , length(gcYearVec));
% for t = 1 : length(gcYearVec)
%     time = (gcYearVec(t) - startYear) * stepsPerYear
%     totalPop = sum(sum(sum(sum(popVec(time , : , : , : , :) , 2) , 3) , 4) , 5);
%     gcRecPop = sum(sum(popVec(time , : , 2 , 2 , :) , 2) , 5); 
%     gcRecPrev(1,t) = (gcRecPop / totalPop) * 100;
% end
% pPos = [pPos; gcRecPrev(:)];
% N = [N ; gcRecPrev_obs(:,3)];
% nPos = [nPos ; gcRecPrev_obs(:,2)];

% Urethral GC
gcUrePrev = zeros(1 , length(gcYearVec));
for t = 1 : length(gcYearVec)
    time = (gcYearVec(t) - startYear) * stepsPerYear;
    totalPop = sum(sum(sum(sum(popVec(time , : , : , : , :) , 2) , 3) , 4) , 5);
    gcUrePop = sum(sum(popVec(time , : , 2 , 3 , :) , 2) , 5); 
    gcUrePrev(1,t) = (gcUrePop / totalPop) * 100;
end
pPos = [pPos; gcUrePrev(:)];
N = [N ; gcUrePrev_obs(:,3)];
nPos = [nPos ; gcUrePrev_obs(:,2)];

% Pharyngeal GC
% gcPhaPrev = zeros(1 , length(gcYearVec));
% for t = 1 : length(gcYearVec)
%     time = (gcYearVec(t) - startYear) * stepsPerYear
%     totalPop = sum(sum(sum(sum(popVec(time , : , : , : , :) , 2) , 3) , 4) , 5);
%     gcPhaPop = sum(sum(popVec(time , : , 2 , 4 , :) , 2) , 5); 
%     gcPhaPrev(1,t) = (gcPhaPop / totalPop) * 100;
% end
% pPos = [pPos; gcPhaPrev(:)];
% N = [N ; gcPhaPrev_obs(:,3)];
% nPos = [nPos ; gcPhaPrev_obs(:,2)];

%% Likelihood function
pPos = pPos ./ 100; % scale percent probabilities to decimals
logL = nPos .* log(pPos) + (N - nPos) .* log(1 - pPos); % log likelihoods for binomial events
negSumLogL = - sum(logL); % negative logL to be minimized
