% Makes matrices for HIV-state transitions for GC-Susceptible, GC-Urethral, ...
% GC-Pharyngeal , GC-Recovered.
% Dimensions: [5 x 5] , [N , I , K , V , P] x [N , I , K , V , P]
% N - negative, I - infectious , K - tested , V - HIV Treated & Virally Suppressed , ...
% P - On PrEP, HIV immune

function makeMats()

% GC-Susceptible HIV transition matrix
gcHivTransMat(: , : , 1) = ...
    [-(lambdaHIV + k_toPrep) , 0 , 0 , 0 , k_prepOut;
    lambdaHIV , -(kTest_i + kTreat_i) , 0 , 0 , 0;
    0 , kTest_i , -kTreat_k , k_treatOut , 0;
    0 , kTreat_i , kTreat_k , k_treatOut , 0;
    k_toPrep , 0 , 0 , 0 , -k_prepOut];

% GC-Urethral HIV transition matrix
gcHivTransMat(: , : , 2) = ...
    [-(lambdaHIV + kprep_u + kprep_u_ps) , 0 , 0 ,0 , k_prepOut;
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

save('gcHivTransMat')
