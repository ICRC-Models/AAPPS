function loadUp()
hivStatus = 5; % (1) negative , (2) infectious , (3) tested , (4) treated...
% (5) on PreP, HIV immune
stiTypes = 4; % (1) none (2) gonorrhea, (3) chlamydia, (4) syphilis
sites = 4; % (1) none (2) rectal, (3) urethral, (4) pharyngeal
risk = 3;
file = 'GC_HIV_ModelParameters.xlsx';
kBorn = xlsread(file , 'Rates' , 'B29');
kDie = xlsread(file , 'Rates' , 'B30');
save('genParams');

% GC parameters
% dimensions for matrices [sites x hivStatus] (rectal, urethral,
% pharyngeal)
file = 'GC_HIV_ModelParameters.xlsx';
pSite = xlsread(file , 'Rates' , 'B4:F6');
gcTreat = xlsread(file , 'Rates' , 'B9:B11');
k_gcPS = xlsread(file , 'Rates' , 'B14:F16');
gcClear = xlsread(file , 'Rates' , 'B19:B21');
kInt = xlsread(file , 'Rates' , 'B24:F26');
perActInf = xlsread(file , 'Rates' , 'C36 : C38');
save('gcParams')

% HIV parameters
file = 'GC_HIV_ModelParameters.xlsx';

% General
rTreat_v = xlsread(file , 'Rates' , 'B44');

% GC-Susceptible
k_toPrep = xlsread(file , 'Rates' , 'B47');
k_prepOut = xlsread(file , 'Rates' , 'B48');
kTest_i = xlsread(file , 'Rates' , 'B49');
kTreat_i = xlsread(file , 'Rates' , 'B50');
kTreat_k = xlsread(file , 'Rates' , 'B51');
k_treatOut = xlsread(file , 'Rates' , 'B52');

% GC-Urethral
kprep_u = xlsread(file , 'Rates' , 'B55');
kprep_u_ps = xlsread(file , 'Rates' , 'B56');
kTest_u = xlsread(file , 'Rates' , 'B57');
kTest_u_ps = xlsread(file , 'Rates' , 'B58');
kTreat_u = xlsread(file , 'Rates' , 'B59');
kTreat_u_ps = xlsread(file , 'Rates' , 'B60');
kTreat_k = xlsread(file , 'Rates' , 'B61');
kTreat_k_ps = xlsread(file , 'Rates' , 'B62');

% GC-Pharyngeal
kprep_p =  xlsread(file , 'Rates' , 'B65');
kprep_p_ps = xlsread(file , 'Rates' , 'B66');
kTest_p = xlsread(file , 'Rates' , 'B67');
kTest_p_ps = xlsread(file , 'Rates' , 'B68');
kTreat_p = xlsread(file , 'Rates' , 'B69');
kTreat_p_ps = xlsread(file , 'Rates' , 'B70');

% GC-Rectal
kprep_r = xlsread(file , 'Rates' , 'B73');
kprep_r_ps = xlsread(file , 'Rates' , 'B74');
kTest_r = xlsread(file , 'Rates' , 'B75');
kTest_r_ps = xlsread(file , 'Rates' , 'B76');
kTreat_r = xlsread(file , 'Rates' , 'B77');
kTreat_r_ps = xlsread(file , 'Rates' , 'B78');

save('gcHivParams')

% partners per year and condom usage
file = 'GC_HIV_ModelParameters.xlsx';
xlsread(file , 'Partners' , 'A3 : B5');
c(1 , : , :) = xlsread(file , 'Partners' , 'A3 : B5');
p_cond(1 , :) = 1 - xlsread(file , 'Partners' , 'B7 : B9');
c(2 , : , :) = xlsread(file , 'Partners' , 'A13 : B15');
p_cond(2 , :) = 1 - xlsread(file , 'Partners' , 'B17 : B19');
c(3 , : , :) = xlsread(file , 'Partners' , 'A24 : B26');
p_cond(3 , :) = 1 - xlsread(file , 'Partners' , 'B28 : B30');
c(4 , : , :) = xlsread(file , 'Partners' , 'A35 : B37');
p_cond(4 , :) = 1 - xlsread(file , 'Partners' , 'B39 : B41');
c(5 , : , :) = xlsread(file , 'Partners' , 'A46 : B48');
p_cond(5 , :) = 1 - xlsread(file , 'Partners' , 'B50 : B52');
acts = 24;
save('partnerParams')

% proportion gc treat
% file = 'GC_HIV_ModelParameters.xlsx';
% xlsread(file , 'Partners' , 'A3 : B5');
% gc_prop(: , : , 1) = xlsread(file , 'Rates' , 'B57 : F59');
% gc_prop(: , : , 2) = xlsread(file , 'Rates' , 'B63 : F65');
% gc_prop(: , : , 3) = xlsread(file , 'Rates' , 'B69 : F71');
% save('gc_prop')
end
