function loadUp(stepSize)
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
pSite = xlsread(file , 'Rates' , 'B4:F6') * stepSize;
gcTreat = xlsread(file , 'Rates' , 'B9:F11') * stepSize;
k_gcPS = xlsread(file , 'Rates' , 'B14:F16') * stepSize;
gcClear = xlsread(file , 'Rates' , 'B19:F21') * stepSize;
kInt = xlsread(file , 'Rates' , 'B24:F26') * stepSize;
perActInf = xlsread(file , 'Rates' , 'C36 : C38') * stepSize;
save('gcParams')

% HIV parameters
file = 'GC_HIV_ModelParameters.xlsx';

% General
rTreat_v = xlsread(file , 'Rates' , 'B44') * stepSize;

% GC-Susceptible
k_toPrep = xlsread(file , 'Rates' , 'B47') * stepSize;
k_prepOut = xlsread(file , 'Rates' , 'B48') * stepSize;
kTest_i = xlsread(file , 'Rates' , 'B49') * stepSize;
kTreat_i = xlsread(file , 'Rates' , 'B50') * stepSize;
kTreat_k = xlsread(file , 'Rates' , 'B51') * stepSize;
k_treatOut = xlsread(file , 'Rates' , 'B52') * stepSize;

% GC-Urethral
kprep_u = xlsread(file , 'Rates' , 'B55') * stepSize;
kprep_u_ps = xlsread(file , 'Rates' , 'B56') * stepSize;
kTest_u = xlsread(file , 'Rates' , 'B57') * stepSize;
kTest_u_ps = xlsread(file , 'Rates' , 'B58') * stepSize;
kTreat_u = xlsread(file , 'Rates' , 'B59') * stepSize;
kTreat_u_ps = xlsread(file , 'Rates' , 'B60') * stepSize;
kTreat_k = xlsread(file , 'Rates' , 'B61') * stepSize;
kTreat_k_ps = xlsread(file , 'Rates' , 'B62') * stepSize;

% GC-Pharyngeal
kprep_p =  xlsread(file , 'Rates' , 'B65') * stepSize;
kprep_p_ps = xlsread(file , 'Rates' , 'B66') * stepSize;
kTest_p = xlsread(file , 'Rates' , 'B67') * stepSize;
kTest_p_ps = xlsread(file , 'Rates' , 'B68') * stepSize;
kTreat_p = xlsread(file , 'Rates' , 'B69') * stepSize;
kTreat_p_ps = xlsread(file , 'Rates' , 'B70') * stepSize;

% GC-Rectal
kprep_r = xlsread(file , 'Rates' , 'B73') * stepSize;
kprep_r_ps = xlsread(file , 'Rates' , 'B74') * stepSize;
kTest_r = xlsread(file , 'Rates' , 'B75') * stepSize;
kTest_r_ps = xlsread(file , 'Rates' , 'B76') * stepSize;
kTreat_r = xlsread(file , 'Rates' , 'B77') * stepSize;
kTreat_r_ps = xlsread(file , 'Rates' , 'B78') * stepSize;

save('gcHivParams')

% partners per year and condom usage
file = 'GC_HIV_ModelParameters.xlsx';
xlsread(file , 'Partners' , 'A3 : B5');
c(1 , : , :) = xlsread(file , 'Partners' , 'A3 : B5') * stepSize;
p_cond(1 , :) = 1 - xlsread(file , 'Partners' , 'B7 : B9');
c(2 , : , :) = xlsread(file , 'Partners' , 'A13 : B15') * stepSize;
p_cond(2 , :) = 1 - xlsread(file , 'Partners' , 'B17 : B19');
c(3 , : , :) = xlsread(file , 'Partners' , 'A24 : B26') * stepSize;
p_cond(3 , :) = 1 - xlsread(file , 'Partners' , 'B28 : B30');
c(4 , : , :) = xlsread(file , 'Partners' , 'A35 : B37') * stepSize;
p_cond(4 , :) = 1 - xlsread(file , 'Partners' , 'B39 : B41');
c(5 , : , :) = xlsread(file , 'Partners' , 'A46 : B48') * stepSize;
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
