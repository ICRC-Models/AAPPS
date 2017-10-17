function loadUp()
hivStatus = 5; % (1) negative , (2) infectious , (3) tested , (4) treated...
% (5) on PreP, HIV immune
stiTypes = 4; % (1) none (2) gonorrhea, (3) chlamydia, (4) syphilis
sites = 3; % (1) rectal, (2) urethral, (3) pharyngeal
risk = 3;
file = 'GC_HIV_ModelParameters.xlsx';
kBorn = xlsread(file , 'Rates' , 'B29');
kDie = xlsread(file , 'Rates' , 'B30');
save('genParams');
clear

% dimensions for matrices [sites x hivStatus]
file = 'GC_HIV_ModelParameters.xlsx';
pSite = xlsread(file , 'Rates' , 'B4:F6');
gcTreat = xlsread(file , 'Rates' , 'B9:F11');
k_gcPS = xlsread(file , 'Rates' , 'B14:F16');
gcClear = xlsread(file , 'Rates' , 'B19:F21');
kInt = xlsread(file , 'Rates' , 'B24:F26');
perActInf = xlsread(file , 'Rates' , 'C36 : C38');
save('gcParams')
clear

% partners per year and condom usage
file = 'GC_HIV_ModelParameters.xlsx';
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
clear

end
