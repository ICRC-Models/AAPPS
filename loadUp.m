function loadUp()
hivStatus = 5; % (1) negative , (2) infectious , (3) tested , (4) treated...
% (5) on PreP, HIV immune
stiTypes = 4; % (1) none (2) gonorrhea, (3) chlamydia, (4) syphilis
sites = 3; % (1) rectal, (2) urethral, (3) pharyngeal
% risk = 3;
file = 'GC_HIV_ModelParameters.xlsx';
kBorn = xlsread(file , 'Rates' , 'B29');
kDie = xlsread(file , 'Rates' , 'B30');
save('genParams');
clear
% dimensions for matrices [sites x hivStatus]
pSite = xlsread(file , 'Rates' , 'B4:F6');
gcTreat = xlsread(file , 'Rates' , 'B9:F11');
k_gcPS = xlsread(file , 'Rates' , 'B14:F16');
gcClear = xlsread(file , 'Rates' , 'B19:F21');
kInt = xlsread(file , 'Rates' , 'B24:F26');
save('gcParams')
end
