% Save all potentially calibrated parameters into structure
function [paramsAll] = genParamStruct()

numParams = 8;
paramsAll = cell(numParams,1);
% partner distribution (anal), [hiv status x risk]
paramsAll{1}.name = 'partnersAnal'; paramsAll{1}.length = 15;
    paramsAll{1}.ic = [[13; 8; 11; 16; 15] ; [8; 5; 9; 6; 7] ; [1; 1; 1; 1; 1]];
    paramsAll{1}.lb = [ones(paramsAll{1}.length*(1/3),1).*10; ...
                       ones(paramsAll{1}.length*(1/3),1).*3; ...
                       ones(paramsAll{1}.length*(1/3),1).*1];
    paramsAll{1}.ub = [ones(paramsAll{1}.length*(1/3),1).*50; ...
                       ones(paramsAll{1}.length*(1/3),1).*9.999; ...
                       ones(paramsAll{1}.length*(1/3),1).*2];
% partner distribution (oral), [hiv status x risk]
paramsAll{2}.name = 'partnersOral'; paramsAll{2}.length = 15;
    paramsAll{2}.ic = [[4; 4; 6; 5; 5] ; [3; 3; 2; 4; 4] ; [1; 1; 1; 1; 1]];
    paramsAll{2}.lb = [ones(paramsAll{1}.length*(1/3),1).*10; ...
                       ones(paramsAll{1}.length*(1/3),1).*3; ...
                       ones(paramsAll{1}.length*(1/3),1).*1];
    paramsAll{2}.ub = [ones(paramsAll{1}.length*(1/3),1).*50; ...
                       ones(paramsAll{1}.length*(1/3),1).*9.999; ...
                       ones(paramsAll{1}.length*(1/3),1).*2];
% acts, [1 x 1]
paramsAll{3}.name = 'acts'; paramsAll{3}.length = 1;
    paramsAll{3}.ic = 24;
    paramsAll{3}.lb = 2;
    paramsAll{3}.ub = 50;
% perActHiv, [1 x 1]
paramsAll{4}.name = 'perActHiv'; paramsAll{4}.length = 1;
    paramsAll{4}.ic = (0.084 * 10 ^ -2);
    paramsAll{4}.lb = 0.00001;
    paramsAll{4}.ub = 0.0186;
% perActGC_anal, [3 x 1], recUre , recPha , ureRec
paramsAll{5}.name = 'perActGC_anal'; paramsAll{5}.length = 3;
    paramsAll{5}.ic = [0.012*0.3; 0.0001*0.3; 0.035*0.3];
    paramsAll{5}.lb = [ones(paramsAll{5}.length,1).*0.00001];
    paramsAll{5}.ub = [ones(paramsAll{5}.length,1).*1.0];
% perActGC_oral, [3 x 1] , urePha , phaRec , phaUre
paramsAll{6}.name = 'perActGC_oral'; paramsAll{6}.length = 3;
    paramsAll{6}.ic = [0.00005*0.3; 0.00005*0.3; 0.0043*0.3];
    paramsAll{6}.lb = [ones(paramsAll{5}.length,1).*0.00001];
    paramsAll{6}.ub = [ones(paramsAll{5}.length,1).*1];
% cAssortTarget, [1 x risk], (0.01 to 1.0)
paramsAll{7}.name = 'cAssortTarget'; paramsAll{7}.length = 3;
    paramsAll{7}.ic = [0.3, 0.3, 0.4];
    paramsAll{7}.lb = [ones(paramsAll{5}.length,1).*0.01];
    paramsAll{7}.ub = [ones(paramsAll{5}.length,1).*0.9999];
% cAssort_init, [1 x risk], (0.01 to 0.99), percentage of cAssortTarget
paramsAll{8}.name = 'cAssort_init'; paramsAll{8}.length = 3;
    paramsAll{8}.ic = [0.2, 0.35, 0.45]; 
    paramsAll{8}.lb = [ones(paramsAll{5}.length,1).*0.01];
    paramsAll{8}.ub = [ones(paramsAll{5}.length,1).*0.99999];
