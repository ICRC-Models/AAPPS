% Run this to calibrate the GC model

function calibrateGCmodel()

close all; clear all; clc

%% Load parameters to be calibrated
% Load structure of all potentially calibrated parameters
[paramsAll] = genParamStruct();

% Choose indices in paramsAll cell array for the parameters you want to calibrate
pIdx = [1,2,3,4,5,6,7,8];
file = 'pIdx_21June19.dat';
csvwrite([pwd , '\' , file] , pIdx)

% Save info into paramsSub cell array for the parameters you want to calibrate
paramsSub = cell(length(pIdx),1);
p = 0;
startIdx = 1;
ic = [];
lb = [];
ub = [];
for s = 1 : length(pIdx)
    paramsSub{s} = paramsAll{pIdx(s)};
    p = p + paramsSub{s}.length; % find size of whole calibration matrix
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1)); % save indices for parameters in whole calibration matrix
    startIdx = startIdx + paramsSub{s}.length;
    ic = [ic; paramsSub{s}.ic]; % vector of all initial conditions
    lb = [lb; paramsSub{s}.lb]; % vector of all lower bounds
    ub = [ub; paramsSub{s}.ub]; % vector of all upper bounds
end

%% Find optimal parameter set
options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display' , 'iter' , 'PlotFcn' , @psplotbestf);
x = patternsearch(@mainCalibrate, ic , [] , [] , [] , [] , lb , ub , [] , options); % find minimum negative summed log-likelihood of function

%% Save calibrated parameters
file = 'hivGc_calib_21June19.dat';
csvwrite([pwd , '\' , file] , x)
