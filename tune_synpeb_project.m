%% ECE 503: staged parameter tuning using TRAIN data only.
% Screening is complete; the saved top-three shortlist is ready for confirmation.
% 1. Run with phase='screen' to make a shortlist on training scans 1 and 2.
% 2. Set phase='confirm' and run again to compare the top 3 configurations
%    on training scans 3:10, all variants, and three random seeds.
% 3. Only confirmation writes selected_parameters.mat for full test mode.
% Figures are disabled during tuning. Test data are never loaded here.
clear; clc;
projectDir=fileparts(mfilename('fullpath'));
addpath(projectDir);
phase='confirm';                       % 'screen' or 'confirm'
cfg=synpeb_config('tune',projectDir);

% Stronger starting proposal bank, fixed sampling density, longer refinement.
cfg.numProposals=1000;
cfg.maxPlanes=100;
cfg.localRadius=4;
cfg.minSupport=6;
cfg.distanceThreshold=0.02;
cfg.lambdaRatio=0.15;
cfg.stride=8;
cfg.maxIterations=15;
cfg.expansionSweeps=3;
cfg.saveFigures=false;
cfg.runAblation=false;

options.phase=phase;
options.outputRoot=fullfile(projectDir,'results','tuning_updated');
options.selectionPath=fullfile(projectDir,'selected_parameters.mat');
options.proposalCounts=[1000 3000 5000];
options.planeCaps=[60 100 150];
options.localRadii=[2 4 8];
options.thresholds=[0.005 0.01 0.02 0.04 0.08];
options.lambdaRatios=[0 0.02 0.05 0.15 0.4];
options.minimumSupports=[6 12 20];
options.screenScans=[1 2];
options.screenVariants={'var1','var4'};
options.screenSeeds=7;
options.confirmScans=3:10;
options.confirmVariants={'var1','var2','var3','var4'};
options.confirmSeeds=[7 19 43];
options.shortlistCount=3;

% Completed configurations are checkpointed. An interrupted configuration
% restarts, but completed ones are reused. If changing algorithm code, use
% a new outputRoot so old measurements are not reused.
synpeb_tuning_search(cfg,options);


