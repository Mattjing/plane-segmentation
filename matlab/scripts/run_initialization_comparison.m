%% ECE 503: sampled PEARL versus RANSAC-initialized PEARL.
% Run AFTER phase='confirm' of tune_synpeb_project has completed.
% This is independent of the tuner and does not change its scripts/results.
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
projectDir=setup_project();
mode='full';                          % 'full' or training-only 'demo'
if strcmp(mode,'full')
    assert(isfile(fullfile(projectDir,'config','selected_parameters.mat')), ...
        ['Confirmed parameters are not ready. Finish phase=''confirm'' ', ...
        'in tune_synpeb_project.m before running this experiment.']);
end
cfg=synpeb_config(mode,projectDir);
cfg.outputRoot=fullfile(projectDir,'results','runs','initialization_comparison',mode);
cfg.mode=['initialization_' mode];
cfg.runAblation=false;                % this experiment compares 3 pipelines
% Full: all 30 test scans, 4 variants, seeds [7 19 43], frozen settings.
% Demo: training scan 1, variants var1/var4, untuned demonstration only.
summary=compare_pearl_initializations(cfg); %#ok<NASGU>

