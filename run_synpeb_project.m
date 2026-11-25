%% ECE 503: sequential RANSAC versus PEARL on SynPEB
% Run AFTER tuning phase='confirm' has written selected_parameters.mat.
% Demo uses TRAIN data only. Full uses all 30 test scans, 4 variants, 3 seeds.
clear; clc; close all;
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);
mode = 'full';                         % 'demo' or 'full'
cfg = synpeb_config(mode, projectDir);
summary = run_synpeb_experiment(cfg); %#ok<NASGU>


