function cfg = synpeb_config(mode, projectDir)
% All numerical defaults are STARTING VALUES, not validated/tuned results.
if nargin < 2, projectDir = fileparts(mfilename('fullpath')); end
cfg.projectDir = projectDir;
cfg.dataRoot = fullfile(projectDir, 'synpeb');
cfg.outputRoot = fullfile(projectDir, 'results', mode);
cfg.mode = mode;
cfg.split = 'train';
cfg.scanIds = 1;
cfg.variants = {'var1','var4'};
cfg.seeds = 7;
cfg.stride = 8;                        % organized-grid sampling, no GT use
cfg.distanceThreshold = 0.025;         % coordinate units; tune on train
cfg.lambdaRatio = 0.15;                % lambda = ratio * threshold^2
cfg.maxPlanes = 60;                    % computational cap, not true count
cfg.minSupport = 12;
cfg.numProposals = 500;
cfg.ransacTrialsPerPlane = 150;
cfg.proposalLocalFraction = 0.7;
cfg.localRadius = 4;                   % radius in sampled grid pixels
cfg.mergeAngleDegrees = 3;
cfg.mergeOffset = 0.015;
cfg.graphDistance = 0.15;              % edge scale, in coordinate units
cfg.graphMaxDistance = 0.30;           % omit edges across depth jumps
cfg.maxIterations = 8;
cfg.expansionSweeps = 2;
cfg.energyTolerance = 1e-9;
cfg.iouThreshold = 0.5;
cfg.fragmentFraction = 0.1;
cfg.runAblation = false;
cfg.saveFigures = true;
cfg.parametersTuned = false;
cfg.parameterFile = fullfile(projectDir,'selected_parameters.mat');
if strcmp(mode,'full')
    assert(isfile(cfg.parameterFile), ...
        'Tune on train first: run tune_synpeb_project.m.');
    s = load(cfg.parameterFile,'selected');
    assert(isfield(s.selected,'algorithmConfig'), ...
        'Parameter file is outdated. Repeat training parameter selection.');
    fields=fieldnames(s.selected.algorithmConfig);
    for j=1:numel(fields),cfg.(fields{j})=s.selected.algorithmConfig.(fields{j});end
    cfg.distanceThreshold = s.selected.distanceThreshold;
    cfg.lambdaRatio = s.selected.lambdaRatio;
    cfg.parametersTuned = true;
    cfg.split = 'test'; cfg.scanIds = 1:30;
    cfg.variants = {'var1','var2','var3','var4'};
    cfg.seeds = [7 19 43];
    cfg.runAblation = true;
elseif ~strcmp(mode,'demo') && ~strcmp(mode,'tune')
    error('Unknown run mode: %s',mode);
end
end

