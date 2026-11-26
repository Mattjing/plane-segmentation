function synpeb_tuning_search(base,options)
% Staged coordinate/grid search, followed by separate training confirmation.
% Every candidate compares BOTH algorithms on identical XYZ input points.
% Primary score: mean matched IoU, averaged equally over the two methods.
assert(strcmp(base.split,'train'),'Tuning must use the train split.');
assert(ismember(options.phase,{'screen','confirm'}),'Unknown tuning phase.');
assert(isempty(intersect(options.screenScans,options.confirmScans)), ...
    'Screening and confirmation training scans must be disjoint.');
assert(all([options.screenScans options.confirmScans]>=1) && ...
    all([options.screenScans options.confirmScans]<=10),'Train scan IDs must be 1:10.');
base.saveFigures=false;base.runAblation=false;base.parametersTuned=false;
if ~isfolder(options.outputRoot),mkdir(options.outputRoot);end
phase=options.phase;
checkpoint=fullfile(options.outputRoot,[phase '_checkpoint.mat']);
shortlistFile=fullfile(options.outputRoot,'screen_shortlist.mat');
state.records=struct([]);state.configs={};state.keys={};
if isfile(checkpoint)
    old=load(checkpoint,'state');state=old.state;
    fprintf('Loaded %d completed %s configurations.\n',numel(state.keys),phase);
end

if strcmp(phase,'screen')
    current=base;
    current.scanIds=options.screenScans;
    current.variants=options.screenVariants;
    current.seeds=options.screenSeeds;
    % Coordinate searches avoid the much larger Cartesian product of all
    % proposal settings. These are heuristic stages, not a global grid search.
    proposalFields={'numProposals','maxPlanes','localRadius'};
    proposalValues={options.proposalCounts,options.planeCaps,options.localRadii};
    for j=1:numel(proposalFields)
        candidates={};
        for value=proposalValues{j}
            candidate=current;candidate.(proposalFields{j})=value;
            candidates{end+1}=candidate; %#ok<AGROW>
        end
        current=choose(candidates,['proposal_' proposalFields{j}]);
    end
    candidates={};
    for threshold=options.thresholds
        for ratio=options.lambdaRatios
            candidate=current;candidate.distanceThreshold=threshold;
            candidate.lambdaRatio=ratio;
            candidates{end+1}=candidate; %#ok<AGROW>
        end
    end
    current=choose(candidates,'threshold_smoothness');
    candidates={};
    for support=options.minimumSupports
        candidate=current;candidate.minSupport=support;
        candidates{end+1}=candidate; %#ok<AGROW>
    end
    choose(candidates,'minimum_support');
    % Rank only configurations evaluated on the current screening cases.
    eligible=[];
    for j=1:numel(state.configs)
        c=state.configs{j};
        if isequal(c.scanIds,options.screenScans) && ...
                isequal(c.variants,options.screenVariants) && isequal(c.seeds,options.screenSeeds)
            eligible(end+1)=j; %#ok<AGROW>
        end
    end
    ranking=rank_records(eligible);
    keep=ranking(1:min(options.shortlistCount,numel(ranking)));
    shortlist.configs=state.configs(keep);
    shortlist.records=state.records(keep);
    shortlist.options=options;
    save(shortlistFile,'shortlist');
    writetable(struct2table(shortlist.records),fullfile(options.outputRoot,'screen_shortlist.csv'));
    fprintf('\nScreening complete. No final parameter file was written.\n');
    fprintf('Set phase=''confirm'' in tune_synpeb_project.m and run again.\n');
    fprintf('Shortlist: %s\n',shortlistFile);
else
    assert(isfile(shortlistFile),'Run phase=''screen'' first.');
    loaded=load(shortlistFile,'shortlist');
    assert(isequal(loaded.shortlist.options.screenScans,options.screenScans) && ...
        isequal(loaded.shortlist.options.screenVariants,options.screenVariants) && ...
        isequal(loaded.shortlist.options.screenSeeds,options.screenSeeds), ...
        'Screen settings changed. Use a new outputRoot and screen again.');
    candidates=loaded.shortlist.configs;
    for j=1:numel(candidates)
        candidates{j}.scanIds=options.confirmScans;
        candidates{j}.variants=options.confirmVariants;
        candidates{j}.seeds=options.confirmSeeds;
        candidates{j}.dataRoot=base.dataRoot;
    end
    [winner,bestRecord]=choose(candidates,'confirmation');
    selected.distanceThreshold=winner.distanceThreshold;
    selected.lambdaRatio=winner.lambdaRatio;
    selected.selectionCriterion='mean training matched IoU averaged over both methods; ties: F1, runtime';
    selected.trainingScans=winner.scanIds;
    selected.screeningScans=options.screenScans;
    selected.variants=winner.variants;selected.seed=winner.seeds;
    selected.stride=winner.stride;
    selected.algorithmConfig=rmfield(winner,{'projectDir','dataRoot','outputRoot', ...
        'mode','split','scanIds','variants','seeds','runAblation','saveFigures', ...
        'parametersTuned','parameterFile'});
    selected.confirmationMetrics=bestRecord;
    if isfile(options.selectionPath)
        stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
        copyfile(options.selectionPath,[options.selectionPath '.' stamp '.bak']);
    end
    save(options.selectionPath,'selected','state');
    fprintf('\nConfirmed settings saved: %s\n',options.selectionPath);
    disp(selected.algorithmConfig);
    fprintf('Full test mode in run_synpeb_project.m will load these settings.\n');
end

    function [winner,record]=choose(candidates,stage)
        indices=zeros(1,numel(candidates));
        for ci=1:numel(candidates)
            c=candidates{ci};
            keyCfg=rmfield(c,{'outputRoot'});
            key=jsonencode(keyCfg);
            hit=find(strcmp(state.keys,key),1);
            if ~isempty(hit)
                indices(ci)=hit;
                fprintf('Reuse %s candidate %d/%d.\n',stage,ci,numel(candidates));
                continue;
            end
            fprintf('\n%s candidate %d/%d: proposals=%d, cap=%d, radius=%d, support=%d, threshold=%g, ratio=%g\n', ...
                stage,ci,numel(candidates),c.numProposals,c.maxPlanes,c.localRadius, ...
                c.minSupport,c.distanceThreshold,c.lambdaRatio);
            c.outputRoot=fullfile(options.outputRoot,phase,'candidate_runs');
            t=run_synpeb_experiment(c);
            r=strcmp(t.Method,'RANSAC');p=strcmp(t.Method,'PEARL');
            assert(any(r) && any(p),'Both methods must have results.');
            item.Stage=stage;item.NumProposals=c.numProposals;
            item.MaxPlanes=c.maxPlanes;item.LocalRadius=c.localRadius;
            item.MinSupport=c.minSupport;item.Threshold=c.distanceThreshold;
            item.LambdaRatio=c.lambdaRatio;
            item.RANSACIoU=mean(t.MeanMatchedIoU(r));
            item.PEARLIoU=mean(t.MeanMatchedIoU(p));
            item.JointIoU=(item.RANSACIoU+item.PEARLIoU)/2;
            item.JointF1=(mean(t.FScore(r))+mean(t.FScore(p)))/2;
            item.RANSACAccuracy=mean(t.Accuracy(r));
            item.PEARLAccuracy=mean(t.Accuracy(p));
            item.RANSACSeconds=mean(t.Seconds(r));item.PEARLSeconds=mean(t.Seconds(p));
            item.TotalSeconds=item.RANSACSeconds+item.PEARLSeconds;
            item.PEARLIoUSD=std(t.MeanMatchedIoU(p));
            item.CaseCount=nnz(p);
            if isempty(state.records),state.records=item;else,state.records(end+1)=item;end
            state.configs{end+1}=c;state.keys{end+1}=key;
            indices(ci)=numel(state.keys);
            save(checkpoint,'state');
            writetable(struct2table(state.records),fullfile(options.outputRoot,[phase '_results.csv']));
        end
        ranked=rank_records(indices);best=ranked(1);
        winner=state.configs{best};record=state.records(best);
        fprintf('Stage %s selected joint IoU %.4f, PEARL IoU %.4f.\n',stage,record.JointIoU,record.PEARLIoU);
    end

    function indices=rank_records(indices)
        records=state.records(indices);
        score=[-[records.JointIoU]' -[records.JointF1]' [records.TotalSeconds]'];
        assert(all(isfinite(score(:))),'Non-finite tuning scores.');
        [~,order]=sortrows(score,[1 2 3]);indices=indices(order);
    end
end
