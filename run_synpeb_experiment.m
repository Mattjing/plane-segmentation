function summary = run_synpeb_experiment(cfg)
% Execute reproducible cases; failures remain explicit, never silently skipped.
assert(isfolder(cfg.dataRoot),'Dataset missing: %s',cfg.dataRoot);
assert(~isempty(which('digraph')) && exist('matchpairs','file')~=0, ...
    'Base MATLAB maxflow and matchpairs are required.');
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
runStamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
runDir=fullfile(cfg.outputRoot,runStamp);mkdir(runDir);
cfg.runDir=runDir;save(fullfile(runDir,'configuration.mat'),'cfg');
diary(fullfile(runDir,'console_log.txt'));cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('SynPEB %s run | %s | %s\n',cfg.mode,cfg.split,runDir);
fprintf('Defaults are untuned unless ParametersTuned=true.\n');
fprintf('Accuracy: GT IDs >=10 only; one-to-one ID matching; misses penalized.\n');
records=struct([]);audit=struct([]);caseNumber=0;
for scanId=cfg.scanIds
    for v=1:numel(cfg.variants)
        variant=cfg.variants{v};
        data=load_synpeb_case(cfg,scanId,variant);
        % GT coordinates are NEVER passed to algorithms. This check documents
        % why same-named gt files must not blindly supply variant point labels.
        auditItem=audit_reference(cfg,scanId,variant);
        if isempty(audit),audit=auditItem;else,audit(end+1)=auditItem;end %#ok<AGROW>
        for seed=cfg.seeds
            caseNumber=caseNumber+1;
            tag=sprintf('%s_pc%02d_%s_seed%d',cfg.split,scanId,variant,seed);
            fprintf('\n%s | %d sampled points, %d graph edges\n',tag,size(data.xyz,1),size(data.edges,1));
            rng(seed,'twister');
            ransac=sequential_ransac(data.xyz,data.grid,cfg);
            rng(seed,'twister');
            pearl=pearl_segment(data.xyz,data.grid,data.edges,data.weights,cfg);
            results={ransac,pearl};names={'RANSAC','PEARL'};configs={cfg,cfg};
            if cfg.runAblation
                ab=cfg;ab.lambdaRatio=0;rng(seed,'twister');
                ablation=pearl_segment(data.xyz,data.grid,data.edges,data.weights,ab,pearl.initialPlanes);
                ablation.seconds=ablation.seconds+pearl.proposalSeconds;
                ablation.hypotheses=pearl.hypotheses;
                ablation.proposalSeconds=pearl.proposalSeconds;
                results{end+1}=ablation;names{end+1}='PEARL_no_smoothness';configs{end+1}=ab;
            end
            for m=1:numel(results)
                res=results{m};metrics=evaluate_planes(data.reference,res.labels,cfg);
                folder=fullfile(runDir,tag);methodTag=[tag '_' names{m}];
                save_case_results(folder,methodTag,data,res,metrics,configs{m});
                record=make_record(scanId,variant,seed,names{m},res,metrics);
                if isempty(records),records=record;else,records(end+1)=record;end %#ok<AGROW>
                fprintf('%s: accuracy %.2f%% | F1 %.3f | IoU %.3f | planes %d/%d | %.2fs\n', ...
                    names{m},100*metrics.accuracy,metrics.fscore,metrics.meanMatchedIoU, ...
                    metrics.predPlaneCount,metrics.gtPlaneCount,res.seconds);
            end
            writetable(struct2table(records),fullfile(runDir,'metrics_per_case.csv'));
        end
    end
end
summary=struct2table(records);
writetable(struct2table(audit),fullfile(runDir,'reference_audit.csv'));
save(fullfile(runDir,'summary.mat'),'summary','cfg');
write_summary(summary,cfg,runDir);
fprintf('\nFinished %d cases. Results saved in:\n%s\n',caseNumber,runDir);
end

function r=make_record(scan,variant,seed,name,res,m)
r.Scan=scan;r.Variant=variant;r.Seed=seed;r.Method=name;
r.Accuracy=m.accuracy;r.MeanMatchedIoU=m.meanMatchedIoU;
r.Precision=m.precision;r.Recall=m.recall;r.FScore=m.fscore;
r.GTPlanes=m.gtPlaneCount;r.PredictedPlanes=m.predPlaneCount;
r.MissedPlanes=m.missedPlanes;r.Oversegmented=m.oversegmentedPlanes;
r.Undersegmented=m.undersegmentedPlanes;r.EvaluatedPoints=m.evaluatedPoints;
r.IgnoredPoints=m.ignoredPoints;r.Seconds=res.seconds;r.Hypotheses=res.hypotheses;
end

function a=audit_reference(cfg,scan,variant)
file=sprintf('pc_%02d.pcd',scan);
gt=read_synpeb_pcd(fullfile(cfg.dataRoot,cfg.split,'gt',file));
noisy=read_synpeb_pcd(fullfile(cfg.dataRoot,cfg.split,variant,file));
a.Scan=scan;a.Variant=variant;
a.PointCountEqual=size(gt.xyz,1)==size(noisy.xyz,1);
assert(a.PointCountEqual,'Point count differs between gt and variant.');
a.ReferenceIDDisagreement=mean(gt.labels~=noisy.labels);
a.ReferenceSource='variant intensity';
end

function write_summary(t,cfg,folder)
fid=fopen(fullfile(folder,'accuracy_summary.txt'),'w');assert(fid>=0);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'ECE 503 - Sequential RANSAC versus PEARL\n');
fprintf(fid,'Mode: %s; split: %s; parameters tuned: %d\n',cfg.mode,cfg.split,cfg.parametersTuned);
fprintf(fid,'Seeds: %s; stride: %d\n',mat2str(cfg.seeds),cfg.stride);
fprintf(fid,'Accuracy evaluated on plane IDs >=10, one-to-one maximum-overlap mapping.\n');
fprintf(fid,'Ignore/non-plane IDs <10 excluded; unassigned plane points count as errors.\n');
fprintf(fid,'No single aggregate confusion matrix: plane IDs are scene-specific.\n');
fprintf(fid,'Mean +/- SD below is across scan/variant/seed records (not independent samples).\n\n');
methods=unique(t.Method,'stable');
for j=1:numel(methods)
    use=strcmp(t.Method,methods{j});
    micro=sum(t.Accuracy(use).*t.EvaluatedPoints(use))/sum(t.EvaluatedPoints(use));
    fprintf(fid,'%s\n  Mean accuracy: %.4f +/- %.4f %%\n  Point-weighted accuracy: %.4f %%\n', ...
        methods{j},100*mean(t.Accuracy(use)),100*std(t.Accuracy(use)),100*micro);
    fprintf(fid,'  Mean matched IoU: %.4f; mean F-score: %.4f\n  Mean runtime: %.3f s\n\n', ...
        mean(t.MeanMatchedIoU(use)),mean(t.FScore(use)),mean(t.Seconds(use)));
end
fprintf(fid,'Results concern synthetic scans. Demo defaults are not final test results.\n');
end

