function summary=compare_pearl_initializations(cfg)
% Three pipelines on identical points and frozen parameters, with GT evaluation
% only AFTER fitting. Calls existing code without modifying the running tuner.
assert(isfolder(cfg.dataRoot),'Missing SynPEB directory.');
if strcmp(cfg.split,'test')
    assert(cfg.parametersTuned,'Test comparison requires confirmed training parameters.');
end
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
folder=fullfile(cfg.outputRoot,stamp);mkdir(folder);
save(fullfile(folder,'configuration.mat'),'cfg');
% Snapshot confirmed settings so later parameter files cannot change this run.
if cfg.parametersTuned && isfile(cfg.parameterFile)
    copyfile(cfg.parameterFile,fullfile(folder,'confirmed_parameters_snapshot.mat'));
end
diary(fullfile(folder,'console_log.txt'));finish=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('Initialization comparison | split=%s | tuned=%d\n%s\n',cfg.split,cfg.parametersTuned,folder);
records=struct([]);paired=struct([]);
for scan=cfg.scanIds
    for vi=1:numel(cfg.variants)
        variant=cfg.variants{vi};data=load_synpeb_case(cfg,scan,variant);
        for seed=cfg.seeds
            tag=sprintf('%s_pc%02d_%s_seed%d',cfg.split,scan,variant,seed);
            fprintf('\n%s: %d points\n',tag,size(data.xyz,1));
            rng(seed,'twister');
            ransac=sequential_ransac(data.xyz,data.grid,cfg);
            rng(seed,'twister');
            sampled=pearl_segment(data.xyz,data.grid,data.edges,data.weights,cfg);
            sampled.initialization='sampled triplet proposals';
            sampled.initializationSeconds=sampled.proposalSeconds;
            sampled.refinementSeconds=max(0,sampled.seconds-sampled.proposalSeconds);
            rng(seed,'twister');
            hybrid=pearl_refine_ransac(data.xyz,data.edges,data.weights,cfg,ransac);
            % Energy is independently evaluated in original squared-distance units.
            baselineEnergy=energy(data,ransac,cfg);
            assert(abs(hybrid.energyHistory(1)-baselineEnergy)<=1e-8*max(1,baselineEnergy), ...
                'Hybrid must start at the exact baseline energy.');
            assert(energy(data,hybrid,cfg)<=baselineEnergy+1e-8*max(1,baselineEnergy), ...
                'Refinement must not increase the RANSAC-start energy.');
            ransac.initializationSeconds=ransac.seconds;ransac.refinementSeconds=0;
            results={ransac,sampled,hybrid};names={'RANSAC','PEARL_sampled','RANSAC_PEARL'};
            metrics=cell(1,3);
            for mi=1:3
                res=results{mi};m=evaluate_planes(data.reference,res.labels,cfg);metrics{mi}=m;
                save_case_results(fullfile(folder,tag),[tag '_' names{mi}],data,res,m,cfg);
                row=make_record(scan,variant,seed,names{mi},res,m,data,cfg);
                if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
                fprintf('%s: accuracy %.2f%% | matched IoU %.4f | F1 %.4f | total %.3fs | refinement %.3fs\n', ...
                    names{mi},100*m.accuracy,m.meanMatchedIoU,m.fscore,res.seconds,res.refinementSeconds);
            end
            pair.Scan=scan;pair.Variant=variant;pair.Seed=seed;
            pair.RANSACAccuracy=metrics{1}.accuracy;
            pair.SampledAccuracy=metrics{2}.accuracy;
            pair.HybridAccuracy=metrics{3}.accuracy;
            pair.HybridAccuracyGainPP=100*(metrics{3}.accuracy-metrics{1}.accuracy);
            pair.SampledAccuracyGainPP=100*(metrics{2}.accuracy-metrics{1}.accuracy);
            pair.RANSACIoU=metrics{1}.meanMatchedIoU;
            pair.SampledIoU=metrics{2}.meanMatchedIoU;
            pair.HybridIoU=metrics{3}.meanMatchedIoU;
            pair.HybridIoUGain=pair.HybridIoU-pair.RANSACIoU;
            pair.HybridF1Gain=metrics{3}.fscore-metrics{1}.fscore;
            pair.RANSACEnergy=baselineEnergy;
            pair.HybridEnergy=energy(data,hybrid,cfg);
            pair.EnergyDecrease=pair.RANSACEnergy-pair.HybridEnergy;
            pair.HybridRefinementSeconds=hybrid.refinementSeconds;
            if isempty(paired),paired=pair;else,paired(end+1)=pair;end %#ok<AGROW>
            summary=struct2table(records);pairedTable=struct2table(paired);
            writetable(summary,fullfile(folder,'initialization_metrics.csv'));
            writetable(pairedTable,fullfile(folder,'paired_refinement_changes.csv'));
            save(fullfile(folder,'checkpoint.mat'),'summary','pairedTable','cfg');
        end
    end
end
summary=struct2table(records);pairedTable=struct2table(paired);
write_summary(folder,summary,pairedTable,cfg);
if cfg.saveFigures,plot_comparison(folder,summary,pairedTable);end
fprintf('\nFinished. Figures and accuracy reports:\n%s\n',folder);
end

function row=make_record(scan,variant,seed,name,res,m,data,cfg)
row.Scan=scan;row.Variant=variant;row.Seed=seed;row.Method=name;
row.Accuracy=m.accuracy;row.MeanMatchedIoU=m.meanMatchedIoU;
row.Precision=m.precision;row.Recall=m.recall;row.FScore=m.fscore;
row.GTPlanes=m.gtPlaneCount;row.PredictedPlanes=m.predPlaneCount;
row.MissedPlanes=m.missedPlanes;row.Oversegmented=m.oversegmentedPlanes;
row.Undersegmented=m.undersegmentedPlanes;
row.EvaluatedPoints=m.evaluatedPoints;row.IgnoredPoints=m.ignoredPoints;
row.TotalSeconds=res.seconds;row.InitializationSeconds=res.initializationSeconds;
row.RefinementSeconds=res.refinementSeconds;row.Hypotheses=res.hypotheses;
row.FinalEnergy=energy(data,res,cfg);
if isempty(res.energyHistory)
    row.InitialEnergy=row.FinalEnergy;row.InitialPlaneCount=size(res.planes,1);
    row.AtIterationLimit=false;
else
    row.InitialEnergy=res.energyHistory(1);row.InitialPlaneCount=size(res.initialPlanes,1);
    row.AtIterationLimit=(numel(res.energyHistory)-1)/2>=cfg.maxIterations;
end
end

function e=energy(data,result,cfg)
labels=result.labels(:);assigned=labels>0;
cost=cfg.distanceThreshold^2*ones(size(labels));
if any(assigned)
    p=result.planes(labels(assigned),:);
    cost(assigned)=(sum(data.xyz(assigned,:).*p(:,1:3),2)+p(:,4)).^2;
end
e=sum(cost);
if ~isempty(data.edges)
    e=e+cfg.lambdaRatio*cfg.distanceThreshold^2*sum(data.weights.* ...
        (labels(data.edges(:,1))~=labels(data.edges(:,2))));
end
end

function write_summary(folder,t,p,cfg)
fid=fopen(fullfile(folder,'initialization_accuracy_summary.txt'),'w');assert(fid>=0);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'ECE 503 - PEARL initialization comparison\nSplit: %s; tuned: %d\n',cfg.split,cfg.parametersTuned);
fprintf(fid,'All metrics use one-to-one plane matching; truth IDs<10 ignored.\n');
fprintf(fid,'Mean +/- SD is descriptive across correlated scene/variant/seed records.\n');
methods=unique(t.Method,'stable');
for j=1:numel(methods)
    u=strcmp(t.Method,methods{j});
    fprintf(fid,'\n%s\n',methods{j});
    fprintf(fid,'  Accuracy: %.4f +/- %.4f %%\n',100*mean(t.Accuracy(u)),100*std(t.Accuracy(u)));
    fprintf(fid,'  Point-weighted accuracy: %.4f %%\n', ...
        100*sum(t.Accuracy(u).*t.EvaluatedPoints(u))/sum(t.EvaluatedPoints(u)));
    fprintf(fid,'  Mean matched IoU: %.5f; precision: %.5f; recall: %.5f; F1: %.5f\n', ...
        mean(t.MeanMatchedIoU(u)),mean(t.Precision(u)),mean(t.Recall(u)),mean(t.FScore(u)));
    fprintf(fid,'  Mean total: %.4fs; initialization: %.4fs; refinement: %.4fs\n', ...
        mean(t.TotalSeconds(u)),mean(t.InitializationSeconds(u)),mean(t.RefinementSeconds(u)));
    fprintf(fid,'  Mean initial/final planes: %.2f / %.2f; sampled hypotheses: %.2f\n', ...
        mean(t.InitialPlaneCount(u)),mean(t.PredictedPlanes(u)),mean(t.Hypotheses(u)));
    fprintf(fid,'  Runs reaching iteration cap: %d / %d\n',nnz(t.AtIterationLimit(u)),nnz(u));
end
fprintf(fid,'\nPaired RANSAC -> RANSAC_PEARL refinement\n');
fprintf(fid,'  Mean accuracy change: %+.4f percentage points\n',mean(p.HybridAccuracyGainPP));
fprintf(fid,'  Mean matched IoU change: %+.5f; F1 change: %+.5f\n',mean(p.HybridIoUGain),mean(p.HybridF1Gain));
fprintf(fid,'  Accuracy improved/tied/worsened: %d / %d / %d\n', ...
    nnz(p.HybridAccuracyGainPP>1e-10),nnz(abs(p.HybridAccuracyGainPP)<=1e-10),nnz(p.HybridAccuracyGainPP < -1e-10));
fprintf(fid,'  Mean energy decrease: %.8g\n',mean(p.EnergyDecrease));
fprintf(fid,'Hybrid total time INCLUDES the RANSAC initialization, counted once.\n');
fprintf(fid,'Both initialization strategies use the same frozen energy/graph settings.\n');
fprintf(fid,'Model-bank sizes and hypothesis budgets can differ; see per-case CSV.\n');
fprintf(fid,'Energy descent does not guarantee improved reference-label accuracy.\n');
end

function plot_comparison(folder,t,p)
methods=unique(t.Method,'stable');means=zeros(3,1);deviations=means;
for j=1:3,u=strcmp(t.Method,methods{j});means(j)=100*mean(t.Accuracy(u));deviations(j)=100*std(t.Accuracy(u));end
fig=figure('Visible','off','Color','w','Position',[50 50 1050 470]);
if isprop(fig,'Theme'),fig.Theme='light';end
tiledlayout(1,2,'TileSpacing','compact');
nexttile;bar(1:3,means);hold on;errorbar(1:3,means,deviations,'k.','LineWidth',1.2);
xticks(1:3);xticklabels(strrep(methods,'_',' '));xtickangle(20);
ylabel('Matched point accuracy (%)');title('Mean and descriptive standard deviation');grid on;
nexttile;scatter(p.RANSACIoU,p.HybridIoU,24,'filled');hold on;
plot([0 1],[0 1],'k--');axis square;xlim([0 1]);ylim([0 1]);grid on;
xlabel('RANSAC matched IoU');ylabel('RANSAC + PEARL matched IoU');title('Paired refinement effect');
axesList=findall(fig,'Type','axes');set(axesList,'Color','w','XColor','k','YColor','k');
set(findall(fig,'Type','text'),'Color','k');
exportgraphics(fig,fullfile(folder,'initialization_comparison.png'),'Resolution',160,'BackgroundColor','white');
savefig(fig,fullfile(folder,'initialization_comparison.fig'));close(fig);
end

