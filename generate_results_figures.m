function stats = generate_results_figures(projectDir,outputDir,comparisonDir,ablationDir)
% Section 7 figures from SAVED predictions and CSV measurements only.
% No segmentation or parameter tuning is rerun.
if nargin<1,projectDir=fileparts(mfilename('fullpath'));end
if nargin<2,outputDir=fullfile(projectDir,'results','report_figures');end
addpath(projectDir);if ~isfolder(outputDir),mkdir(outputDir);end
if nargin<3,comparisonDir=fullfile(projectDir,'results','initialization_comparison','full','20261124_022344_155');end
comparison=comparisonDir;
if nargin<4,ablationDir=fullfile(projectDir,'results','full','20261124_062337_283');end
ablation=ablationDir;
assert(isfolder(comparison) && isfolder(ablation), ...
    'Supply saved comparison and ablation run directories as arguments 3 and 4.');
T=readtable(fullfile(comparison,'initialization_metrics.csv'),'TextType','string');
A=readtable(fullfile(ablation,'metrics_per_case.csv'),'TextType','string');
frozen=load(fullfile(comparison,'configuration.mat'),'cfg');cfg=frozen.cfg;
cfg.projectDir=projectDir;cfg.dataRoot=fullfile(projectDir,'synpeb');
methods=["RANSAC","PEARL_sampled","RANSAC_PEARL"];
names={'RANSAC','Sampled PEARL','Hybrid'};colors=[0.12 0.37 0.65;0.85 0.37 0.10;0.18 0.56 0.39];
stats=struct('comparisonSource',comparison,'ablationSource',ablation,'matlabVersion',version);
for j=1:3
    U=T(T.Method==methods(j),:);assert(height(U)==360);
    stat.Method=char(methods(j));stat.Accuracy=mean(U.Accuracy);stat.AccuracySD=std(U.Accuracy);
    stat.MeanMatchedIoU=mean(U.MeanMatchedIoU);stat.F1=mean(U.FScore);
    stat.TotalSeconds=mean(U.TotalSeconds);stat.InitializationSeconds=mean(U.InitializationSeconds);
    stat.RefinementSeconds=mean(U.RefinementSeconds);stat.FinalPlanes=mean(U.PredictedPlanes);
    stat.InitialPlanes=mean(U.InitialPlaneCount);stat.IterationCapCount=nnz(U.AtIterationLimit);
    stats.comparison(j)=stat;
end
R=sortrows(T(T.Method==methods(1),:),{'Scan','Variant','Seed'});
S=sortrows(T(T.Method==methods(2),:),{'Scan','Variant','Seed'});
H=sortrows(T(T.Method==methods(3),:),{'Scan','Variant','Seed'});
assert(isequal(R(:,1:3),S(:,1:3))&&isequal(R(:,1:3),H(:,1:3)));
assert(all(H.FinalEnergy<=R.FinalEnergy+1e-8));
stats.pairedAccuracyImproved=nnz(H.Accuracy>R.Accuracy+1e-12);
stats.pairedAccuracyWorsened=nnz(H.Accuracy<R.Accuracy-1e-12);

% Figure 4: three-way descriptive comparisons, all from one timing session.
fig=makefig([1380 460]);tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
nexttile;bars(100*[stats.comparison.Accuracy],100*[stats.comparison.AccuracySD],names,colors);
ylabel('Matched point accuracy (%)');ylim([0 100]);title('(a) Mean and descriptive SD');
nexttile;bars([stats.comparison.MeanMatchedIoU],[],names,colors);
ylabel('Mean matched IoU');ylim([0 0.55]);title('(b) Equal weight per reference plane');
nexttile;bars([stats.comparison.TotalSeconds],[],names,colors);
ylabel('Algorithm runtime (s)');ylim([0 23]);title('(c) Mean total runtime');
savefigures(fig,outputDir,'results_three_way');

% Fixed, transparent failure selection rules, computed from the metrics.
joint=(R.Accuracy+S.Accuracy+H.Accuracy)/3;
joint(R.Seed~=7)=inf;[~,first]=min(joint);
[~,second]=min(H.Accuracy-R.Accuracy);selected=[first second];
stats.failureSelection={'lowest three-way mean point accuracy among seed-7 cases', ...
    'largest hybrid point-accuracy decrease over all 360 paired cases'};
fig=makefig([1450 780]);tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for row=1:2
    item=R(selected(row),:);scan=item.Scan;variant=char(item.Variant);seed=item.Seed;
    tag=sprintf('test_pc%02d_%s_seed%d',scan,variant,seed);
    data=load_synpeb_case(cfg,scan,variant);
    saved=cell(1,3);
    for j=1:3,saved{j}=load(fullfile(comparison,tag,[tag '_' char(methods(j)) '_results.mat']),'result','metrics');end
    nexttile;showmap(data,saved{1}.metrics,false);
    title(sprintf('Reference: %d planes\npc %02d, %s, seed %d',heightids(saved{1}.metrics),scan,variant,seed),'FontSize',13);
    for j=1:3
        nexttile;showmap(data,saved{j}.metrics,true);m=saved{j}.metrics;
        title(sprintf('%s: %.2f%% accuracy\nIoU %.3f | %d planes',names{j},100*m.accuracy,m.meanMatchedIoU,m.predPlaneCount),'FontSize',13);
    end
    stats.failureCases(row).scan=scan;stats.failureCases(row).variant=variant;stats.failureCases(row).seed=seed;
end
savefigures(fig,outputDir,'results_failure_cases');

% Figure 6: paired change and a counterexample to energy = accuracy.
fig=makefig([1080 440]);tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;good=H.Accuracy>=R.Accuracy;
scatter(100*R.Accuracy(good),100*H.Accuracy(good),15,colors(3,:),'filled');
scatter(100*R.Accuracy(~good),100*H.Accuracy(~good),48,[0.8 0.1 0.1],'filled');
plot([0 100],[0 100],'k--');axis equal;xlim([30 100]);ylim([30 100]);grid on;
xlabel('RANSAC accuracy (%)');ylabel('Hybrid accuracy (%)');
title('(a) Paired refinement effect');
legend({'Improved or tied','Worsened','Equal accuracy'},'Location','northwest','FontSize',10);
nexttile;tag=sprintf('test_pc%02d_%s_seed%d',R.Scan(second),R.Variant(second),R.Seed(second));
E=readtable(fullfile(comparison,tag,[tag '_RANSAC_PEARL_energy.csv']));
plot(E.Update,E.Energy,'-o','Color',colors(3,:),'LineWidth',1.8,'MarkerSize',4);hold on;
yline(R.FinalEnergy(second),'--');grid on;
legend({'Hybrid energy','Initial RANSAC energy'},'Location','northeast','FontSize',10);
xlabel('Optimizer update (expand or refit)');ylabel('Energy (m^2)');
title('(b) Energy decreases despite lower accuracy','FontSize',12);
stats.regressionEnergyInitial=E.Energy(1);stats.regressionEnergyFinal=E.Energy(end);
savefigures(fig,outputDir,'results_paired_energy');

% Figure 7: eta=0.4 versus eta=0, using the same saved initial plane bank.
P=sortrows(A(A.Method=="PEARL",:),{'Scan','Variant','Seed'});
N=sortrows(A(A.Method=="PEARL_no_smoothness",:),{'Scan','Variant','Seed'});
assert(height(P)==360&&isequal(P(:,1:3),N(:,1:3)));
seedRows=find(P.Seed==7);[~,order]=sort(P.Accuracy(seedRows));chosen=seedRows(order(ceil(numel(order)/2)));
scan=P.Scan(chosen);variant=char(P.Variant(chosen));seed=P.Seed(chosen);
tag=sprintf('test_pc%02d_%s_seed%d',scan,variant,seed);data=load_synpeb_case(cfg,scan,variant);
p=load(fullfile(ablation,tag,[tag '_PEARL_results.mat']),'result','metrics');
n=load(fullfile(ablation,tag,[tag '_PEARL_no_smoothness_results.mat']),'result','metrics');
assert(isequal(p.result.initialPlanes,n.result.initialPlanes),'Ablation plane banks differ.');
stats.ablationExample=struct('scan',scan,'variant',variant,'seed',seed,'sameInitialPlanes',true);
fig=makefig([1360 770]);tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
abNames={'Spatial prior','No spatial prior'};abColors=colors([2 1],:);
nexttile;bars(100*[mean(P.Accuracy),mean(N.Accuracy)],100*[std(P.Accuracy),std(N.Accuracy)],abNames,abColors);
ylabel('Matched point accuracy (%)');ylim([0 100]);title('(a) Mean and descriptive SD');
nexttile;bars([mean(P.MeanMatchedIoU),mean(N.MeanMatchedIoU)],[],abNames,abColors);
ylabel('Mean matched IoU');ylim([0 0.55]);title('(b) Reference-plane overlap');
nexttile;bars([mean(P.PredictedPlanes),mean(N.PredictedPlanes)],[],abNames,abColors);
yline(mean(P.GTPlanes),'k--','Reference mean','LabelHorizontalAlignment','left');
ylabel('Number of extracted planes');ylim([0 175]);title('(c) Active plane labels');
nexttile;showmap(data,p.metrics,false);title(sprintf('(d) Reference: %d planes',p.metrics.gtPlaneCount));
nexttile;showmap(data,p.metrics,true);title(sprintf('(e) Prior: %.2f%% | IoU %.3f',100*p.metrics.accuracy,p.metrics.meanMatchedIoU),'FontSize',13);
nexttile;showmap(data,n.metrics,true);title(sprintf('(f) No prior: %.2f%% | IoU %.3f',100*n.metrics.accuracy,n.metrics.meanMatchedIoU),'FontSize',13);
savefigures(fig,outputDir,'results_spatial_ablation');

% Retain the exact source tables and selection metadata for reproducibility.
writetable(T,fullfile(outputDir,'three_way_raw_metrics.csv'));
writetable(A,fullfile(outputDir,'ablation_raw_metrics.csv'));
writetable(struct2table(stats.comparison),fullfile(outputDir,'three_way_summary.csv'));
fid=fopen(fullfile(outputDir,'results_figure_statistics.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(stats,PrettyPrint=true),'char');
fprintf('Generated four result figures from saved predictions. No segmentation was rerun.\n%s\n',outputDir);
end
function n=heightids(m),n=m.gtPlaneCount;end
function fig=makefig(sz)
fig=figure('Visible','off','Color','w','Position',[50 50 sz]);
if isprop(fig,'Theme'),fig.Theme='light';end
end
function bars(values,sd,names,colors)
b=bar(1:numel(values),values,'FaceColor','flat');b.CData=colors;hold on;
if ~isempty(sd),errorbar(1:numel(values),values,sd,'k.','LineWidth',1.2);end
xticks(1:numel(values));xticklabels(replace(string(names),'Sampled PEARL','PEARL'));xtickangle(0);grid on;
for k=1:numel(values)
    label=sprintf('%.2f',values(k));if max(values)<1,label=sprintf('%.3f',values(k));end
    text(k,values(k),label,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12,'BackgroundColor','w');
end
end
function showmap(data,m,prediction)
if prediction,labels=m.mappedPredicted;else,labels=m.referenceCompact;end
cmap=[0.7 0.7 0.7;labelcolors(m.gtIds);0 0 0];
pointColors=repmat([0.88 0.88 0.88],size(data.xyz,1),1);
pointColors(m.evalMask,:)=cmap(labels+1,:);
img=ones([size(data.grid) 3]);valid=data.grid>0;idx=data.grid(valid);
for ch=1:3,slice=img(:,:,ch);slice(valid)=pointColors(idx,ch);img(:,:,ch)=slice;end
image(img);axis image;axis off;
end
function colors=labelcolors(ids)
hue=mod(double(ids(:))*0.618033988749895,1);
colors=hsv2rgb([hue 0.64*ones(numel(ids),1) 0.86*ones(numel(ids),1)]);
end
function savefigures(fig,outputDir,name)
set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',13,'Color','w','XColor','k','YColor','k');
set(findall(fig,'Type','text'),'FontName','Times New Roman','Color','k');
exportgraphics(fig,fullfile(outputDir,[name '.png']),'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,fullfile(outputDir,[name '.pdf']),'ContentType','image','Resolution',300,'BackgroundColor','white');
savefig(fig,fullfile(outputDir,[name '.fig']));close(fig);
end
