function save_case_results(folder,tag,data,result,metrics,cfg)
if ~isfolder(folder),mkdir(folder);end
save(fullfile(folder,[tag '_results.mat']),'result','metrics','cfg','-v7');
% Retain raw, original IDs for independently reproducing the evaluation.
labels=table(data.originalIndices,data.reference,result.labels, ...
    'VariableNames',{'OriginalPointIndex','ReferenceID','PredictedID'});
writetable(labels,fullfile(folder,[tag '_labels.csv']));
fid=fopen(fullfile(folder,[tag '_accuracy.txt']),'w');
assert(fid>=0,'Cannot write report.');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\nSource: %s\n',tag,data.source);
fprintf(fid,'Parameters tuned on train: %d\n',cfg.parametersTuned);
fprintf(fid,'Point accuracy (one-to-one overlap matching): %.6f%%\n',100*metrics.accuracy);
fprintf(fid,'Evaluated plane points: %d; ignored label<10 points: %d\n', ...
    metrics.evaluatedPoints,metrics.ignoredPoints);
fprintf(fid,'GT planes: %d; estimated planes: %d\n',metrics.gtPlaneCount,metrics.predPlaneCount);
fprintf(fid,'Precision: %.6f; recall: %.6f; F-score: %.6f (IoU >= %.2f)\n', ...
    metrics.precision,metrics.recall,metrics.fscore,cfg.iouThreshold);
fprintf(fid,'Mean matched plane IoU: %.6f (unmatched GT=0)\n',metrics.meanMatchedIoU);
fprintf(fid,'Missed: %d; oversegmented GT: %d; undersegmented predicted: %d\n', ...
    metrics.missedPlanes,metrics.oversegmentedPlanes,metrics.undersegmentedPlanes);
fprintf(fid,'Fragmentation rule: >=2 overlaps, each >=%.2f of respective region\n',cfg.fragmentFraction);
fprintf(fid,'Predicted outlier fraction on annotated planes: %.6f\n',metrics.predictedOutlierRate);
fprintf(fid,'Algorithm time including proposals: %.6f seconds\n',result.seconds);
fprintf(fid,'Sampled hypotheses: %d\n',result.hypotheses);
fprintf(fid,'Distance threshold: %.6g; lambda: %.6g; stride: %d\n', ...
    cfg.distanceThreshold,cfg.lambdaRatio*cfg.distanceThreshold^2,cfg.stride);
fprintf(fid,'GT labels are evaluation only. Unmatched predicted planes count as errors.\n');
fprintf(fid,'Confusion rows=reference; columns=matched prediction.\n');
names=[{'Outlier'};arrayfun(@(v)sprintf('GT %d',v),metrics.gtIds,'UniformOutput',false);{'Unmatched'}];
cm=metrics.confusion;
writematrix(cm,fullfile(folder,[tag '_confusion_counts.csv']));
fidLabels=fopen(fullfile(folder,[tag '_confusion_labels.txt']),'w');
for j=1:numel(names),fprintf(fidLabels,'%d: %s\n',j,names{j});end
fclose(fidLabels);
if ~isempty(result.energyHistory)
    et=table((0:numel(result.energyHistory)-1)',result.energyStages, ...
        result.energyHistory,result.energyTimes, ...
        'VariableNames',{'Update','Stage','Energy','ElapsedSeconds'});
    writetable(et,fullfile(folder,[tag '_energy.csv']));
end
if ~cfg.saveFigures,return;end
% Export graphics without depending on Statistics/Computer Vision Toolbox.
fig=figure('Visible','off','Color','w','Position',[50 50 1200 540]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;imagesc(cm);axis image;colorbar;
title(sprintf('Counts | accuracy %.2f%%',100*metrics.accuracy));
xlabel('Matched predicted plane');ylabel('Reference plane');
if numel(names)<=24
    tickIds=1:numel(names);
else
    tickIds=unique([1,4:ceil(numel(names)/24):numel(names)-3,numel(names)]);
end
xticks(tickIds);yticks(tickIds);
xticklabels(names(tickIds));yticklabels(names(tickIds));xtickangle(90);
nexttile;normCm=cm./max(sum(cm,2),1);imagesc(normCm,[0 1]);axis image;colorbar;
title('Row-normalized confusion');xlabel('Matched predicted plane');ylabel('Reference plane');
xticks(tickIds);yticks(tickIds);
xticklabels(names(tickIds));yticklabels(names(tickIds));xtickangle(90);
if numel(names)<=20
    for r=1:size(cm,1)
        for c=1:size(cm,2)
            if cm(r,c)>0
                text(c,r,sprintf('%.2f',normCm(r,c)),'HorizontalAlignment','center', ...
                    'FontSize',7,'Color',contrast_color(normCm(r,c)),'Tag','CellValue');
            end
        end
    end
end
sgtitle(strrep(tag,'_',' '));
style_figure(fig);
exportgraphics(fig,fullfile(folder,[tag '_confusion.png']),'Resolution',160,'BackgroundColor','white');
savefig(fig,fullfile(folder,[tag '_confusion.fig']));close(fig);
fig=figure('Visible','off','Color','w','Position',[50 50 1250 470]);
tiledlayout(1,2,'TileSpacing','compact');
% Matched GT colors are shared; outlier=gray and unmatched=black.
cmap=[.65 .65 .65;lines(numel(metrics.gtIds));0 0 0];
gtColor=zeros(size(data.reference));gtColor(metrics.evalMask)=metrics.referenceCompact;
predColor=zeros(size(data.reference));predColor(metrics.evalMask)=metrics.mappedPredicted;
nexttile;scatter3(data.xyz(:,1),data.xyz(:,2),data.xyz(:,3),5,cmap(gtColor+1,:),'.');
axis equal;grid on;view(35,20);title('Reference (ignored points gray)');xlabel('X');ylabel('Y');zlabel('Z');
nexttile;scatter3(data.xyz(:,1),data.xyz(:,2),data.xyz(:,3),5,cmap(predColor+1,:),'.');
axis equal;grid on;view(35,20);title('Prediction (unmatched planes black)');xlabel('X');ylabel('Y');zlabel('Z');
sgtitle(strrep(tag,'_',' '));
style_figure(fig);
exportgraphics(fig,fullfile(folder,[tag '_segmentation.png']),'Resolution',150,'BackgroundColor','white');close(fig);
if ~isempty(result.energyHistory)
    fig=figure('Visible','off','Color','w');
    plot(0:numel(result.energyHistory)-1,result.energyHistory,'-o','LineWidth',1.5);
    xlabel('Update (expand/refit)');ylabel('Energy (squared coordinate units)');grid on;
    title('PEARL energy descent');
    style_figure(fig);
    exportgraphics(fig,fullfile(folder,[tag '_energy.png']),'Resolution',150,'BackgroundColor','white');close(fig);
end
end

function c=contrast_color(v)
if v<.25,c=[1 1 1];else,c=[0 0 0];end
end

function style_figure(fig)
% Explicit colors avoid inheriting a user's dark desktop graphics theme.
if isprop(fig,'Theme'),fig.Theme='light';end
set(fig,'Color','w');
axesList=findall(fig,'Type','axes');
for j=1:numel(axesList)
    set(axesList(j),'Color','w','XColor','k','YColor','k','ZColor','k','FontSize',9);
end
texts=findall(fig,'Type','text');
for j=1:numel(texts)
    if ~strcmp(get(texts(j),'Tag'),'CellValue'),set(texts(j),'Color','k');end
end
set(findall(fig,'Type','colorbar'),'Color','k');
end

