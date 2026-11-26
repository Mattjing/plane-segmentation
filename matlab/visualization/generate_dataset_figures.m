function stats = generate_dataset_figures(projectDir, outputDir)
% Reproducible DATASET illustrations for Section 6. No algorithms are rerun.
% Run from Final Project: generate_dataset_figures
% Requires the project's read_synpeb_pcd, load_synpeb_case and fit_plane_svd.
% All visualization examples are TRAIN scans chosen by index, not test scores.
if nargin<1, projectDir=fileparts(fileparts(fileparts(mfilename('fullpath')))); end
if nargin<2, outputDir=fullfile(projectDir,'results','figures','dataset'); end
addpath(genpath(fullfile(projectDir,'matlab')));
if ~isfolder(outputDir), mkdir(outputDir); end
cfg=synpeb_config('full',projectDir);cfg.split='train';
stats=struct;stats.matlabVersion=version;stats.stride=cfg.stride;
stats.graphDistance=cfg.graphDistance;stats.graphMaxDistance=cfg.graphMaxDistance;
stats.examples='train pc_01 and pc_02; variants of train pc_01';
stats.selectedParameters=struct('threshold',cfg.distanceThreshold, ...
    'lambdaRatio',cfg.lambdaRatio,'numProposals',cfg.numProposals, ...
    'maxPlanes',cfg.maxPlanes,'minSupport',cfg.minSupport, ...
    'localRadius',cfg.localRadius,'maxIterations',cfg.maxIterations, ...
    'expansionSweeps',cfg.expansionSweeps);

% Count every installed file without decoding the entire benchmark.
records=struct([]);
for split={'train','test'}
    for name={'gt','var1','var2','var3','var4'}
        files=dir(fullfile(cfg.dataRoot,split{1},name{1},'pc_*.pcd'));
        a.Split=split{1};a.Directory=name{1};a.FileCount=numel(files);
        if isempty(records),records=a;else,records(end+1)=a;end %#ok<AGROW>
    end
end
writetable(struct2table(records),fullfile(outputDir,'dataset_inventory.csv'));
stats.inventory=records;

% Figure 1: two examples; geometry uses exactly the experiment stride.
fig=make_figure([1420 490]);tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for scan=1:2
    pc=read_synpeb_pcd(fullfile(cfg.dataRoot,'train','gt',sprintf('pc_%02d.pcd',scan)));
    d=load_synpeb_case(cfg,scan,'gt');
    nexttile;cloud(d.xyz,d.reference);title(sprintf('(%c) Train pc %02d: reference geometry','a'+scan-1,scan));
    if scan==1,firstPC=pc;end
    stats.overview(scan).scan=scan;stats.overview(scan).width=pc.width;
    stats.overview(scan).height=pc.height;stats.overview(scan).rawPoints=size(pc.xyz,1);
    stats.overview(scan).sampledPoints=size(d.xyz,1);
    stats.overview(scan).sampledPlaneCount=numel(unique(d.reference(d.reference>=10)));
end
nexttile;image(label_image(firstPC));axis image;axis ij;
xlabel('Column');ylabel('Row');title('(c) Train pc 01: organized reference labels');
title(tl,'SynPEB examples: colored by supplied reference plane membership','FontSize',18);
save_figure(fig,outputDir,'dataset_overview');

% Figure 2: all four measurement variants, each with its own annotations.
variants={'var1','var2','var3','var4'};cases=cell(1,4);residuals=cell(1,4);
for v=1:4
    d=load_synpeb_case(cfg,1,variants{v});cases{v}=d;
    ids=unique(d.reference(d.evalMask));r=nan(size(d.reference));
    for id=ids(:)'
        support=d.reference==id;p=fit_plane_svd(d.xyz(support,:));
        if all(isfinite(p)),r(support)=abs(d.xyz(support,:)*p(1:3)'+p(4));end
    end
    residuals{v}=r(isfinite(r));
    stats.variants(v).variant=variants{v};
    stats.variants(v).sampledPoints=size(d.xyz,1);
    stats.variants(v).evaluatedPoints=nnz(d.evalMask);
    stats.variants(v).ignoredPoints=nnz(~d.evalMask);
    stats.variants(v).planeCount=numel(ids);
    stats.variants(v).medianAnnotatedFitResidual=median(residuals{v});
    stats.variants(v).graphEdges=size(d.edges,1);
end
fig=make_figure([1420 490]);tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
shared=[cases{1}.xyz;cases{4}.xyz];lims=[min(shared,[],1);max(shared,[],1)];
for v=[1 4]
    nexttile;cloud(cases{v}.xyz,cases{v}.reference);
    xlim(lims(:,1)');ylim(lims(:,2)');zlim(lims(:,3)');
    title(sprintf('(%c) Train pc 01: %s',char('a'+(v==4)),variants{v}));
end
nexttile;hold on;colors=lines(4);cap=0;
for v=1:4
    values=sort(1000*residuals{v});
    plot(values,(1:numel(values))'/numel(values)*100,'Color',colors(v,:),'LineWidth',2.2);
    cap=max(cap,values(max(1,ceil(0.99*numel(values)))));
end
xlim([0 cap]);ylim([0 100]);grid on;
xlabel('Orthogonal residual (mm)');ylabel('Cumulative fraction (%)');
legend(variants,'Location','southeast','FontSize',12);
title('(c) Within annotated plane supports');
title(tl,'Supplied measurement variants of one training scan','FontSize',18);
save_figure(fig,outputDir,'dataset_variants');

% Figure 3: full image, experimental subsampling, graph crop.
pc=read_synpeb_pcd(fullfile(cfg.dataRoot,'train','var1','pc_01.pcd'));
d=cases{1};img=label_image(pc);rows=1:cfg.stride:pc.height;cols=1:cfg.stride:pc.width;
sampled=img(rows,cols,:);cropRows=23:40;cropCols=23:40;
fig=make_figure([1400 490]);tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
nexttile;image(img);axis image;axis ij;xlabel('Column');ylabel('Row');
title(sprintf('(a) Organized labels: %d by %d',pc.height,pc.width));
nexttile;image(sampled);axis image;axis ij;hold on;
rectangle('Position',[cropCols(1)-0.5 cropRows(1)-0.5 numel(cropCols) numel(cropRows)], ...
    'EdgeColor','k','LineWidth',2);
xlabel('Sampled column');ylabel('Sampled row');
title(sprintf('(b) Stride %d: %d valid points',cfg.stride,size(d.xyz,1)));
nexttile;hold on;
[rr,cc]=find(d.grid>0);validIds=d.grid(d.grid>0);
position=zeros(size(d.xyz,1),2);position(validIds,:)=[cc rr];
inside=position(:,1)>=cropCols(1)&position(:,1)<=cropCols(end)& ...
    position(:,2)>=cropRows(1)&position(:,2)<=cropRows(end);
edges=d.edges(inside(d.edges(:,1))&inside(d.edges(:,2)),:);
xx=[position(edges(:,1),1) position(edges(:,2),1) nan(size(edges,1),1)]';
yy=[position(edges(:,1),2) position(edges(:,2),2) nan(size(edges,1),1)]';
plot(xx(:),yy(:),'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
scatter(position(inside,1),position(inside,2),32,label_colors(d.reference(inside)),'filled');
axis equal;axis ij;xlim([cropCols(1)-1 cropCols(end)+1]);ylim([cropRows(1)-1 cropRows(end)+1]);
xlabel('Sampled column');ylabel('Sampled row');box on;
title('(c) Geometry-based graph in boxed region');
title(tl,'Experimental preprocessing and spatial neighborhood','FontSize',18);
save_figure(fig,outputDir,'dataset_preprocessing');

stats.graphExample=struct('scan',1,'variant','var1','nodes',size(d.xyz,1), ...
    'edges',size(d.edges,1),'rawPoints',d.rawPointCount, ...
    'invalidSampledPoints',d.excludedInvalid,'gridRows',size(d.grid,1),'gridColumns',size(d.grid,2));
fid=fopen(fullfile(outputDir,'dataset_figure_statistics.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(stats,PrettyPrint=true),'char');
fprintf('Saved three MATLAB figures and data audit to %s\n',outputDir);
end

function fig=make_figure(position)
fig=figure('Visible','off','Color','w','Position',[50 50 position]);
if isprop(fig,'Theme'),fig.Theme='light';end
end
function cloud(x,labels)
scatter3(x(:,1),x(:,2),x(:,3),10,label_colors(labels),'filled');
axis equal;axis tight;grid on;view(35,25);
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');
end
function colors=label_colors(labels)
labels=double(labels(:));colors=repmat([0.78 0.78 0.78],numel(labels),1);
ok=isfinite(labels)&labels>=10;
hue=mod(labels(ok)*0.618033988749895,1);
colors(ok,:)=hsv2rgb([hue 0.64*ones(nnz(ok),1) 0.86*ones(nnz(ok),1)]);
end
function img=label_image(pc)
colors=label_colors(pc.labels);img=zeros(pc.height,pc.width,3);
for channel=1:3,img(:,:,channel)=reshape(colors(:,channel),pc.width,pc.height)';end
end
function save_figure(fig,outputDir,name)
ax=findall(fig,'Type','axes');set(ax,'FontName','Times New Roman','FontSize',14, ...
    'Color','w','XColor','k','YColor','k','ZColor','k');
set(findall(fig,'Type','text'),'FontName','Times New Roman','Color','k');
exportgraphics(fig,fullfile(outputDir,[name '.png']),'Resolution',240,'BackgroundColor','white');
% Dense 3-D scatter is exported as a high-resolution image inside the PDF;
% native MATLAB objects remain editable in the accompanying FIG file.
exportgraphics(fig,fullfile(outputDir,[name '.pdf']),'ContentType','image', ...
    'Resolution',300,'BackgroundColor','white');
savefig(fig,fullfile(outputDir,[name '.fig']));close(fig);
end
