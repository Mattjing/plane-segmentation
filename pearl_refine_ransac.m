function result=pearl_refine_ransac(x,edges,weights,cfg,ransac)
% PEARL's basic data+Potts optimizer starting from EXACT RANSAC labels/models.
% Independent helper: does not change pearl_segment.m while tuning runs.
% No GT labels, extra proposals, label costs, or final support filtering.
t=tic;n=size(x,1);planes=ransac.planes;labels=ransac.labels(:);
assert(numel(labels)==n && all(isfinite(labels)) && ...
    all(labels==round(labels)) && all(labels>=0 & labels<=size(planes,1)), ...
    'RANSAC labels must index the supplied plane bank, with outlier=0.');
assert(size(planes,2)==4 && all(isfinite(planes(:))),'Invalid initial plane bank.');
assert(cfg.distanceThreshold>0 && cfg.lambdaRatio>=0,'Invalid energy settings.');
if ~isempty(planes)
    assert(all(abs(sqrt(sum(planes(:,1:3).^2,2))-1)<1e-6), ...
        'Initial plane normals must have unit length.');
end
result.initialPlanes=planes;result.initialLabels=labels;
result.initialization='exact RANSAC planes and assignments';
scale=cfg.distanceThreshold^2;w=cfg.lambdaRatio*weights;
D=costs(x,planes,scale);energy=value(labels,D,edges,w);
history=energy*scale;stages={'RANSAC_start'};times=0;
for iteration=1:cfg.maxIterations
    startEnergy=energy;
    for sweep=1:cfg.expansionSweeps
        sweepEnergy=energy;
        for alpha=randperm(size(D,2))-1
            proposed=alpha_expansion_move(labels,D,edges,w,alpha);
            newEnergy=value(proposed,D,edges,w);
            if newEnergy<energy-cfg.energyTolerance
                labels=proposed;energy=newEnergy;
            end
        end
        if sweepEnergy-energy<=cfg.energyTolerance,break;end
    end
    history(end+1,1)=energy*scale;stages{end+1,1}='expand';times(end+1,1)=toc(t); %#ok<AGROW>
    for k=1:size(planes,1)
        in=labels==k;if nnz(in)<3,continue;end
        plane=fit_plane_svd(x(in,:));
        if any(~isfinite(plane)),continue;end
        after=sum((x(in,:)*plane(1:3)'+plane(4)).^2)/scale;
        before=sum(D(in,k+1));
        if after<=before+1e-12,planes(k,:)=plane;end
    end
    D=costs(x,planes,scale);energy=value(labels,D,edges,w);
    history(end+1,1)=energy*scale;stages{end+1,1}='refit';times(end+1,1)=toc(t); %#ok<AGROW>
    if startEnergy-energy<=cfg.energyTolerance,break;end
end
assert(all(diff(history)<=1e-8*max(1,history(1))), ...
    'Hybrid PEARL energy unexpectedly increased.');
% Compact used labels; preserve all assignments and the optimized energy.
active=unique(labels(labels>0));map=zeros(size(planes,1)+1,1);
map(active+1)=1:numel(active);
result.labels=map(labels+1);result.planes=planes(active,:);
result.refinementSeconds=toc(t);
result.initializationSeconds=ransac.seconds;
result.seconds=ransac.seconds+result.refinementSeconds;
result.proposalSeconds=ransac.seconds;
result.hypotheses=ransac.hypotheses;
result.energyHistory=history(:);result.energyStages=stages(:);
result.energyTimes=ransac.seconds+times(:);
end

function D=costs(x,planes,scale)
D=ones(size(x,1),size(planes,1)+1);
if ~isempty(planes)
    D(:,2:end)=(x*planes(:,1:3)'+planes(:,4)').^2/scale;
end
end

function e=value(labels,D,edges,w)
e=sum(D(sub2ind(size(D),(1:numel(labels))',labels+1)));
if ~isempty(edges),e=e+sum(w.*(labels(edges(:,1))~=labels(edges(:,2))));end
end

