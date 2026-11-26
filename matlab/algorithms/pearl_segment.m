function result = pearl_segment(x,grid,edges,weights,cfg,initialPlanes)
% PEARL basic data+Potts energy, NO label costs; base MATLAB maxflow.
t=tic;hypotheses=0;
if nargin<6
    [planes,hypotheses]=propose(x,grid,cfg);
else
    planes=initialPlanes;
end
result.proposalSeconds=toc(t);
result.initialPlanes=planes;
scale=cfg.distanceThreshold^2;
% Internally divide the complete energy by threshold^2 for stable capacities.
% Returned histories are multiplied back, so reported E uses squared distances.
D=data_cost(x,planes,scale);
[~,labels]=min(D,[],2);labels=labels-1;
w=cfg.lambdaRatio*weights;
energy=compute_energy(labels,D,edges,w);
history=energy*scale;stages={'initial'};times=toc(t);
for iter=1:cfg.maxIterations
    startEnergy=energy;
    for sweep=1:cfg.expansionSweeps
        sweepEnergy=energy;
        for alpha=randperm(size(D,2))-1
            candidate=alpha_expansion_move(labels,D,edges,w,alpha);
            newEnergy=compute_energy(candidate,D,edges,w);
            if newEnergy<energy-cfg.energyTolerance
                labels=candidate;energy=newEnergy;
            end
        end
        if energy>=sweepEnergy-cfg.energyTolerance,break;end
    end
    history(end+1,1)=energy*scale; stages{end+1,1}='expand'; times(end+1,1)=toc(t); %#ok<AGROW>
    for k=1:size(planes,1)
        in=labels==k;
        if nnz(in)<3,continue;end
        p=fit_plane_svd(x(in,:));
        if any(~isfinite(p)),continue;end
        before=sum(D(in,k+1));
        after=sum((x(in,:)*p(1:3)'+p(4)).^2)/scale;
        if after<=before+1e-12,planes(k,:)=p;end
    end
    D=data_cost(x,planes,scale);
    energy=compute_energy(labels,D,edges,w);
    history(end+1,1)=energy*scale; stages{end+1,1}='refit';times(end+1,1)=toc(t); %#ok<AGROW>
    if startEnergy-energy<=cfg.energyTolerance,break;end
end
assert(all(diff(history)<=1e-8*max(1,history(1))), ...
    'PEARL energy increased unexpectedly.');
% Remove unused labels ONLY; do not reassign or merge after optimization.
active=unique(labels(labels>0));mapping=zeros(size(planes,1)+1,1);
mapping(active+1)=1:numel(active); labels=mapping(labels+1);
result.labels=labels;result.planes=planes(active,:);
result.seconds=toc(t);result.hypotheses=hypotheses;
result.energyHistory=history(:);result.energyStages=stages(:);
result.energyTimes=times(:);
end

function [planes,count]=propose(x,grid,cfg)
allPoints=(1:size(x,1))'; candidates=zeros(0,4);support=[];count=0;
for j=1:cfg.numProposals
    p=sample_plane(x,grid,allPoints,cfg);count=count+1;
    if any(~isfinite(p)),continue;end
    in=abs(x*p(1:3)'+p(4))<=cfg.distanceThreshold;
    if nnz(in)<cfg.minSupport,continue;end
    candidates(end+1,:)=p;support(end+1)=nnz(in); %#ok<AGROW>
end
[~,order]=sort(support,'descend');planes=zeros(0,4);
for j=order
    p=candidates(j,:);
    if ~isempty(planes)
        dotp=planes(:,1:3)*p(1:3)';
        angle=acosd(min(1,abs(dotp)));
        offset=abs(planes(:,4)-sign(dotp)*p(4));
        if any(angle<cfg.mergeAngleDegrees & offset<cfg.mergeOffset),continue;end
    end
    planes(end+1,:)=p; %#ok<AGROW>
    if size(planes,1)>=cfg.maxPlanes,break;end
end
end

function D=data_cost(x,planes,scale)
D=ones(size(x,1),size(planes,1)+1); % constant outlier cost=threshold^2
if ~isempty(planes)
    D(:,2:end)=(x*planes(:,1:3)'+planes(:,4)').^2/scale;
end
end

function e=compute_energy(labels,D,edges,w)
e=sum(D(sub2ind(size(D),(1:numel(labels))',labels+1)));
if ~isempty(edges),e=e+sum(w.*(labels(edges(:,1))~=labels(edges(:,2))));end
end
