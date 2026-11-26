function test_synpeb_implementation()
% Algorithmic checks: min-cut vs exhaustive binary energies, plane geometry,
% label-permutation invariance, and penalties for misses/fragmentation.
rng(103,'twister');
for trial=1:25
    n=6;D=rand(n,4);edges=nchoosek(1:n,2);w=rand(size(edges,1),1);
    labels=randi([0 3],n,1);alpha=randi([0 3]);
    best=inf;
    for code=0:2^n-1
        cand=labels;cand(logical(bitget(code,1:n)))=alpha;
        best=min(best,energy(cand,D,edges,w));
    end
    actual=alpha_expansion_move(labels,D,edges,w,alpha);
    assert(abs(energy(actual,D,edges,w)-best)<1e-8, ...
        'Expansion cut did not match exhaustive binary optimum.');
end
x=[randn(100,2) 2*ones(100,1)];plane=fit_plane_svd(x);
assert(max(abs(x*plane(1:3)'+plane(4)))<1e-10,'Plane SVD incorrect.');
cfg=synpeb_config('demo');
gt=[10;10;10;10;11;11;11;11;3];
pr=[2;2;2;2;1;1;1;1;0];
m=evaluate_planes(gt,pr,cfg);
assert(m.accuracy==1 && m.fscore==1 && m.ignoredPoints==1);
fragmented=[1;1;2;2;3;3;3;3;0];m=evaluate_planes(gt,fragmented,cfg);
assert(m.accuracy==0.75,'One-to-one matching must penalize fragmentation.');
assert(m.oversegmentedPlanes==1);
m=evaluate_planes(gt,zeros(size(gt)),cfg);
assert(m.accuracy==0 && m.recall==0 && m.missedPlanes==2);
% End-to-end PEARL, including SVD updates, without reference information.
[xx,yy]=meshgrid(linspace(-1,1,12));
x=[xx(:) yy(:) 0.5*(xx(:)>0)+0.002*randn(numel(xx),1)];
cfg.maxPlanes=8;cfg.numProposals=60;cfg.maxIterations=3;
cfg.minSupport=8;grid=reshape(1:size(x,1),size(xx));
a=grid(:,1:end-1);b=grid(:,2:end);edges=[a(:) b(:)];
p=pearl_segment(x,grid,edges,ones(size(edges,1),1),cfg);
assert(all(diff(p.energyHistory)<=1e-8),'Energy must not increase.');
assert(numel(p.labels)==size(x,1));
cfg.lambdaRatio=0;
p0=pearl_segment(x,grid,edges,ones(size(edges,1),1),cfg,p.initialPlanes);
assert(all(diff(p0.energyHistory)<=1e-8),'Ablation energy must not increase.');
empty=pearl_segment(x,grid,edges,ones(size(edges,1),1),cfg,zeros(0,4));
assert(all(empty.labels==0),'Empty proposal set should produce outliers.');
fprintf('PASS: 25 exhaustive alpha-expansion checks; SVD; matching; PEARL descent.\n');
end

function e=energy(labels,D,edges,w)
e=sum(D(sub2ind(size(D),(1:numel(labels))',labels+1))) ...
    +sum(w.*(labels(edges(:,1))~=labels(edges(:,2))));
end

