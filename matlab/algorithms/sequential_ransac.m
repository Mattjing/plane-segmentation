function result = sequential_ransac(x,grid,cfg)
t=tic; n=size(x,1); labels=zeros(n,1); planes=zeros(0,4);
remaining=(1:n)'; hypotheses=0;
while numel(remaining)>=cfg.minSupport && size(planes,1)<cfg.maxPlanes
    best=false(numel(remaining),1); bestPlane=nan(1,4);
    bestResidual=inf;
    for k=1:cfg.ransacTrialsPerPlane
        p=sample_plane(x,grid,remaining,cfg); hypotheses=hypotheses+1;
        if any(~isfinite(p)),continue;end
        residual=abs(x(remaining,:)*p(1:3)'+p(4));
        in=residual<=cfg.distanceThreshold;
        if nnz(in)>nnz(best) || (nnz(in)==nnz(best) && sum(residual(in).^2)<bestResidual)
            best=in;bestPlane=p;bestResidual=sum(residual(in).^2);
        end
    end
    if nnz(best)<cfg.minSupport,break;end
    for refine=1:3
        p=fit_plane_svd(x(remaining(best),:));
        if any(~isfinite(p)),break;end
        in=abs(x(remaining,:)*p(1:3)'+p(4))<=cfg.distanceThreshold;
        if nnz(in)<cfg.minSupport,break;end
        best=in;bestPlane=p;
    end
    planes(end+1,:)=bestPlane; %#ok<AGROW>
    labels(remaining(best))=size(planes,1);
    remaining=remaining(~best);
end
result.labels=labels;result.planes=planes;
result.seconds=toc(t);result.hypotheses=hypotheses;
result.energyHistory=[];
end

