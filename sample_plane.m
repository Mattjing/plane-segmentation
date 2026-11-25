function plane = sample_plane(x,grid,available,cfg)
% Mixed global and local non-collinear triplets, no labels.
plane=nan(1,4);
if numel(available)<3,return;end
if rand<cfg.proposalLocalFraction
    anchor=available(randi(numel(available)));
    [r,c]=find(grid==anchor,1);
    if ~isempty(r)
        patch=grid(max(1,r-cfg.localRadius):min(end,r+cfg.localRadius), ...
            max(1,c-cfg.localRadius):min(end,c+cfg.localRadius));
        nearby=intersect(patch(patch>0),available);
    else
        nearby=[];
    end
    if numel(nearby)>=3
        take=nearby(randperm(numel(nearby),3));
    else
        take=available(randperm(numel(available),3));
    end
else
    take=available(randperm(numel(available),3));
end
plane=fit_plane_svd(x(take,:));
end

