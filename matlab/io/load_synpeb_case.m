function data = load_synpeb_case(cfg,scanId,variant)
name = sprintf('pc_%02d.pcd',scanId);
pc = read_synpeb_pcd(fullfile(cfg.dataRoot,cfg.split,variant,name));
assert(pc.height>1 && pc.width*pc.height==size(pc.xyz,1), ...
    'Expected organized SynPEB cloud.');
rows=1:cfg.stride:pc.height; cols=1:cfg.stride:pc.width;
gridIndex=reshape(1:pc.width*pc.height,pc.width,pc.height)';
sampled=gridIndex(rows,cols);
ids=sampled(:);
valid=all(isfinite(pc.xyz(ids,:)),2) & any(pc.xyz(ids,:)~=0,2);
data.xyz=pc.xyz(ids(valid),:);
data.reference=pc.labels(ids(valid)); % variant's OWN labels, not gt pixel IDs
data.originalIndices=ids(valid);
data.grid=zeros(size(sampled)); data.grid(valid)=1:nnz(valid);
data.source=pc.filename;
data.rawPointCount=size(pc.xyz,1);
data.excludedInvalid=nnz(~valid);
assert(all(abs(data.reference-round(data.reference))<1e-5), ...
    'Non-integer reference labels: check dataset encoding.');
data.reference=round(data.reference);
% Graph uses geometry only. Never use reference labels to create edges.
a=data.grid(:,1:end-1); b=data.grid(:,2:end);
c=data.grid(1:end-1,:); d=data.grid(2:end,:);
u=[a(:);c(:)]; v=[b(:);d(:)]; ok=u>0 & v>0;
u=u(ok);v=v(ok);
dist=sqrt(sum((data.xyz(u,:)-data.xyz(v,:)).^2,2));
ok=dist<=cfg.graphMaxDistance;
data.edges=[u(ok) v(ok)];
data.weights=exp(-dist(ok).^2/(2*cfg.graphDistance^2));
data.evalMask=isfinite(data.reference) & data.reference>=10;
assert(any(data.evalMask),'No annotated plane points (IDs >=10).');
end

