function labels = alpha_expansion_move(labels,D,edges,weights,alpha)
% Exact binary min-cut for ONE Potts alpha expansion; D col1=outlier.
% Binary state 0=keep, 1=switch. Source side=0, sink side=1.
% Original pair costs: E00,E01,E10,E11(=0). Their submodular
% decomposition is const + ai*xi + aj*xj + c*abs(xi-xj).
n=numel(labels); idx=sub2ind(size(D),(1:n)',labels+1);
u0=D(idx);u1=D(:,alpha+1);
if ~isempty(edges)
    i=edges(:,1);j=edges(:,2);
    e00=weights.*(labels(i)~=labels(j));
    e01=weights.*(labels(i)~=alpha);
    e10=weights.*(alpha~=labels(j));
    c=(e01+e10-e00)/2;
    assert(all(c>=-1e-10),'Non-submodular expansion move.'); c=max(c,0);
    u1=u1+accumarray(i,e10-e00-c,[n 1])+accumarray(j,e01-e00-c,[n 1]);
else
    i=[];j=[];c=[];
end
shift=min(u0,u1);u0=max(0,u0-shift);u1=max(0,u1-shift);
source=n+1;sink=n+2;nodes=(1:n)';
s=[repmat(source,n,1);nodes;i;j];
t=[nodes;repmat(sink,n,1);j;i];
capacity=[u1;u0;c;c];good=capacity>0;
g=digraph(s(good),t(good),capacity(good),n+2);
[~,~,~,sinkSide]=maxflow(g,source,sink);
switchNodes=sinkSide(sinkSide<=n);
labels(switchNodes)=alpha;
end

