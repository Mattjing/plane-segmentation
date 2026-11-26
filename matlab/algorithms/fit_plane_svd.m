function plane = fit_plane_svd(x)
% Unit normal, signed offset; NaNs mark degenerate samples.
plane=nan(1,4);
if size(x,1)<3, return; end
mu=mean(x,1); centered=x-mu;
[~,s,v]=svd(centered,0);
if size(s,1)<2 || s(2,2)<1e-10, return; end
n=v(:,end); n=n/norm(n);
[~,j]=max(abs(n)); if n(j)<0,n=-n;end
plane=[n' -mu*n];
end

