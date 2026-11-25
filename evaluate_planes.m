function m = evaluate_planes(reference,predicted,cfg)
% Reference IDs <10 are ignore/non-plane; they NEVER influence fitting.
% Point accuracy uses maximum-overlap one-to-one matching (not many-to-one).
% Unmatched predictions stay wrong; predicted outliers on GT planes are wrong.
valid=isfinite(reference) & reference>=10;
gt=reference(valid);pr=predicted(valid);
[gtIds,~,g]=unique(gt);gCount=numel(gtIds);
% Include every extracted plane in detection denominators, even if it only
% covers ignored points. Such a plane cannot obtain an annotated GT match.
predIds=unique(predicted(predicted>0));pCount=numel(predIds);
[~,p]=ismember(pr,predIds);
C=zeros(gCount,pCount);
isPlane=p>0;
if pCount>0
    C=accumarray([g(isPlane) p(isPlane)],1,[gCount pCount]);
end
ng=accumarray(g,1,[gCount 1]);
np=zeros(pCount,1);
for k=1:pCount,np(k)=nnz(predicted==predIds(k));end
% IoU/instance detection considers only evaluable points in both areas.
npEval=sum(C,1)';
iou=C./max(ng+npEval'-C,1);
mapping=zeros(pCount,1);
if pCount>0
    pairs=matchpairs(-C,0); % minimize negative overlap = maximize correct count
    for k=1:size(pairs,1),mapping(pairs(k,2))=pairs(k,1);end
else
    pairs=zeros(0,2);
end
mapped=zeros(size(pr));mapped(isPlane)=mapping(p(isPlane));
unmatched=isPlane & mapped==0;mapped(unmatched)=gCount+1;
% Rows/cols: outlier, GT plane IDs, unmatched. GT has no outlier rows here.
cm=accumarray([g+1 mapped+1],1,[gCount+2 gCount+2]);
correct=nnz(mapped==g);
m.accuracy=correct/numel(g);
m.evaluatedPoints=numel(g);m.ignoredPoints=nnz(~valid);
m.gtIds=gtIds;m.predIds=predIds;m.mapping=mapping;
m.confusion=cm;m.referenceCompact=g;m.mappedPredicted=mapped;
m.evalMask=valid;m.overlap=C;m.iou=iou;
m.gtPlaneCount=gCount;m.predPlaneCount=pCount;
if pCount>0
    detectionCost=-(1+iou).*(iou>=cfg.iouThreshold);
    detectionPairs=matchpairs(detectionCost,0);
else
    detectionPairs=zeros(0,2);
end
tp=size(detectionPairs,1);
m.precision=tp/max(pCount,1);m.recall=tp/max(gCount,1);
m.fscore=2*m.precision*m.recall/max(m.precision+m.recall,eps);
m.missedPlanes=gCount-tp;
if pCount>0
    bestIoU=max(iou,[],2);
    m.oversegmentedPlanes=nnz(sum(C./max(ng,1)>=cfg.fragmentFraction,2)>=2);
    m.undersegmentedPlanes=nnz(sum(C./max(npEval',1)>=cfg.fragmentFraction,1)>=2);
else
    bestIoU=zeros(gCount,1);m.oversegmentedPlanes=0;m.undersegmentedPlanes=0;
end
m.meanBestIoU=mean(bestIoU); % diagnostic allows multiple GT to share a pred
matchedIoU=zeros(gCount,1);
for k=1:size(pairs,1),matchedIoU(pairs(k,1))=iou(pairs(k,1),pairs(k,2));end
m.meanMatchedIoU=mean(matchedIoU); % zero for each unmatched GT plane
m.mappingPairs=pairs;m.detectionPairs=detectionPairs;
m.predictedOutlierRate=nnz(pr==0)/numel(pr);
m.predSupport=np;
end

