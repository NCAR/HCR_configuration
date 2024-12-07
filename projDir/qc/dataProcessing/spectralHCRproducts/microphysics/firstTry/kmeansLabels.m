function micLabels=kmeansLabels(momentsTime,momentsSpecSmoothCorr)

dbzLims=[-50,20];
velLims=[-18,18];
widthLims=[0,3];
skewLims=[-3,3];
kurtLims=[-8,8];
lrwidthLims=[0,15];
lslopeLims=[0,25];
rslopeLims=[-25,0];

dbzIn=momentsTime.dbz;
dbzIn(dbzIn<dbzLims(1))=dbzLims(1);
dbzIn(dbzIn>dbzLims(2))=dbzLims(2);
dbzScaled=1/(dbzLims(2)-dbzLims(1)).*(dbzIn-dbzLims(1));

velIn=momentsSpecSmoothCorr.velRaw;
velIn(velIn<velLims(1))=velLims(1);
velIn(velIn>velLims(2))=velLims(2);
velScaled=1/(velLims(2)-velLims(1)).*(velIn-velLims(1));

widthIn=momentsSpecSmoothCorr.width;
widthIn(widthIn<widthLims(1))=widthLims(1);
widthIn(widthIn>widthLims(2))=widthLims(2);
widthScaled=1/(widthLims(2)-widthLims(1)).*(widthIn-widthLims(1));

skewIn=momentsSpecSmoothCorr.skew;
skewIn(skewIn<skewLims(1))=skewLims(1);
skewIn(skewIn>skewLims(2))=skewLims(2);
skewScaled=1/(skewLims(2)-skewLims(1)).*(skewIn-skewLims(1));

kurtIn=momentsSpecSmoothCorr.kurt;
kurtIn(kurtIn<kurtLims(1))=kurtLims(1);
kurtIn(kurtIn>kurtLims(2))=kurtLims(2);
kurtScaled=1/(kurtLims(2)-kurtLims(1)).*(kurtIn-kurtLims(1));

lrwidthIn=momentsSpecSmoothCorr.lrwidth;
lrwidthIn(lrwidthIn<lrwidthLims(1))=lrwidthLims(1);
lrwidthIn(lrwidthIn>lrwidthLims(2))=lrwidthLims(2);
lrwidthScaled=1/(lrwidthLims(2)-lrwidthLims(1)).*(lrwidthIn-lrwidthLims(1));

lslopeIn=momentsSpecSmoothCorr.lslope;
lslopeIn(lslopeIn<lslopeLims(1))=lslopeLims(1);
lslopeIn(lslopeIn>lslopeLims(2))=lslopeLims(2);
lslopeScaled=1/(lslopeLims(2)-lslopeLims(1)).*(lslopeIn-lslopeLims(1));

rslopeIn=momentsSpecSmoothCorr.rslope;
rslopeIn(rslopeIn<rslopeLims(1))=rslopeLims(1);
rslopeIn(rslopeIn>rslopeLims(2))=rslopeLims(2);
rslopeScaled=1/(rslopeLims(2)-rslopeLims(1)).*(rslopeIn-rslopeLims(1));

stackVars=cat(3,dbzScaled,velScaled,widthScaled,skewScaled,kurtScaled,lrwidthScaled,lslopeScaled,rslopeScaled);

kmeansIn=reshape(stackVars,size(momentsSpecSmoothCorr.velRaw,1)*size(momentsSpecSmoothCorr.velRaw,2),[],1);

kInds=kmeans(kmeansIn,10,'MaxIter',1000);
micLabels=reshape(kInds,size(momentsSpecSmoothCorr.velRaw,1),size(momentsSpecSmoothCorr.velRaw,2));

end