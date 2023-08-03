function minInds=findMinInds(powerSpecSmooth,maskSpec)
minInds=nan(size(powerSpecSmooth,1),1);

nonNanInds=find(maskSpec==1);
minIndsNonNan=nan(length(nonNanInds)+1,1);
minIndsNonNan(1)=round(size(powerSpecSmooth,2)/2);

for ii=1:length(nonNanInds)
    powerGate=powerSpecSmooth(nonNanInds(ii),:);

    minPower=min(powerGate);
    minSpecInds=find(powerGate==minPower);

    diffMins=minSpecInds-minIndsNonNan(ii);
    [~,minDiffInd]=min(abs(diffMins));
    minIndsNonNan(ii+1)=minSpecInds(minDiffInd);
end
    
minIndsNonNan(1)=[];
minInds(nonNanInds)=minIndsNonNan;

end