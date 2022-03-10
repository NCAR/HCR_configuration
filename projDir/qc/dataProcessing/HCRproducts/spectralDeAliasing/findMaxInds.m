function maxInds=findMaxInds(powerSpecSmooth,maskSpec)
maxInds=nan(size(powerSpecSmooth,1),1);

nonNanInds=find(maskSpec==1);
maxIndsNonNan=nan(length(nonNanInds)+1,1);
maxIndsNonNan(1)=round(size(powerSpecSmooth,2)/2);

for ii=1:length(nonNanInds)
    powerGate=powerSpecSmooth(nonNanInds(ii),:);

    maxPower=max(powerGate);
    maxSpecInds=find(powerGate==maxPower);

    diffMax=maxSpecInds-maxIndsNonNan(ii);
    [~,minDiffInd]=min(abs(diffMax));
    maxIndsNonNan(ii+1)=maxSpecInds(minDiffInd);
end
    
maxIndsNonNan(1)=[];
maxInds(nonNanInds)=maxIndsNonNan;

end