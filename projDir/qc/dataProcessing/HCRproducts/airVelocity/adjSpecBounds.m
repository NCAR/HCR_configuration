function [powerAdj,specVelAdj]=adjSpecBounds(powerSpec,velDeAlias,sampleNum)

%% Filter

duplicateSpec=5;

% Add spectra side by side
powerSpecLarge=repmat(powerSpec,1,duplicateSpec+2);
powerSpecSmoothLarge=movmedian(powerSpecLarge,round(size(powerSpec,2))/5,2);
powerSpecSmoothLarge=powerSpecSmoothLarge(:,sampleNum+1:end-sampleNum);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1);

powerSpecLarge=powerSpecLarge(:,sampleNum+1:end-sampleNum);

powerSpecSmooth=powerSpecSmoothLarge(:,floor(duplicateSpec/2)*sampleNum+1:end-floor(duplicateSpec/2)*sampleNum);

%% Find maximum

maxInds=findMaxIndsSpec(powerSpecSmooth,velDeAlias);

maxInds=maxInds+floor(duplicateSpec/2)*sampleNum;

%% Build adjusted spectrum

leftAdd=floor(sampleNum/2);
if mod(sampleNum,2)==0
    rightAdd=leftAdd-1;
else
    rightAdd=leftAdd;
end

powerAdj=nan(size(powerSpec));
specVelAdj=nan(size(powerSpec));

for kk=1:size(powerSpec,1)
    try
        powerAdj(kk,:)=powerSpecLarge(kk,maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
        specVelAdj(kk,:)=velSpecLarge(maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
    end
end

powerAdj(isnan(velDeAlias),:)=nan;

end