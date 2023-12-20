function [powerAdj,specVelAdj]=adjSpecBoundsV(specDB,velIn,data)

%% Filter

sampleNum=length(data.time);

duplicateSpec=5;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec+2);
powerSpecSmoothLarge=movmedian(powerSpecLarge,round(size(specDB,2))/5,2);
powerSpecSmoothLarge=powerSpecSmoothLarge(:,sampleNum+1:end-sampleNum);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

powerSpecLarge=powerSpecLarge(:,sampleNum+1:end-sampleNum);

powerSpecSmooth=powerSpecSmoothLarge(:,floor(duplicateSpec/2)*sampleNum+1:end-floor(duplicateSpec/2)*sampleNum);

%% Find maximum

maxInds=findMaxIndsSpec(powerSpecSmooth,velIn);

maxInds=maxInds+floor(duplicateSpec/2)*sampleNum;

%% Build adjusted spectrum

leftAdd=floor(sampleNum/2);
if mod(sampleNum,2)==0
    rightAdd=leftAdd-1;
else
    rightAdd=leftAdd;
end

powerAdj=nan(size(specDB));
specVelAdj=nan(size(specDB));

for kk=1:size(specDB,1)
    try
        powerAdj(kk,:)=powerSpecLarge(kk,maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
        specVelAdj(kk,:)=velSpecLarge(maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
    end
end

powerAdj(isnan(velIn),:)=nan;

end