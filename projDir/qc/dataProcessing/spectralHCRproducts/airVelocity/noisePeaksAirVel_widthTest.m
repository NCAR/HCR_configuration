function [powerRMnoiseAvRM,powerRMnoiseAvRMS,velOut,velOutS]=noisePeaksAirVel_widthTest(specDB,velIn,data,widthC)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries

sampleNum=length(data.time);
halfSN=round(sampleNum/2);

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

largeInds=1:length(velSpecLarge);

powerRMnoiseAvRM=nan(size(specDB));
powerRMnoiseAvRMS=nan(size(specDB));
velOut=nan(size(specDB));
velOutS=nan(size(specDB));

meanNoiseAll=nan(size(specDB,1),1);
noiseThreshAll=meanNoiseAll;

%% Remove noise
% Moving average
meanOverPoints=3; % Average over this number of points
secondMean=round(sampleNum/3);
movAv=movmedian(powerSpecLarge,meanOverPoints,2);
movAv2=movmedian(movAv,secondMean,2);

movAv(:,1:round(sampleNum/10))=nan;
movAv(:,end-round(sampleNum/10):end)=nan;

movAv2(:,1:round(sampleNum/3))=nan;
movAv2(:,end-round(sampleNum/3):end)=nan;

loopInds=find(any(~isnan(specDB),2));

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    thisMov=movAv(ii,:);

    % Get one whole spectrum
    [~,minIndTest]=min(movAv2(ii,:));

    testPow=thisMov(minIndTest:minIndTest+sampleNum-1);
    testVel=velSpecLarge(minIndTest:minIndTest+sampleNum-1);

    % Find noise floor and noise threshold    
    [noiseThreshAll(ii),meanNoiseAll(ii)]=findNoiseThresh(testPow,meanOverPoints);

    % Correct for aircraft width
    [powWidthCorr,powFilt]=aircraftWidthCorr(testPow,widthC,testVel,noiseThreshAll(ii),sampleNum);
    powWClarge=repmat(powWidthCorr,1,duplicateSpec+2);
    powShift=powWClarge(sampleNum*2-minIndTest:sampleNum*2-minIndTest+length(velSpecLarge)-1);

    powSlarge=repmat(powFilt,1,duplicateSpec+2);
    powShiftS=powSlarge(sampleNum*2-minIndTest:sampleNum*2-minIndTest+length(velSpecLarge)-1);

    % Remove noise below threshold for corrected
    thisMovRM=powShift;
    thisMovRM(thisMovRM<noiseThreshAll(ii))=nan;
    thisMask1=~isnan(thisMovRM);
    thisMovRM(thisMask1==0)=nan;

    % Find indices of maxima
    maxMov=max(thisMov,[],2,'omitmissing');
    maxInds=find(thisMov==maxMov);

    % Find correct velocity
    velAtInds=velSpecLarge(maxInds);
    [~,velDiffMin]=min(abs(velAtInds-velIn(ii)));
    velDiffMinInd=maxInds(velDiffMin);

    smallerInds=largeInds(velDiffMinInd-halfSN:velDiffMinInd+halfSN);
    smallerMask=thisMask1(smallerInds);
    smallerMov=thisMov(smallerInds);
    maxInd1=find(smallerMov==maxMov);

    specRegs=bwconncomp(smallerMask);

    for jj=1:specRegs.NumObjects
        if ~ismember(maxInd1,specRegs.PixelIdxList{jj})
            smallerMask(specRegs.PixelIdxList{jj})=0;
        end
    end

    maskOneInd=find(smallerMask,1,'first');
    finalInds=smallerInds(maskOneInd):smallerInds(maskOneInd)+sampleNum-1;

     if ~isempty(finalInds)
        largeMask=zeros(size(largeInds));
        largeMask(smallerInds)=smallerMask;

        thisMovRM(largeMask==0)=nan;
        thisMovRM=thisMovRM(finalInds);
       
        thisVel=velSpecLarge(finalInds);

        powerRMnoiseAvRM(ii,:)=thisMovRM;
        velOut(ii,:)=thisVel;
    else
        continue
    end

    % Remove noise below threshold for smooth
    thisMovRMS=powShiftS;
    thisMovRMS(thisMovRMS<noiseThreshAll(ii))=nan;
    thisMask1=~isnan(thisMovRMS);
    thisMovRMS(thisMask1==0)=nan;

    % Find correct velocity
    smallerMask=thisMask1(smallerInds);
    smallerMov=thisMov(smallerInds);
    maxInd1=find(smallerMov==maxMov);

    specRegs=bwconncomp(smallerMask);

    for jj=1:specRegs.NumObjects
        if ~ismember(maxInd1,specRegs.PixelIdxList{jj})
            smallerMask(specRegs.PixelIdxList{jj})=0;
        end
    end

    maskOneInd=find(smallerMask,1,'first');
    finalInds=smallerInds(maskOneInd):smallerInds(maskOneInd)+sampleNum-1;

    if ~isempty(finalInds)
        largeMask=zeros(size(largeInds));
        largeMask(smallerInds)=smallerMask;

        thisMovRMS(largeMask==0)=nan;
        thisMovRMS=thisMovRMS(finalInds);

        thisVelS=velSpecLarge(finalInds);

        velOutS(ii,:)=thisVelS;
       
        powerRMnoiseAvRMS(ii,:)=thisMovRMS;
    else
        continue
    end
end
end