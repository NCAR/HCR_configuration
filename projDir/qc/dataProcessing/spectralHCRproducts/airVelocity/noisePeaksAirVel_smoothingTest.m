function [err]=noisePeaksAirVel_smoothingTest(specDB,velIn,data,widthC,err,figdir)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries

sampleNum=length(data.time);

if sampleNum~=987
    return
end

halfSN=round(sampleNum/2);

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

largeInds=1:length(velSpecLarge);

powerRMnoise=nan(size(specDB));
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

%% Create random split indices
ind1=[];
ind2=[];

for bb=1:10
    randP=randperm(100);
    ind1=cat(2,ind1,randP(1:2:100)+100*(bb-1));
    ind2=cat(2,ind2,randP(2:2:100)+100*(bb-1));
end

ind1=sort(ind1);
ind2=sort(ind2);

ind1(ind1>sampleNum-1)=[];
ind2(ind2>sampleNum-1)=[];

oneMinTwo=length(ind1)-length(ind2);
moveInds=abs(oneMinTwo/2);
if oneMinTwo<0
    ind1=cat(2,ind1,ind2(end-moveInds+1:end));
    ind2(end-moveInds+1:end)=[];
elseif oneMinTwo>0
    ind2=cat(2,ind2,ind1(end-moveInds+1:end));
    ind1(end-moveInds+1:end)=[];
end

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
    [err]=smoothingTest(testPow,widthC,testVel,noiseThreshAll(ii),sampleNum,err,ind1,ind2,figdir);
end
end