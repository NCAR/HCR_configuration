function [err]=noisePeaks_smoothingTestTime(specDB,specDB1,specDB2,data,err,figdir)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries

sampleNum=length(data.time);
sampleNum1=floor(sampleNum/2);

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);
powerSpecLarge1=repmat(specDB1,1,duplicateSpec);
powerSpecLarge2=repmat(specDB2,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

velSpecLarge1=-duplicateSpec*pi:2*pi/(sampleNum1):duplicateSpec*pi;
velSpecLarge1=velSpecLarge1(1:end-1).*mode(data.lambda)./(4*pi.*mode(data.prt*2));

%% Remove noise
% Moving average
movAv=movmedian(powerSpecLarge,3,2);
secondMean=round(sampleNum/3);
movAvPlus=movmedian(movAv,secondMean,2);

movAvPlus(:,1:round(sampleNum/3))=nan;
movAvPlus(:,end-round(sampleNum/3):end)=nan;

movAv1=movmedian(powerSpecLarge1,3,2);
secondMean1=round(sampleNum1/3);
movAvPlus1=movmedian(movAv1,secondMean1,2);

movAvPlus1(:,1:round(sampleNum1/3))=nan;
movAvPlus1(:,end-round(sampleNum1/3):end)=nan;

loopInds=find(any(~isnan(specDB),2));

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    thisMov=powerSpecLarge(ii,:);
    thisMov1=powerSpecLarge1(ii,:);
    thisMov2=powerSpecLarge2(ii,:);

    % Get one whole spectrum
    [~,minIndTest]=min(movAvPlus(ii,:));
    [~,minIndTest1]=min(movAvPlus1(ii,:));

    testPow=thisMov(minIndTest:minIndTest+sampleNum-1);
    testVel=velSpecLarge(minIndTest:minIndTest+sampleNum-1);

    testPow1=thisMov1(minIndTest1:minIndTest1+sampleNum1-1);
    testPow2=thisMov2(minIndTest1:minIndTest1+sampleNum1-1);
    testVel1=velSpecLarge1(minIndTest1:minIndTest1+sampleNum1-1);

    % Correct for aircraft width
    [err]=smoothingTestTime(testPow,testVel,testPow1,testVel1,testPow2,err,figdir);
end
end