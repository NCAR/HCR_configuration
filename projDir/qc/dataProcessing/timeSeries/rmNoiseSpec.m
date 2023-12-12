function [powerSmoothOut,powerRMnoiseOut]=rmNoiseSpec(powerAdj)
sampleNum=size(powerAdj,2);

% Find maxima and minima in spectra

powerSmoothOut=nan(size(powerAdj));
powerRMnoiseOut=nan(size(powerAdj));

loopInds=find(any(~isnan(powerAdj),2));

x=linspace(-1,1,sampleNum);
[~,powerSmoothAll,~]= forsythe_polyfit(x,powerAdj',25);

for aa=1:size(loopInds,1)
    ii=loopInds(aa);

    powerSmooth=powerSmoothAll(:,ii);
    powerSmooth=powerSmooth';
    powerSmooth(1:round(sampleNum/100))=nan;
    powerSmooth(sampleNum-round(sampleNum/100)+1:end)=nan;

    [locsMax,prom]=islocalmax(powerSmooth,'MinSeparation',round(sampleNum/6),'MinProminence',1);
    locsMax=find(locsMax==1);
    prom=prom(locsMax);
    locsMin=islocalmin(powerSmooth);
    locsMin=find(locsMin==1);

    if isempty(locsMin)
        [~,locsMin]=min(powerSmooth);
    end

    locsMin=cat(2,1+round(sampleNum/20),locsMin,length(powerSmooth)-round(sampleNum/20));

    noisePower=powerAdj(ii,:);

    if ~isempty(locsMax)
        peakPower=powerSmooth(locsMax);
        valPower=powerSmooth(locsMin);

        for jj=1:length(locsMax)
            thisMax=locsMax(jj);
            promThis=prom(jj);
            peakPowerThis=peakPower(jj);

            powerDiff=peakPowerThis-valPower;
            locsMinThis=locsMin;
            locsMinThis(powerDiff<promThis-promThis/3)=[];

            diffMaxMin=thisMax-locsMinThis;

            leftMin=max(find(diffMaxMin>0));
            leftInd=locsMinThis(leftMin);
            rightMin=min(find(diffMaxMin<0));
            rightInd=locsMinThis(rightMin);

            if isempty(leftInd)
                leftInd=1;
            end
            if isempty(rightInd)
                rightInd=length(noisePower);
            end

            noisePower(leftInd:rightInd)=nan;
        end
    end

    powerRMnoise=powerSmooth;
    noisePerc=prctile(noisePower,60);

    powerRMnoise(powerSmooth<noisePerc)=nan;

    lineMask=double(~isnan(powerRMnoise));
    lineMask(lineMask==0)=nan;
    lineMask=movmedian(lineMask,round(sampleNum/20),'includemissing');
    lineMask=movmedian(lineMask,round(sampleNum/20),'omitmissing');

    powerRMnoise(isnan(lineMask))=nan;

    powerWithNoise=powerSmooth;
    powerWithNoise(~isnan(powerRMnoise))=nan;
    powerRMnoise(powerSmooth<=max(powerWithNoise,[],'omitmissing'))=nan;

    lineMask=double(~isnan(powerRMnoise));
    lineMask(lineMask==0)=nan;
    lineMask=movmedian(lineMask,round(sampleNum/20),'includemissing');
    lineMask=movmedian(lineMask,round(sampleNum/20),'omitmissing');

    powerRMnoise(isnan(lineMask))=nan;

    noisePower2=powerAdj(ii,:);
    noisePower2(~isnan(powerRMnoise))=nan;
    noisePerc2=prctile(noisePower2,60);

    powerRMnoise(powerSmooth<noisePerc2)=nan;

    lineMask=double(~isnan(powerRMnoise));
    lineMask(lineMask==0)=nan;
    lineMask=movmedian(lineMask,round(sampleNum/20),'includemissing');
    lineMask=movmedian(lineMask,round(sampleNum/20),'omitmissing');

    powerRMnoise(isnan(lineMask))=nan;

    lineMask=~isnan(powerRMnoise);
    diffLine=diff(lineMask);

    startLine=find(diffLine==1);
    endLine=find(diffLine==-1);

    if ~isempty(startLine) | ~isempty(endLine)
        if isempty(startLine) & ~isempty(endLine)
            startLine=1;
        end
        if ~isempty(startLine) & isempty(endLine)
            endLine=length(powerRMnoise);
        end
        if startLine(1)>endLine(1)
            startLine=[1,startLine];
        end
        if startLine(end)>endLine(end)
            endLine=[endLine,length(powerRMnoise)];
        end
    end

    rmLine=zeros(size(startLine));
    lineMax=nan(size(startLine));
    for kk=1:length(startLine)
        lineTest=powerRMnoise(startLine(kk):endLine(kk));
        spread=max(lineTest,[],'omitnan')-min(lineTest,[],'omitnan');
        if spread<2
            rmLine(kk)=1;
        end
        lineMax(kk)=max(lineTest);
    end

    rmLine(lineMax~=max(lineMax))=1;
    rmInds=find(rmLine==1);
    for kk=1:length(rmInds)
        powerRMnoise(startLine(rmInds(kk)):endLine(rmInds(kk)))=nan;
        powerSmooth(startLine(rmInds(kk)):endLine(rmInds(kk)))=nan;
    end

    powerSmoothOut(ii,:)=powerSmooth;
    powerRMnoiseOut(ii,:)=powerRMnoise;

end
end