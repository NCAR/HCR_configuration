function [powerSmoothOut,powerRMnoiseOut]=findPeaksValleysSpec_test(powerAdj,velAdj,sampleNum)
% Find maxima and minima in spectra

powerSmoothOut=nan(size(powerAdj));
powerRMnoiseOut=nan(size(powerAdj));

for ii=1:size(powerAdj,1)

    powerOrig=powerAdj(ii,:);
    B=dop(velAdj(ii,:)',25);
    powerSmooth=B*B'*powerOrig';
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

    plot(powerAdj(ii,:))
    hold on
    plot(powerSmooth,'-r','linewidth',2);
    scatter(locsMin,powerSmooth(locsMin),'filled')
    scatter(locsMax,powerSmooth(locsMax),'filled')
    
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
   
    lineMask=~isnan(powerRMnoise);
    lineMask=bwareaopen(lineMask,round(sampleNum/20));

    powerRMnoise(lineMask==0)=nan;

    powerWithNoise=powerSmooth;
    powerWithNoise(~isnan(powerRMnoise))=nan;
    powerRMnoise(powerSmooth<=max(powerWithNoise,[],'omitmissing'))=nan;

    lineMask=~isnan(powerRMnoise);
    lineMask=bwareaopen(lineMask,round(sampleNum/20));

    powerRMnoise(lineMask==0)=nan;

    noisePower2=powerAdj(ii,:);
    noisePower2(~isnan(powerRMnoise))=nan;
    noisePerc2=prctile(noisePower2,60);

    powerRMnoise(powerSmooth<noisePerc2)=nan;

    lineMask=~isnan(powerRMnoise);
    lineMask=bwareaopen(lineMask,round(sampleNum/20));

    powerRMnoise(lineMask==0)=nan;

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

    for kk=1:length(startLine)
        lineTest=powerRMnoise(startLine(kk):endLine(kk));
        spread=max(lineTest,[],'omitnan')-min(lineTest,[],'omitnan');
        if spread<1
            powerRMnoise(startLine(kk):endLine(kk))=nan;
        end
    end

    plot(powerRMnoise,'-g');
    hold off

    powerSmoothOut(ii,:)=powerSmooth;
    powerRMnoiseOut(ii,:)=powerRMnoise;

end
end