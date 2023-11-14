function [powerSmoothOut,powerRMnoiseOut,locsMinOut,locsMaxOut]=findPeaksValleysSpec_test(powerAdj,velAdj,sampleNum)
% Find maxima and minima in spectra

powerLarge=repmat(powerAdj,1,3);
% powerLarge=nan(size(powerAdj,1),size(powerAdj,2)*3);
% powerLarge(:,sampleNum+1:2*sampleNum)=powerAdj;
% powerLarge=fillmissing(powerLarge,'linear',2,'EndValues',-61.301);

% powerMed=movmedian(powerLarge,51,2);
% rmOut=abs(powerMed-powerLarge);
% 
% powerMedRM=powerLarge;
% powerMedRM(rmOut>5)=nan;
% 
% powerMedRM=fillmissing(powerMedRM,'linear',2,'EndValues','nearest');
% powerSmoothBig=movmedian(powerMedRM,51,2);
% powerSmoothBig=movmedian(powerSmoothBig,21,2);

sampleNumOdd=round(sampleNum/6);
if ~mod(sampleNumOdd,2)
    sampleNumOdd=sampleNumOdd-1;
end

velMin=min(velAdj,[],2);
velMax=max(velAdj,[],2);
velDiff=velAdj(1,2)-velAdj(1,1);
velAdjLarge=cat(2,velAdj-(velMax-velMin)-velDiff,velAdj,velAdj+(velMax-velMin)+velDiff);

powerSmoothOut=nan(size(powerAdj));
powerRMnoiseOut=nan(size(powerAdj));

for ii=1:size(powerAdj,1)

    powerOrig=powerLarge(ii,:);
    regrFit=polyfit(velAdjLarge(ii,:),powerOrig,101);
    powerSmooth=polyval(regrFit,velAdjLarge(ii,:));
    
    % powerSmooth=powerSmoothBig(ii,:);
    % powerOrig=powerLarge(ii,:);

    %plot(powerOrig)
    %hold on
    %plot(powerMedRM(ii,:));
    %plot(powerSmooth,'-r','LineWidth',2)
    %hold off
    powerSmooth=powerSmooth(sampleNum+1:2*sampleNum);
    
    [locsMax,prom]=islocalmax(powerSmooth);
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
    
    %noisePower=powerSmooth;
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
    noisePerc=prctile(noisePower,50);

    powerRMnoise(powerSmooth<noisePerc)=nan;
    %powerRMnoise(powerSmooth<=max(noisePower))=nan;

    lineMask=~isnan(powerRMnoise);
    lineMask=bwareaopen(lineMask,round(sampleNum/6));

    powerRMnoise(lineMask==0)=nan;

    spread=max(powerRMnoise,[],'omitnan')-min(powerRMnoise,[],'omitnan');
    if spread<3
        powerRMnoise(:)=nan;
    end

    plot(powerRMnoise,'-g');
    hold off

    powerSmoothOut(ii,:)=powerSmooth;
    powerRMnoiseOut(ii,:)=powerRMnoise;
    locsMinOut.(['min',num2str(ii)])=locsMin;
    locsMaxOut.(['max',num2str(ii)])=locsMax;

end

end