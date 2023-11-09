function [powerSmoothOut,powerRMnoiseOut,locsMinOut,locsMaxOut]=findPeaksValleysSpec(powerAdj,velAdj,sampleNum)
% Find maxima and minima in spectra

powerLarge=repmat(powerAdj,1,3);

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

    % Find maxima and minima
    % thisVel=velAdj(ii,:);
    % if max(~isnan(thisVel))==0
    %     continue
    % end
    
    powerOrig=powerLarge(ii,:);
    regrFit=polyfit(velAdjLarge(ii,:),powerOrig,25);
    powerSmooth=polyval(regrFit,velAdjLarge(ii,:));
    powerSmooth=powerSmooth(sampleNum+1:2*sampleNum);
    
    [locsMax,prom]=islocalmax(powerSmooth,'MinSeparation',sampleNum/10,'MinProminence',1);
    locsMax=find(locsMax==1);
    prom=prom(locsMax);
    locsMin=islocalmin(powerSmooth,'MinProminence',0.5);
    locsMin=find(locsMin==1);

    if isempty(locsMin)
        [~,locsMin]=min(powerSmooth);
    end

    locsMin=cat(2,1+round(sampleNum/20),locsMin,length(powerSmooth)-round(sampleNum/20));

    noisePower=powerSmooth;

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

            noisePower(leftInd:rightInd)=nan;
         end
    end

    powerRMnoise=powerSmooth;
    noisePerc=prctile(noisePower,90);

    powerRMnoise(powerSmooth<noisePerc)=nan;

    lineMask=~isnan(powerRMnoise);
    lineMask=bwareaopen(lineMask,round(sampleNum/6));

    powerRMnoise(lineMask==0)=nan;

    spread=max(powerRMnoise,[],'omitnan')-min(powerRMnoise,[],'omitnan');
    if spread<3
        powerRMnoise(:)=nan;
    end

    powerSmoothOut(ii,:)=powerSmooth;
    powerRMnoiseOut(ii,:)=powerRMnoise;
    locsMinOut.(['min',num2str(ii)])=locsMin;
    locsMaxOut.(['max',num2str(ii)])=locsMax;

    % Plot
    % if ii==rangeInd
    %     figure
    %     plot(thisVel,powerAdj(ii,:));
    %     hold on
    %     plot(thisVel,powerSmooth);
    %     plot(thisVel,powerRMnoise,'LineWidth',1.5);
    %     scatter(thisVel(locsMax),powerSmooth(locsMax),'filled','MarkerFaceColor','m');
    %     scatter(thisVel(locsMin),powerSmooth(locsMin),'filled','MarkerFaceColor','c');
    %     hold off
    %     xlim([thisVel(1),thisVel(end)]);
    %     ylim([-80,10]);
    %     stopHere=1;
    % end
end

end