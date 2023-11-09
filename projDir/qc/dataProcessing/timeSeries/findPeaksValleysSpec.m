function [locsMin,locsMax]=findPeaksValleysSpec(powerAdj,velAdj,sampleNum,plotRangeKM,plotInd,range)
% Find maxima and minima in spectra

% Find index of specified range
rangeInd=min(find((range./1000)>=plotRangeKM(plotInd)));

powerLarge=repmat(powerAdj,1,3);

sampleNumOdd=round(sampleNum/6);
if ~mod(sampleNumOdd,2)
    sampleNumOdd=sampleNumOdd-1;
end

%powerSmooth=sgolayfilt(powerLarge,3,sampleNumOdd,[],2);
%powerSmooth=powerSmooth(:,sampleNum+1:2*sampleNum);
x=repmat(velAdj,1,3);
regrFit=polyfit(x,powerLarge,21);
powerSmooth=polyval(regrFit,x);

for ii=1:size(powerAdj,1)

    % Find maxima and minima
    thisVel=velAdj(ii,:);
    if max(~isnan(thisVel))==0
        continue
    end
    %powerLine=powerSmooth(ii,:);
    powerOrig=powerLarge(ii,:);
    regrFit=polyfit(x(ii,:),powerOrig,25);
    powerLine=polyval(regrFit,x(ii,:));
    powerLine=powerLine(sampleNum+1:2*sampleNum);
    
    [locsMax,prom]=islocalmax(powerLine,'MinSeparation',sampleNum/10,'MinProminence',1);
    locsMax=find(locsMax==1);
    prom=prom(locsMax);
    locsMin=islocalmin(powerLine,'MinProminence',0.5);
    locsMin=find(locsMin==1);

    if isempty(locsMin)
        [~,locsMin]=min(powerLine);
    end

    locsMin=cat(2,1+round(sampleNum/20),locsMin,length(powerLine)-round(sampleNum/20));

    noisePower=powerLine;

    if ~isempty(locsMax)
        peakPower=powerLine(locsMax);
        valPower=powerLine(locsMin);

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

    powerFilt=powerLine;
    noisePerc=prctile(noisePower,90);

    powerFilt(powerLine<noisePerc)=nan;

    lineMask=~isnan(powerFilt);
    lineMask=bwareaopen(lineMask,round(sampleNum/6));

    powerFilt(lineMask==0)=nan;

    spread=max(powerFilt,[],'omitnan')-min(powerFilt,[],'omitnan');
    if spread<3
        powerFilt(:)=nan;
    end

    % Plot
    if ii==rangeInd
        figure
        plot(thisVel,powerAdj(ii,:));
        hold on
        plot(thisVel,powerLine);
        plot(thisVel,powerFilt,'LineWidth',1.5);
        scatter(thisVel(locsMax),powerLine(locsMax),'filled','MarkerFaceColor','m');
        scatter(thisVel(locsMin),powerLine(locsMin),'filled','MarkerFaceColor','c');
        hold off
        xlim([thisVel(1),thisVel(end)]);
        ylim([-80,10]);
        stopHere=1;
    end
end

end