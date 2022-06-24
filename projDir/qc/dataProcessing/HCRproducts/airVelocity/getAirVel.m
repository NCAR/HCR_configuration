function [airVel,smallerVel,largerVel,maxVel,traceRefl]=getAirVel(powerAdj,phaseAdj,elev,sampleNum,lambda,prt,range,dbz1km,noiseIn)
% Get air velocity from spectral data

airVel=nan(size(powerAdj,1),1);
smallerVel=nan(size(powerAdj,1),1);
largerVel=nan(size(powerAdj,1),1);
maxVel=nan(size(powerAdj,1),1);
traceRefl=nan(size(powerAdj,1),1);

powerLarge=repmat(powerAdj,1,3);

%powerSmooth=movmedian(powerLarge,round(sampleNum/6),2);
powerSmooth=sgolayfilt(powerLarge,3,round(sampleNum/6),[],2);
powerSmooth=powerSmooth(:,sampleNum+1:2*sampleNum);

for ii=1:length(airVel)

    % Find maxima and minima
    thisPhase=phaseAdj(ii,:);
    if max(~isnan(thisPhase))==0
        continue
    end
    powerLine=powerSmooth(ii,:);
    powerRaw=powerAdj(ii,:);

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

    % Calculate left and right vel
    peakPowInd=cat(2,powerFilt(locsMax)',locsMax');
    peakPowInd(any(isnan(peakPowInd),2),:)=[];

    if ~isempty(peakPowInd)
        peakPowSort=sortrows(peakPowInd,'descend');

        maxVel(ii)=lambda/(4*pi*prt)*thisPhase(peakPowSort(1,2));
        if size(peakPowSort)>1
            twoPeakInds=peakPowSort(1:2,2);
            if elev<0
                smallerVel(ii)=lambda/(4*pi*prt)*thisPhase(min(twoPeakInds));
                largerVel(ii)=lambda/(4*pi*prt)*thisPhase(max(twoPeakInds));
            else
                smallerVel(ii)=lambda/(4*pi*prt)*thisPhase(max(twoPeakInds));
                largerVel(ii)=lambda/(4*pi*prt)*thisPhase(min(twoPeakInds));
            end
        end
    end

    % Calculate air velocity and tracer reflectivity

    if elev<0
        airInd=min(find(~isnan(powerFilt)));
    else
        airInd=max(find(~isnan(powerFilt)));
    end

    if ~isempty(airInd)
        airVel(ii)=lambda/(4*pi*prt)*thisPhase(airInd);

        % DBZ
        range(range<0)=nan;
        traceRefl(ii)=powerFilt(airInd)-noiseIn+20*log10(range(ii)./1000)+dbz1km;
    end

    % Plot
    plotYes=0;
    if plotYes
        plot(thisPhase,powerAdj(ii,:));
        hold on
        plot(thisPhase,powerLine);
        plot(thisPhase,powerFilt,'LineWidth',1.5);
        scatter(thisPhase(locsMax),powerLine(locsMax),'filled','MarkerFaceColor','m');
        scatter(thisPhase(locsMin),powerLine(locsMin),'filled','MarkerFaceColor','c');
        scatter(thisPhase(airInd),powerLine(airInd),'filled','MarkerFaceColor','r');
        hold off
        ylim([-80,10]);
        stopHere=1;
    end
end

end