function [airVel,traceRefl]=getAirVel(powerAdj,phaseAdj,elev,sampleNum,lambda,prt,range,dbz1km,noiseIn)
% Get air velocity from spectral data

airVel=nan(size(powerAdj,1),1);
traceRefl=nan(size(powerAdj,1),1);

powerLarge=repmat(powerAdj,1,3);

powerMed=movmedian(powerLarge,round(sampleNum/6),2);
powerMed=powerMed(:,sampleNum+1:2*sampleNum);

for ii=1:length(airVel)

    % Find maxima and minima
    thisPhase=phaseAdj(ii,:);
    if max(~isnan(thisPhase))==0
        continue
    end
    powerLine=powerMed(ii,:);
    powerRaw=powerAdj(ii,:);

    [locsMax,prom]=islocalmax(powerLine,'MinSeparation',sampleNum/5,'MinProminence',2);
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

        leftMin=[];
        rightMin=[];
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
    if spread<1
        powerFilt(:)=nan;
    end

    if elev<0
        airInd=min(find(~isnan(powerFilt)));
    else
        airInd=max(find(~isnan(powerFilt)));
    end
    
    % Calculate air velocity and tracer reflectivity
    if ~isempty(airInd)
        airVel(ii)=lambda/(4*pi*prt)*thisPhase(airInd);

        % SNR
        powerLin=10^(powerFilt(airInd)/10);
        %noiseLin=10^((noiseMed+noiseStd)/10);
        noiseLin=10^(noiseIn/10);
        snrLin=(powerLin-noiseLin)/noiseLin;
        snrLin(snrLin<0)=nan;
        snr=10*log10(snrLin);

        % DBZ
        range(range<0)=nan;
        traceRefl(ii)=snr+20*log10(range(ii)./1000)+dbz1km;
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