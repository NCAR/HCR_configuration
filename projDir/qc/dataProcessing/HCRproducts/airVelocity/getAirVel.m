function [airVel,traceRefl]=getAirVel(powerAdj,phaseAdj,sampleNum,lambda,prt,range,dbz1km)
% Get air velocity from spectral data

airVel=nan(size(powerAdj,1),1);
traceRefl=nan(size(powerAdj,1),1);

powerLarge=cat(2,powerAdj,powerAdj,powerAdj);

powerMed=movmedian(powerLarge,round(sampleNum/6),2);
% powerMedSmooth=movmedian(powerMed,round(sampleNum/5),2);
% 
% powerMin=min(powerMedSmooth(:,sampleNum+1:2*sampleNum),[],2);
% 
% powerFilt=powerMed;
% powerFilt(powerFilt<powerMin)=nan;
% 
% powerMask=~isnan(powerFilt);
% 

powerMed=powerMed(:,sampleNum+1:2*sampleNum);
% for ii=1:length(airVel)
%     powerLine=powerMed(ii,:);
%     [pks,locs]=findpeaks(powerLine,'NPeaks',2,'MinPeakHeight',5);
% end

for ii=1:length(airVel)

    % Find maxima and minima
    thisPhase=phaseAdj(ii,:);
    powerLine=powerMed(ii,:);

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

%             if isempty(leftInd) | isempty(rightInd)
%                 if isempty(leftInd)
%                     leftInd=rightInd-round(sampleNum/10);
%                     rightInd=min([length(noisePower),rightInd+round(sampleNum/10)]);
%                 elseif isempty(rightInd)
%                     rightInd=leftInd+round(sampleNum/10);
%                     leftInd=max([1,leftInd-round(sampleNum/10)]);
%                     end
%                     noisePower(1:leftInd)=nan;
%             else
                noisePower(leftInd:rightInd)=nan;
            %end
        end
    end

    noiseMed=median(noisePower,'omitnan');
    noiseStd=std(noisePower,'omitnan');

    powerFilt=powerLine;
    powerFilt(powerLine<noiseMed+noiseStd)=nan;

    lineMask=~isnan(powerFilt);
    lineMask=bwareaopen(lineMask,round(sampleNum/6));

    powerFilt(lineMask==0)=nan;

    spread=max(powerFilt,[],'omitnan')-min(powerFilt,[],'omitnan');
    if spread<1
        powerFilt(:)=nan;
    end

    airInd=min(find(~isnan(powerFilt)));
    
    % Calculate air velocity and racer reflectivity
    if ~isnan(airInd)
        airVel(ii)=lambda/(4*pi*prt)*phaseAdj(ii,airInd);

        powerLin=10^(powerFilt(airInd)/10);
        noiseLin=10^((noiseMed+noiseStd)/10);
        snrLin=(powerLin-noiseLin)/noiseLin;
        snrLin(snrLin<0)=nan;
        snr=10*log10(snrLin);

        % DBZ
        range(range<0)=nan;
        traceRefl(ii)=snr+20*log10(range(ii)./1000)+dbz1km;
    end

    % Plot
%     plot(thisPhase,powerAdj(ii,:));
%     hold on
%     plot(thisPhase,powerLine);
%     plot(thisPhase,powerFilt,'LineWidth',1.5);
%     scatter(thisPhase(locsMax),powerLine(locsMax),'filled');
%     scatter(thisPhase(locsMin),powerLine(locsMin),'filled');
%     hold off
%     ylim([-70,0]);
end

end