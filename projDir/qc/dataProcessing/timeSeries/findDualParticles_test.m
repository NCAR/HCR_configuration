function dualParticles=findDualParticles_test(powerIn,specVelIn)
% Find maxima and minima in spectra
close all
dualParticles=[];
sampleNum=size(powerIn,2);

dataInds=find(any(~isnan(powerIn),2));

% Remove transmitter pulse
dataInds(dataInds<=12)=[];

for jj=1:length(dataInds)
    ii=dataInds(jj);

    powerOrig=powerIn(ii,:);

    % Find start and end of line segments
    lineMask=~isnan(powerOrig);
    diffLine=diff(lineMask);

    startLine=find(diffLine==1)+1;
    endLine=find(diffLine==-1);

    if isempty(startLine) & ~isempty(endLine)
        startLine=1;
    end
    if ~isempty(startLine) & isempty(endLine)
        endLine=length(powerOrig);
    end
    if startLine(1)>endLine(1)
        startLine=[1,startLine];
    end
    if startLine(end)>endLine(end)
        endLine=[endLine,length(powerOrig)];
    end

    % Find maxima and minima
    [locsMax,prom]=islocalmax(powerOrig,'MaxNumExtrema',2);
    locsMax=find(locsMax==1);
    prom=prom(locsMax);
    locsMin=islocalmin(powerOrig);
    locsMin=find(locsMin==1);

    locsMin=cat(2,startLine,endLine,locsMin);
    locsMin=unique(locsMin);

    maxEnd=ismember(locsMax,[startLine,endLine]);
    locsMax(maxEnd==1)=[];

    %% Find second peak if it doesn't exist

    % Find change points
    diffCurve=cat(2,nan,diff(powerOrig));
    [diffMax,promD]=islocalmax(diffCurve);
    diffMax=find(diffMax==1);
    promD=promD(diffMax);
    diffMin=islocalmin(diffCurve);
    diffMin=find(diffMin==1);

    % Find points to test
    testX=[];
    testY=[];
    addPeak=[];
    if length(locsMax)==1

        % Left side
        val1=diffMax(diffMax<locsMax);
        peak1=diffMin(diffMin<locsMax);

        for kk=1:length(peak1)
            thisPeak=peak1(kk);
            peakValDiff=thisPeak-val1;

            nextLeftInd=val1(max(find(peakValDiff>0)));
            nextRightInd=val1(min(find(peakValDiff<0)));
            if ~isempty(nextLeftInd) & ~isempty(nextRightInd)
                testX=cat(1,testX,[nextLeftInd,nextRightInd]);
            end
        end

        % Right side
        val2=diffMin(diffMin>locsMax);
        peak2=diffMax(diffMax>locsMax);

        for kk=1:length(peak2)
            thisPeak=peak2(kk);
            peakValDiff=thisPeak-val2;

            nextLeftInd=val2(max(find(peakValDiff>0)));
            nextRightInd=val2(min(find(peakValDiff<0)));
            if ~isempty(nextLeftInd) & ~isempty(nextRightInd)
                testX=cat(1,testX,[nextLeftInd,nextRightInd]);
            end
        end

        testPeak=[];
        for kk=1:size(testX,1)
            testXthis=testX(kk,:);
            testYthis=powerOrig(testXthis);
            testY=cat(1,testY,testYthis);
            fitLineMod=polyfit(testXthis,testYthis,1);
            fitLine=polyval(fitLineMod,testX(kk,1):testX(kk,2));
            powerMinusLine=powerOrig(testXthis(1):testXthis(2))-fitLine;
            [maxPeak,maxPeakInd]=max(powerMinusLine);
            maxPeakInd=maxPeakInd+testXthis(1);
            if maxPeak>3
                testPeak=cat(1,testPeak,[maxPeak,maxPeakInd]);
            end
        end
        if ~isempty(testPeak)
            [~,maxAddInd]=max(testPeak(:,1));
            goodPeak=testPeak(maxAddInd,2);

            if goodPeak<locsMax
                peakValDiff=goodPeak-peak1;
                addPeak=peak1(min(find(peakValDiff<0)));
            else
                peakValDiff=goodPeak-peak2;
                addPeak=peak2(max(find(peakValDiff>0)));
            end
        end
    end

    plot(powerOrig,'-b','linewidth',2);
    hold on
    scatter(locsMin,powerOrig(locsMin),'filled','MarkerFaceColor','r')
    scatter(locsMax,powerOrig(locsMax),80,'filled','MarkerFaceColor','g')
    scatter(diffMin,powerOrig(diffMin),'filled','MarkerFaceColor','y')
    scatter(diffMax,powerOrig(diffMax),'filled','MarkerFaceColor','k')
    if ~isempty(testX)
        for mm=1:size(testX,1)
            scatter(testX(mm,:),testY(mm,:),'filled','MarkerFaceColor','m')
        end
    end
    if ~isempty(addPeak)
        scatter(addPeak,powerOrig(addPeak),80,'filled','MarkerFaceColor','c');
    end

    yyaxis right
    plot(diffCurve,'-c','linewidth',1.5);
    hold on
    scatter(diffMin,diffCurve(diffMin),'filled','MarkerFaceColor','y')
    scatter(diffMax,diffCurve(diffMax),'filled','MarkerFaceColor','k')
    ylim([-0.5,0.5]);
    grid on

    cla
    yyaxis left
    cla
end
end