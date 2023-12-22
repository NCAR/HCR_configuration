function momentsVelDual=findMultiVels(powerIn,specVelIn,powerRaw,momentsVelDual,nn)
% Find maxima and minima in spectra
velDual=nan(size(powerIn,1),35);

dataInds=find(any(~isnan(powerIn),2));

% Remove transmitter pulse
dataInds(dataInds<=12)=[];

locsMaxAll=islocalmax(powerIn,2);

diffCurve=cat(2,nan(size(powerIn,1),1),diff(powerIn,1,2));
diffMaxAll=islocalmax(diffCurve,2);
diffMinAll=islocalmin(diffCurve,2);

for jj=1:length(dataInds)
    ii=dataInds(jj);
    
    powerOrig=powerIn(ii,:);

    % Find start and end of line segments
    lineMask=~isnan(powerOrig);
    diffLine=diff(lineMask);

    startLine=find(diffLine==1)+1;
    endLine=find(diffLine==-1);

    if isempty(startLine) & isempty(endLine)
        startLine=1;
        endLine=length(powerOrig);
    end

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

    % Find maxima % and minima
    locsMax=find(locsMaxAll(ii,:)==1);
   
    % Find change points
    diffMax=find(diffMaxAll(ii,:)==1);
    diffMin=find(diffMinAll(ii,:)==1);

    % Find points to test
    testX=[];
    testY=[];

    for kk=1:length(diffMax)-1
        if powerOrig(diffMax(kk+1))>powerOrig(diffMax(kk))
            findMin=sum(diffMin>diffMax(kk) & diffMin<diffMax(kk+1));
            if findMin==1
                hasPeak=sum(locsMax>diffMax(kk) & locsMax<diffMax(kk+1));
                if hasPeak==0
                    testX=cat(1,testX,[diffMax(kk),diffMax(kk+1)]);
                end
            end
        end
    end
    for kk=1:length(diffMin)-1
        if powerOrig(diffMin(kk))>powerOrig(diffMin(kk+1))
            findMax=sum(diffMax>diffMin(kk) & diffMax<diffMin(kk+1));
            if findMax==1
                hasPeak=sum(locsMax>diffMin(kk) & locsMax<diffMin(kk+1));
                if hasPeak==0
                    testX=cat(1,testX,[diffMin(kk),diffMin(kk+1)]);
                end
            end
        end
    end
        
    addPeak=[];
    for kk=1:size(testX,1)
        smallPeakInds=[];
        testXthis=testX(kk,:);
        testYthis=powerOrig(testXthis);
        testY=cat(1,testY,testYthis);
        fitLineMod=polyfit(testXthis,testYthis,1);
        fitLine=polyval(fitLineMod,testX(kk,1):testX(kk,2));
        powerMinusLine=powerOrig(testXthis(1):testXthis(2))-fitLine;
        cutProm=3;
        testProm=max(powerMinusLine)-min(powerMinusLine);
        if testProm>cutProm
            [~,smallPeaks]=max(powerMinusLine);
            smallPeakInds=smallPeaks+testXthis(1);
        end
        if powerOrig(testXthis(1))<powerOrig(testXthis(2))
            for ll=1:length(smallPeakInds)
                goodPeak=smallPeakInds(ll);
                peakValDiff=goodPeak-diffMin;
                getPeak=diffMin(min(find(peakValDiff<0)));
                if getPeak<testXthis(2)
                    addPeak=cat(2,addPeak,getPeak);
                end
            end
        else
            for ll=1:length(smallPeakInds)
                goodPeak=smallPeakInds(ll);
                peakValDiff=goodPeak-diffMax;
                getPeak=diffMax(max(find(peakValDiff>0)));
                if getPeak>testXthis(1)
                    addPeak=cat(2,addPeak,getPeak);
                end
            end
        end
    end

        locsMax=cat(2,locsMax,addPeak);

        %% Add to output

    % Check against noise floor
    if sum(~isnan(locsMax))>2
        minPow=min(powerOrig,[],'omitmissing');
        testPow=powerOrig(locsMax);
        absDiff=abs(testPow-minPow);
        locsMax(absDiff<5)=[];
        if isempty(locsMax)
            [~,locsMax]=max(powerOrig,[],'omitmissing');
        end
    end

    velDual(ii,1:length(locsMax))=specVelIn(ii,locsMax);

    if nn==inf
        scatter(specVelIn(ii,:),powerRaw(ii,:),'MarkerEdgeColor','c');
        xlim([-8,8]);
        hold on
        plot(specVelIn(ii,:),powerOrig,'-b','linewidth',2);
        scatter(specVelIn(ii,locsMax(~isnan(locsMax))),powerOrig(locsMax(~isnan(locsMax))),80,'filled','MarkerFaceColor','r')
        scatter(specVelIn(ii,diffMin),powerOrig(diffMin),'filled','MarkerFaceColor','y')
        scatter(specVelIn(ii,diffMax),powerOrig(diffMax),'filled','MarkerFaceColor','k')
        if ~isempty(testX)
            for mm=1:size(testX,1)
                scatter(specVelIn(ii,testX(mm,:)),testY(mm,:),'filled','MarkerFaceColor','g')
            end
        end
        cla
    end
end
findEmpty=all(isnan(velDual),1);
velDual(:,findEmpty)=[];

checkDims=size(momentsVelDual,3)-size(velDual,2);
if checkDims<0
    momentsVelDual=padarray(momentsVelDual,[0,0,abs(checkDims)],nan,'post');
end

dim3=size(velDual,2);
momentsVelDual(:,nn,1:dim3)=single(velDual);
end