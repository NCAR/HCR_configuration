function [momentsVelDual,momentsPowDual,leftSvel,rightSvel,leftPvel,rightPvel]= ...
    findMultiVels(powerIn,specVelIn,powerRaw,momentsVelDual,momentsPowDual,leftSvel,rightSvel,leftPvel,rightPvel,nn)
% Find maxima and minima in spectra
velDual=nan(size(powerIn,1),35);
powDual=nan(size(powerIn,1),35);

% Remove lower part
maxIn=max(powerIn,[],2,'omitmissing');
minIn=min(powerIn,[],2,'omitmissing');

powSpread=(maxIn-minIn);
powerIn(powerIn<minIn+powSpread/3)=nan;

dataInds=find(any(~isnan(powerIn),2));

% Remove transmitter pulse
dataInds(dataInds<=12)=[];

[locsMaxAll,maxPromAll]=islocalmax(powerIn,2);

diff1Curve=cat(2,nan(size(powerIn,1),1),diff(powerIn,1,2));
diff1MaxAll=islocalmax(diff1Curve,2);
diff1MinAll=islocalmin(diff1Curve,2);

% For shoulder points
diff2Curve=cat(2,diff(diff1Curve,1,2),nan(size(powerIn,1),1));
maxIn2=max(powerIn,[],2,'omitmissing');
minIn2=min(powerIn,[],2,'omitmissing');

powSpread2=(maxIn2-minIn2);
diff2Curve(powSpread2<10,:)=nan;

diff2MinAll=islocalmin(diff2Curve,2);

for jj=1:length(dataInds)
    ii=dataInds(jj);
    
    powerOrig=powerIn(ii,:);

    % Find maxima % and minima
    locsMax=find(locsMaxAll(ii,:)==1);
    maxProm=maxPromAll(ii,locsMax);
       
    % Find change points
    diff1Max=find(diff1MaxAll(ii,:)==1);
    diff1Min=find(diff1MinAll(ii,:)==1);

    % Find steep points
    diff2Min=find(diff2MinAll(ii,:)==1);
    leftShoulder=[];
    rightShoulder=[];
    if ~isempty(diff2Min)
        leftShoulder=diff2Min(1);
        rightShoulder=diff2Min(end);
        if leftShoulder==rightShoulder
            leftShoulder=[];
            rightShoulder=[];
        end
    end
    if ~isempty(leftShoulder)
        maxProm(locsMax<=leftShoulder)=[];
        maxProm(locsMax>=rightShoulder)=[];
        locsMax(locsMax<=leftShoulder)=[];
        locsMax(locsMax>=rightShoulder)=[];
        
        % Output
        leftSvel(ii,nn)=specVelIn(ii,leftShoulder);
        rightSvel(ii,nn)=specVelIn(ii,rightShoulder);
        leftPvel(ii,nn)=powerIn(ii,leftShoulder);
        rightPvel(ii,nn)=powerIn(ii,rightShoulder);
    end

    % Find points to test
    testX=[];
    testY=[];

    for kk=1:length(diff1Max)-1
        if powerOrig(diff1Max(kk+1))>powerOrig(diff1Max(kk))
            findMin=sum(diff1Min>diff1Max(kk) & diff1Min<diff1Max(kk+1));
            if findMin==1
                hasPeak=sum(locsMax>diff1Max(kk) & locsMax<diff1Max(kk+1));
                if hasPeak==0
                    testX=cat(1,testX,[diff1Max(kk),diff1Max(kk+1)]);
                end
            end
        end
    end
    for kk=1:length(diff1Min)-1
        if powerOrig(diff1Min(kk))>powerOrig(diff1Min(kk+1))
            findMax=sum(diff1Max>diff1Min(kk) & diff1Max<diff1Min(kk+1));
            if findMax==1
                hasPeak=sum(locsMax>diff1Min(kk) & locsMax<diff1Min(kk+1));
                if hasPeak==0
                    testX=cat(1,testX,[diff1Min(kk),diff1Min(kk+1)]);
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
                peakValDiff=goodPeak-diff1Min;
                getPeak=diff1Min(min(find(peakValDiff<0)));
                if getPeak<testXthis(2)
                    addPeak=cat(2,addPeak,getPeak);
                end
            end
        else
            for ll=1:length(smallPeakInds)
                goodPeak=smallPeakInds(ll);
                peakValDiff=goodPeak-diff1Max;
                getPeak=diff1Max(max(find(peakValDiff>0)));
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
    powDual(ii,1:length(locsMax))=powerIn(ii,locsMax);

    if nn==1
        scatter(specVelIn(ii,:),powerRaw(ii,:),'MarkerEdgeColor','k');
        xlim([-8,8]);
        hold on
        plot(specVelIn(ii,:),powerOrig,'-b','linewidth',2);
        scatter(specVelIn(ii,locsMax(~isnan(locsMax))),powerOrig(locsMax(~isnan(locsMax))),80,'filled','MarkerFaceColor','r')
        scatter(specVelIn(ii,leftShoulder),powerOrig(leftShoulder),80,'filled','MarkerFaceColor','m')
        scatter(specVelIn(ii,rightShoulder),powerOrig(rightShoulder),80,'filled','MarkerFaceColor','m')
        scatter(specVelIn(ii,diff1Min),powerOrig(diff1Min),'filled','MarkerFaceColor','y')
        scatter(specVelIn(ii,diff1Max),powerOrig(diff1Max),'filled','MarkerFaceColor','k')
        if ~isempty(testX)
            for mm=1:size(testX,1)
                scatter(specVelIn(ii,testX(mm,:)),testY(mm,:),'filled','MarkerFaceColor','g')
            end
        end

        % yyaxis right
        % plot(specVelIn(ii,:),diff2Curve(ii,:),'-c','linewidth',1.5);
        % hold on
        % scatter(diffMin,diffCurve(diffMin),'filled','MarkerFaceColor','y')
        % scatter(diffMax,diffCurve(diffMax),'filled','MarkerFaceColor','k')
        %ylim([-0.5,0.5]);
        grid on

        cla
        % yyaxis left
        % cla
    end
end
findEmpty=all(isnan(velDual),1);
velDual(:,findEmpty)=[];
powDual(:,findEmpty)=[];

checkDims=size(momentsVelDual,3)-size(velDual,2);
if checkDims<0
    momentsVelDual=padarray(momentsVelDual,[0,0,abs(checkDims)],nan,'post');
    momentsPowDual=padarray(momentsPowDual,[0,0,abs(checkDims)],nan,'post');
end

dim3=size(velDual,2);
momentsVelDual(:,nn,1:dim3)=single(velDual);
momentsPowDual(:,nn,1:dim3)=single(powDual);
end