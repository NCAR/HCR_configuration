function [finalRay] = deAliasSingleRay(velRay,velPrev,nyq,plotStart,time)
% De-alias single ray
%% Check if folding occurs

diffPrevInit=velRay-velPrev;
diffRay=diff(velRay);

if all(isnan(velRay)) | (max(abs(diffRay),[],'omitnan')<2*nyq-0.5*nyq & median(abs(diffPrevInit))<2*nyq-0.5*nyq)
    finalRay=velRay;
else

    %% Make ray consistent within itself

    % Find folds
    diffVelGates=diff(velRay);
    foldInds=find(abs(diffVelGates)>2*nyq-0.5*nyq);
    bigJumpMag=diffVelGates(foldInds);

    vertConsRay=velRay;

    for jj=1:length(foldInds)
        if bigJumpMag(jj)>0
            vertConsRay(foldInds(jj)+1:end)=vertConsRay(foldInds(jj)+1:end)-2*nyq;
        else
            vertConsRay(foldInds(jj)+1:end)=vertConsRay(foldInds(jj)+1:end)+2*nyq;
        end
    end

    %% Make consistent with previous

    diffPrev=vertConsRay-velPrev;
    timeConsRay=vertConsRay;

    maxFolding=ceil(max(abs(diffPrev))/(2*nyq));

    for kk=1:maxFolding
        indsPos=find(diffPrev>(kk*2*nyq-0.5*nyq));
        timeConsRay(indsPos)=vertConsRay(indsPos)-kk*2*nyq;
        indsNeg=find(diffPrev<-(kk*2*nyq-0.5*nyq));
        timeConsRay(indsNeg)=vertConsRay(indsNeg)+kk*2*nyq;
    end

    %% Moving median test

    movMedRay=timeConsRay;
    medTimeCons=movmedian(movMedRay,75,'omitnan'); % Was 50
    diffMed=movMedRay-medTimeCons;
    movMedRay(diffMed>nyq*0.8)=movMedRay(diffMed>nyq*0.8)-2*nyq;
    movMedRay(diffMed<-(nyq*0.8))=movMedRay(diffMed<-(nyq*0.8))+2*nyq;

    medTimeCons2=movmedian(movMedRay,5,'omitnan');
    diffMed2=movMedRay-medTimeCons2;
    movMedRay(diffMed2>nyq*0.8)=movMedRay(diffMed2>nyq*0.8)-2*nyq;
    movMedRay(diffMed2<-(nyq*0.8))=movMedRay(diffMed2<-(nyq*0.8))+2*nyq;

    %% Check for outliers

    outliersRay=movMedRay;
    diffPrevTime=outliersRay-velPrev;

    nanInds=isnan(outliersRay);
    nanDiff=diff(nanInds);
    nanEndInds=find(nanDiff~=0);

    % Check if jumps occur
    diffTimeConsGates=diff(outliersRay);

    countWhile=0;

    while max(abs(diffTimeConsGates))>nyq

        bigJumpInds=find(abs(diffTimeConsGates)>nyq);

        smallJumpInds=find(abs(diffTimeConsGates)>nyq*0.5);
        allInds=cat(1,smallJumpInds,nanEndInds);
        allInds=sort(allInds);

        % Check if the problem is before or after the jump
        beforePrev=outliersRay(bigJumpInds(1))-velPrev(bigJumpInds(1));
        afterPrev=outliersRay(bigJumpInds(1)+1)-velPrev(bigJumpInds(1)+1);

        if abs(beforePrev)>=abs(afterPrev)
            endStretch=bigJumpInds(1);
            indInAll=find(allInds==endStretch);
            startStretch=allInds(indInAll-1)+1;

            meanDiffPrev=mean(diffPrevTime(startStretch:endStretch),'omitnan');

            if beforePrev>0 & abs(meanDiffPrev)>nyq*0.5
                outliersRay(startStretch:endStretch)=outliersRay(startStretch:endStretch)-2*nyq;
            elseif beforePrev<0 & abs(meanDiffPrev)>nyq*0.5
                outliersRay(startStretch:endStretch)=outliersRay(startStretch:endStretch)+2*nyq;
            end
        else
            startStretch=bigJumpInds(1)+1;
            indInAll=find(allInds==startStretch-1);
            if indInAll<length(allInds)
                endStretch=allInds(indInAll+1);
            else
                endStretch=length(outliersRay);
            end

            meanDiffPrev=mean(diffPrevTime(startStretch:endStretch),'omitnan');

            if afterPrev>0 & abs(meanDiffPrev)>nyq*0.5
                outliersRay(startStretch:endStretch)=outliersRay(startStretch:endStretch)-2*nyq;
            elseif afterPrev<0 & abs(meanDiffPrev)>nyq*0.5
                outliersRay(startStretch:endStretch)=outliersRay(startStretch:endStretch)+2*nyq;
            end
        end

        testRay=outliersRay;
        testRay(1:endStretch)=nan;

        countWhile=countWhile+1;

        diffTimeConsGates=diff(testRay);

        if countWhile>100
            warning('Jumping out of while loop ...');
            diffTimeConsGates=0;
        end
    end

    %% Correct extremes
    testExtPos=outliersRay;

    % Positive
    testExtPos=correctExtremes(testExtPos,nyq,velPrev);

    % Negative
    testExtNeg=-testExtPos;
    testExtNeg=correctExtremes(testExtNeg,nyq,-velPrev);
    finalRay=-testExtNeg;

    %% Test plot

    if ~isempty(plotStart) & time>plotStart
        plot(velPrev)
        hold on
        plot(velRay);
        plot(vertConsRay)
        plot(timeConsRay)
        plot(movMedRay)
        plot(outliersRay)
        plot(testExtPos)
        plot(finalRay,'-k','LineWidth',2)

        hold off
        xlim([1 500])
        ylim([-25 25])
        title(datestr(time));
        legend('velPrev','orig','vertCons','timeCons','movMed','outliers','finalPos','final')
        stopHere=1;
    end
end
end

