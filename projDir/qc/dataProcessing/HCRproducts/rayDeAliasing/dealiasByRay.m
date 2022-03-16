function velDeAliased=dealiasByRay(velIn,elev,nyq,dataFreq,time,plotStart)
velDeAliased=nan(size(velIn));

%velIn(:,elev>0)=-velIn(:,elev>0);

nonNanInds=find(any(~isnan(velIn),1));

velPrev=repmat(2,size(velIn,1),1);
prevCount=zeros(size(velPrev));
prevKeep=nan(size(velPrev));

for ii=1:length(nonNanInds)
    velRay=velIn(:,nonNanInds(ii));

    %% Check if folding occurs

    diffPrevInit=velRay-velPrev;
    maxPrev=max(velPrev(~isnan(velRay)),[],'omitnan');
    diffRay=diff(velRay);

    if max(abs(diffRay),[],'omitnan')<2*nyq-0.5*nyq & maxPrev<0
        finalRay=velRay;
    elseif max(abs(diffRay),[],'omitnan')<2*nyq-0.5*nyq & median(abs(diffPrevInit))<2*nyq-0.5*nyq
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

        maxFolding=floor(max(abs(diffPrev))/(2*nyq));

        for kk=1:maxFolding
            indsPos=find(diffPrev>(kk*2*nyq-0.5*nyq));
            timeConsRay(indsPos)=vertConsRay(indsPos)-kk*2*nyq;
            indsNeg=find(diffPrev<-(kk*2*nyq-0.5*nyq));
            timeConsRay(indsNeg)=vertConsRay(indsNeg)+kk*2*nyq;
        end

        %% Moving median test

        movMedRay=timeConsRay;
        medTimeCons=movmedian(movMedRay,50,'omitnan');
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
                endStretch=allInds(indInAll+1);

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
        finalRay=outliersRay;

        % Positive
        finalRay=correctExtremes(finalRay,nyq);

        % Negative
        finalRay=-finalRay;
        finalRay=correctExtremes(finalRay,nyq);
        finalRay=-finalRay;

        %% Test plot

        if time(nonNanInds(ii))>plotStart
            plot(velPrev)
            hold on
            plot(velRay);
            plot(vertConsRay)
            plot(timeConsRay)
            plot(movMedRay)
            plot(outliersRay)
            plot(finalRay,'LineWidth',2)

            hold off
            xlim([1 500])
            ylim([-25 25])
            title(datestr(time(nonNanInds(ii))));
            legend('velPrev','orig','vertCons','timeCons','movMed','outliers','final')
            stopHere=1;
        end
    end
    %% Add to output

    velDeAliased(:,nonNanInds(ii))=finalRay;

    %% Set up time consistency check

    % Decide which velocities to keep for time consistency check: join
    % close regions and remove small isolated regions
    velForPrev=finalRay;
    velForPrev=movmean(velForPrev,5,'omitnan');
    velForPrev=movmean(velForPrev,5,'includenan');
    velForPrevMask=~isnan(velForPrev);
    velForPrevMask=bwareaopen(velForPrevMask,15);

    % Create new velPrev
    velPrev=movmedian(finalRay,20,'omitnan');
    velPrev(isnan(finalRay))=nan;
    velPrev(~velForPrevMask)=nan;

    % Add new values to prevKeep and handle counts
    prevKeep(velForPrevMask)=velPrev(velForPrevMask);
    prevCount(velForPrevMask)=0;
    prevCount(~velForPrevMask)=prevCount(~velForPrevMask)+1;
    prevKeep(prevCount>=dataFreq*5)=nan;

    % Add old velocities
    velPrev(isnan(velPrev))=prevKeep(isnan(velPrev));

    % Interpolate prev
    prevMedLarge=movmedian(velPrev,100,'omitnan');

    velPrev(isnan(prevMedLarge))=2;
    velPrev=fillmissing(velPrev,'linear','EndValues','nearest');

end
%velDeAliased(:,elev>0)=-velDeAliased(:,elev>0);
end