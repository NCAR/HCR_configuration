function velDeAliased=dealiasByRay(velIn,elev,nyq,dataFreq,time,plotStart)
velDeAliased=nan(size(velIn));

%velIn(:,elev>0)=-velIn(:,elev>0);

nonNanInds=find(any(~isnan(velIn),1));

velPrev=repmat(2,size(velIn,1),1);
prevCount=zeros(size(velPrev));
prevKeep=nan(size(velPrev));

for ii=1:length(nonNanInds)
    velRay=velIn(:,nonNanInds(ii));

    %% Make ray consistent within itself

    % Find folds
    diffVelGates=diff(velRay);
    foldInds=find(abs(diffVelGates)>2*nyq-0.5*nyq);
    diffFoldInds=diffVelGates(foldInds);

    vertConsRay=velRay;

    for jj=1:length(foldInds)
        if diffFoldInds(jj)>0
            vertConsRay(foldInds(jj)+1:end)=vertConsRay(foldInds(jj)+1:end)-2*nyq;
        else
            vertConsRay(foldInds(jj)+1:end)=vertConsRay(foldInds(jj)+1:end)+2*nyq;
        end
    end

    %% Make consistent with previous

    diffPrev=vertConsRay-velPrev;
    horConsRay=vertConsRay;

    maxFolding=floor(max(abs(diffPrev))/(2*nyq));

    for kk=1:maxFolding
        indsPos=find(diffPrev>(kk*2*nyq-0.5*nyq));
        horConsRay(indsPos)=vertConsRay(indsPos)-kk*2*nyq;
        indsNeg=find(diffPrev<-(kk*2*nyq-0.5*nyq));
        horConsRay(indsNeg)=vertConsRay(indsNeg)+kk*2*nyq;
    end

    %if time(nonNanInds(ii))>plotStart
    plot(velRay);
    hold on
    plot(vertConsRay)
    plot(velPrev)
    plot(horConsRay)
    hold off
    xlim([1 500])
    title(datestr(time(nonNanInds(ii))));
    stopHere=1;
    %end

    % Handle previous velocity
    velForPrev=horConsRay;
    velForPrev=movmean(velForPrev,5,'omitnan');
    velForPrev=movmean(velForPrev,5,'includenan');
    velForPrevMask=~isnan(velForPrev);
    velForPrevMask=bwareaopen(velForPrevMask,10);

    % Replenish old max indeces
    prevKeep(velForPrevMask)=horConsRay(velForPrevMask);
    prevCount(velForPrevMask)=0;
    prevCount(~velForPrevMask)=prevCount(~velForPrevMask)+1;
    prevKeep(prevCount>=dataFreq*5)=nan;

    % Add old max indeces
    velPrev=movmedian(horConsRay,20,'omitnan');
    velPrev(isnan(horConsRay))=nan;
    velPrev(~velForPrevMask)=nan;

    velPrev(isnan(velPrev))=prevKeep(isnan(velPrev));

    % Interpolate prev
    %prevMed=movmedian(velPrev,50,'omitnan');
    prevMedLarge=movmedian(velPrev,100,'omitnan');

    velPrev(isnan(prevMedLarge))=2;
    velPrev=fillmissing(velPrev,'linear','EndValues','nearest');

    velDeAliased(:,nonNanInds(ii))=horConsRay;

end
%velDeAliased(:,elev>0)=-velDeAliased(:,elev>0);
end