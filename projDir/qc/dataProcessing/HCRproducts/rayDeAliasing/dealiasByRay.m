function velDeAliased=dealiasByRay(velIn,elev,nyq,dataFreq,time,plotStart)
velDeAliased=nan(size(velIn));

defaultPrev=nyq;

velPrev=repmat(defaultPrev,size(velIn,1),1);
prevCount=zeros(size(velPrev));
prevKeep=nan(size(velPrev));
flipYes=0;

for ii=1:size(velIn,2)
    if elev(ii)>0
        velRay=-velIn(:,ii);
    else
        velRay=velIn(:,ii);
    end

    finalRay=deAliasSingleRay(velRay,velPrev,nyq,plotStart,time(ii));

    %% Set up time consistency check

    [velPrev,prevCount,prevKeep,flipYes]=setUpPrev(finalRay,velPrev,prevCount,prevKeep,flipYes,elev(ii),dataFreq,defaultPrev);

    %% Add to output
    if elev(ii)>0
        velDeAliased(:,ii)=-finalRay;
    else
        velDeAliased(:,ii)=finalRay;
    end
end
end