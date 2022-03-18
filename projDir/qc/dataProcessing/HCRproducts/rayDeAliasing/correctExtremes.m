function finalRay=correctExtremes(finalRay,nyq,prevRay)

% Check if extremes exist
if max(finalRay,[],'omitnan')<2*nyq-0.5*nyq
    return
end

% Find extreme regions
extThresh=2*nyq-0.5*nyq;
jumpThresh=nyq;

extremes=finalRay>nyq+nyq*0.5;

velExt=finalRay;
velExt(~extremes)=nan;
velExt=movmean(velExt,5,'omitnan');
velExt=movmean(velExt,5,'includenan');
extremes(~isnan(velExt))=1;

diffExt=diff(extremes);
startExt=find(diffExt==1)+1;
endExt=find(diffExt==-1);

if length(startExt)>length(endExt)
    endExt=cat(1,endExt,length(extremes));
end

for ll=1:length(startExt)
    % Check if end points are extreme
    thisStart=finalRay(startExt(ll));
    thisEnd=finalRay(endExt(ll));

    if (thisStart<extThresh & (thisEnd<extThresh | isnan(endExt(ll)))) | ...
            (thisEnd<extThresh & (thisStart<extThresh | isnan(startExt(ll))))
        continue
    end

    startJump=finalRay(startExt(ll))-finalRay(startExt(ll)-1);
    endJump=finalRay(endExt(ll))-finalRay(endExt(ll)+1);

    medStretch=median(finalRay(startExt(ll):endExt(ll)),1,'omitnan');

    if medStretch>nyq+nyq*0.5 & (startJump>jumpThresh | endJump>jumpThresh | isnan(startJump) | isnan(endJump))
        finalRay(startExt(ll):endExt(ll))=finalRay(startExt(ll):endExt(ll))-2*nyq;
    end

    diffPrev=abs(finalRay(startExt(ll):endExt(ll))-prevRay(startExt(ll):endExt(ll)));
    minDiff=min(diffPrev,[],'omitnan');

    if minDiff>nyq+0.5*nyq
        finalRay(startExt(ll):endExt(ll))=finalRay(startExt(ll):endExt(ll))+2*nyq;
    end
end