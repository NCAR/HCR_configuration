function finalRay=correctExtremes(finalRay,nyq)

% Check if extremes exist
if max(finalRay,[],'omitnan')<2*nyq-0.5*nyq
    return
end

% Find extreme regions
extThresh=2*nyq-0.5*nyq;
jumpThresh=nyq;

extremes=finalRay>nyq;

velExt=finalRay;
velExt(~extremes)=nan;
velExt=movmean(velExt,9,'omitnan');
velExt=movmean(velExt,9,'includenan');
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
    endJump=finalRay(endExt(ll)+1)-finalRay(endExt(ll));

    if startJump>jumpThresh | endJump>jumpThresh
        finalRay(startExt(ll):endExt(ll))=finalRay(startExt(ll):endExt(ll))-2*nyq;
    end
end