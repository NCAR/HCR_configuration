function testExt=correctExtremes(testExt,nyq,prevRay)

% Check if extremes exist
if max(testExt,[],'omitnan')<2*nyq-0.5*nyq
    return
end

% Find extreme regions
extThresh=2*nyq-0.5*nyq;
jumpThresh=nyq;

extremes=testExt>nyq+nyq*0.5;

velExt=testExt;
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
    thisStart=testExt(startExt(ll));
    thisEnd=testExt(endExt(ll));

    if (thisStart<extThresh & (thisEnd<extThresh | isnan(endExt(ll)))) | ...
            (thisEnd<extThresh & (thisStart<extThresh | isnan(startExt(ll))))
        continue
    end

    startJump=testExt(startExt(ll))-testExt(startExt(ll)-1);
    endJump=testExt(endExt(ll))-testExt(min([endExt(ll)+1,length(velExt)]));

    medStretch=median(testExt(startExt(ll):endExt(ll)),1,'omitnan');

    if medStretch>nyq+nyq*0.5 & (startJump>jumpThresh | endJump>jumpThresh | isnan(startJump) | isnan(endJump))
        testExt(startExt(ll):endExt(ll))=testExt(startExt(ll):endExt(ll))-2*nyq;
    end

    diffPrev=abs(testExt(startExt(ll):endExt(ll))-prevRay(startExt(ll):endExt(ll)));
    imp=length(find(diffPrev<nyq+0.5*nyq));

    if imp<3
        testExt(startExt(ll):endExt(ll))=testExt(startExt(ll):endExt(ll))+2*nyq;
    end
end