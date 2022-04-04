function velOut=filterGatesTD(velIn)
% Filter velocity on a gate by gate basis

stdVel=movstd(velIn,3);

velOut=velIn;

velMask=velIn;
velMask(stdVel>1)=nan;

% Close gaps
velMask=movmean(velMask,5,'omitnan');
velMask=movmean(velMask,5,'includenan');

% Remove small objects
velMask=~isnan(velMask);
velMask=bwareaopen(velMask,20);

velOut=velIn;
velOut(1:16)=nan;
velOut(velMask==0)=nan;

plotYes=0;
if plotYes
    subplot(2,1,1)
    plot(velIn)
    hold on
    plot(velOut)

    subplot(2,1,2)
    plot(stdVel)
    ylim([0 10])
    stopHere=1;
end

end