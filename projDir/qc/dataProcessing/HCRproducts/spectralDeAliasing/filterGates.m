function maskSpec=filterGates(powerSpec)
% Filter power spectrum on a gate by gate basis
sumSpec=sum(powerSpec,2);
medSpec=movmedian(sumSpec,50);

minSpec=min(medSpec);
sumSpec(sumSpec<minSpec+0.1)=nan;

maskSpec=~isnan(sumSpec);
maskSpec=bwareaopen(maskSpec,20);

plotYes=0;
if plotYes
    subplot(2,1,1)
    plot(sumSpec)
    hold on
    plot(medSpec)

    subplot(2,1,2)
    plot(maskSpec)
    ylim([-1 2])
end

end