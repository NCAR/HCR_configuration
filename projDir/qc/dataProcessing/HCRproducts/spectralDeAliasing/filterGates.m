function maskSpec=filterGates(powerSpec)
% Filter power spectrum on a gate by gate basis
sumSpec=sum(powerSpec,2);
medSpec=movmedian(sumSpec,50);

minSpec=min(medSpec);
sumSpec(sumSpec<minSpec+0.1)=nan;

maskSpec=~isnan(sumSpec);
maskSpec=bwareaopen(maskSpec,20);

end