function [powerSpecFilt] = filterPowerSpecPerc(powerSpecLarge,sampleNum)
% Remove lowest 10% of median

powerSpecMed=movmedian(powerSpecLarge,sampleNum/10,2);

powerSpecMin=min(powerSpecMed,[],2,'omitnan');

spread=max(powerSpecMed,[],2,'omitnan')-powerSpecMin;

powerSpecFilt=powerSpecMed;
powerSpecFilt(powerSpecMed<powerSpecMin+spread./2)=nan;
powerSpecFilt(movmedian(spread,50)<5,:)=nan;

end