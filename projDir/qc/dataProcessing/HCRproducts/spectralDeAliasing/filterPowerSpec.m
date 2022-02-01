function [powerSpecFilt,powerSpecMed,powerSpecMed2] = filterPowerSpec(powerSpecLarge,sampleNum)
% Remove noise from power spectrum

powerSpecMed=movmedian(powerSpecLarge,sampleNum/10,2);

powerSpecMin=min(powerSpecMed,[],2,'omitnan');

powerSpecLow=powerSpecLarge;
powerSpecLow(powerSpecMed>powerSpecMin+5)=nan;

powerSpecMed2=movmedian(powerSpecLow,sampleNum/10,2,'includenan');
powerSpecMed2=fillmissing(powerSpecMed2,'nearest',2);

medDiff=powerSpecMed-powerSpecMed2;

powerSpecFilt=powerSpecLarge;
powerSpecFilt(medDiff<1)=nan;
powerSpecFilt(powerSpecLarge<powerSpecMin)=nan;

end