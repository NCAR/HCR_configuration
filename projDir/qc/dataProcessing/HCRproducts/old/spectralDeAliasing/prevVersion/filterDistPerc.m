function [distFilt] = filterDistPerc(distIn,sampleNum)
% Remove lowest 10% of median

distMed=movmedian(distIn,sampleNum/10,2);

distMin=min(distMed,[],2,'omitnan');
distMax=max(distMed,[],2,'omitnan');

spread=distMax-distMin;

distFilt=distMed;
distFilt(distMed<distMax-spread./5)=nan;

end