function [refSig0,surfFlag,refFlag]=makeRefSig0(sig0meas,sig0model,surfFlag)
% Combine measured and model sig0. Use measured whenever possible and model
% when no measured close by.

refSig0=nan(size(sig0meas));
refFlag=ones(size(sig0meas));

sig0meas(surfFlag~=2)=nan;
% In CSET we have cases where the power was lowered close to the surface
sig0meas(sig0meas<0)=nan;

% Clean up data
% First remove outliers
sig0medLarge=movmedian(sig0meas,2000,'omitnan');

% Remove data where distance to movmedLarge is too big
diffMeasMed=abs(sig0meas-sig0medLarge);
sig0meas(diffMeasMed>1)=nan;

% Remove short data stretches
sig0mask=zeros(size(sig0meas));
sig0mask(~isnan(movmedian(sig0meas,3,'omitnan')))=1;
sig0mask=bwareaopen(sig0mask,10);

sig0meas(~isnan(sig0meas) & ~sig0mask)=nan;

% Find areas with no measurement close by
sig0med=movmedian(sig0meas,3000,'omitnan');

refSig0=movmedian(sig0meas,100,'omitnan');
refSig0(isnan(sig0meas))=nan;
refSig0(isnan(sig0med))=sig0model(isnan(sig0med));

refFlag(isnan(sig0med))=3;
refFlag(isnan(refSig0))=2;

surfFlag(surfFlag==2 & refFlag~=1)=0;

refSig0=fillmissing(refSig0,'linear','EndValues','nearest');
end