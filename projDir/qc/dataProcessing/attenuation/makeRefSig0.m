function [refSig0,surfFlag]=makeRefSig0(sig0meas,sig0model,surfFlag)
% Combine measured and model sig0. Use measured whenever possible and model
% when no measured close by.

sig0measClear=sig0meas;
sig0measClear(surfFlag~=2)=nan;
% In CSET we have cases where the power was lowered close to the surface
sig0measClear(sig0measClear<0)=nan;

% Clean up measured sig0 data
% First remove outliers
sig0medLarge=movmedian(sig0measClear,2001,'omitnan');

% Remove data where distance to movmedLarge is too big
diffMeasMed=abs(sig0measClear-sig0medLarge);
sig0measClear(diffMeasMed>1)=nan;

% Remove short data stretches
sig0mask=zeros(size(sig0measClear));
sig0mask(~isnan(movmedian(sig0measClear,3,'omitnan')))=1;
sig0mask=bwareaopen(sig0mask,10);

sig0measClear(~isnan(sig0measClear) & ~sig0mask)=nan;

% Smooth slightly
refSig01=movmedian(sig0measClear,101,'omitnan');
refSig01(isnan(sig0measClear))=nan;

% Find areas with no measurement close by
sig0med=movmedian(sig0measClear,6001,'omitnan');
refSig01(isnan(sig0med))=sig0model(isnan(sig0med));

% refFlag
% 1=measured
% 2=interpolated
% 3=model
refFlag=ones(size(sig0meas));
refFlag(isnan(sig0med))=3;
refFlag(isnan(refSig01))=2;

refSig01=fillmissing(refSig01,'linear','EndValues','nearest');

% Use maximum of measured and model
measMod=cat(1,refSig01,sig0model);
[refSig0,maxInd]=max(measMod,[],1,'omitnan');

%refFlag(maxInd==2)=3;

% Remove "clear air" surface flag where not measured
surfFlag(surfFlag==2 & refFlag~=1)=0;
end