function minOut=hitschfeldBordan_optimize(zIn,pia,rangeIn,alpha,beta,firstInd)
% Calculate attenuation corrected reflectivity with the Hitschfeld Bordan
% solution.
% Fairall et al. (2018) https://doi.org/10.1175/JTECH-D-17-0025.1 (Eqs. 20-22)

zInLin=10.^(zIn./10);

As=10.^(-pia./10);

q=0.2*beta*log(10);

sMat=alpha.*zInLin.^beta.*(rangeIn(2)-rangeIn(1))./1000;

Sdiff=double(cumsum(sMat,1,'reverse','omitnan'));
Sdiff(isnan(zIn))=nan;

zHBlin=zInLin./((As.^beta+q.*Sdiff).^(1/beta));

zOut=10.*log10(zHBlin);

% Compute minimum of first non-nans
zDiff=zOut-zIn;
goodInds=find(~isnan(firstInd));
topZdiff=[];
for ii=1:length(goodInds)
    testZ=zDiff(firstInd(goodInds(ii)),ii);
    if ~isnan(testZ)
        topZdiff=[topZdiff,zDiff(firstInd(goodInds(ii)),ii)];
    end
end
topZdiffAbs=abs(topZdiff);
minOut=sum(topZdiffAbs);
end