function zOut=hitschfeldBordan_surfUp(zIn,pia,rangeIn)
%Calculate attenuation corrected reflectivity with the Hitschfeld Bordan
%solution

zInLin=10.^(zIn./10);

alpha=0.05;
beta=1;

As=10.^(-pia./10);

q=0.2*beta*log(10);

sMat=alpha.*zInLin.^beta.*(rangeIn(2)-rangeIn(1));

Sdiff=double(cumsum(sMat,1,'reverse','omitnan'));
Sdiff(isnan(zIn))=nan;

%zHB=zIn.*(1-q.*Sr).^-(1/beta);
%zHB=zIn.*exp(-0.2.*log(10).*Sr);

zHBlin=zInLin./((As.^beta+q.*Sdiff).^(1/beta));

zOut=10.*log10(zHBlin);
%zOut=real(zOut);
end