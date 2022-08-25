function zOut=hitschfeldBordan(zIn,rangeIn)
%Calculate attenuation corrected reflectivity with the Hitschfeld Bordan
%solution

%zInLin=10.^(zIn./10);

alpha=0.002;
beta=0.808;

q=0.2*beta*log(10);

sMat=alpha.*zIn.^beta.*(rangeIn(2)-rangeIn(1));

Sr=double(cumsum(sMat,1,'omitnan'));

%zHB=zIn.*(1-q.*Sr).^-(1/beta);
zHB=zIn.*exp(-0.2.*log(10).*Sr);

zOut=zHB;
zOut=real(zOut);
end