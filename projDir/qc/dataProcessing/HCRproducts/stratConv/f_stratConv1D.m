function stratConv1D=f_stratConv1D(stratConv,asl,icingLevel)
% Convert 2D stratiform convective field to 1D
stratConv1D=nan(1,size(stratConv,2));
stratConv1D(:)=2;

stratConvBelow=stratConv;
stratCount=zeros(size(stratConv1D));
for ii=1:length(stratConv1D)
    aslCol=asl(:,ii);
    stratConvBelow(aslCol>icingLevel(ii),ii)=nan;
    
    stratCount(ii)=sum(stratConv(:,ii)==2);
end

stratConv1D(any(stratConvBelow==1,1))=1;
stratConv1D(stratConv1D==2 & stratCount==0)=1;
stratConv1D(all(isnan(stratConv),1))=nan;
end