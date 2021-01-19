function stratConv1D=f_stratConv1D(stratConv,meltingLayer)
% Convert 2D stratiform convective field to 1D
stratConv1D=nan(1,size(stratConv,2));
stratConv1D(:)=2;

stratConvBelow=stratConv;
stratConvBelow(meltingLayer~=10)=nan;
stratCount=sum(stratConv>19,1);

stratConv1D(any(stratConvBelow<12,1))=1;
stratConv1D(stratConv1D==2 & stratCount==0)=1;
stratConv1D(all(isnan(stratConv),1))=nan;
end