function stratConv1D=f_stratConv1Dperc(stratConv)
% Convert 2D stratiform convective field to 1D
stratConv1D=nan(1,size(stratConv,2));

sumConv=sum(stratConv==10,1)+sum(stratConv==11,1)+sum(stratConv==12,1)+sum(stratConv==13,1);
sumStrat=sum(stratConv==20,1)+sum(stratConv==21,1);

sumTot=sumStrat+sumConv;

stratConv1D=sumConv./sumTot;

% Remove rays with no data
stratConv1D(all(isnan(stratConv),1))=nan;
end