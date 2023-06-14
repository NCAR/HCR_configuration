function LDR=preProcessLDR(LDR)

LDR(1:20,:)=nan;

% Remove small areas
maskLDRtemp=~isnan(LDR);
maskLDRtemp=bwareaopen(maskLDRtemp,5000);

LDR(~isnan(LDR) & maskLDRtemp==0)=nan;

maskLDRtemp=~isnan(LDR);

stats=regionprops('table',maskLDRtemp,'Solidity','PixelIdxList');

for ll=1:size(stats,1)
    if stats.Solidity(ll)<0.5
        maskLDRtemp(stats.PixelIdxList{ll})=0;
    end
end
LDR(~isnan(LDR) & maskLDRtemp==0)=nan;
end