function [LDRgrad,smoothLDR,iceMask]=findIce(data)

%% LDR gradient

disp('Processing LDR ...')
smoothLDR=movmean(data.LDR,40,1,'includenan');

LDRgrad=diff(smoothLDR,1);

%% Ice regions

iceRegs=LDRgrad;
iceRegs(data.asl(1:end-1,:)<data.ICING_LEVEL+1000)=nan;

iceMask=iceRegs>0.1;
iceMask=imclose(iceMask,strel('disk',7));
iceMask=imfill(iceMask,'holes');

iceMask=bwareaopen(iceMask,1000);
end