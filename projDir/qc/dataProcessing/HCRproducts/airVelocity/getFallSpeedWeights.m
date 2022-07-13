function [wr,wi]=getFallSpeedWeights(echo,icingIn,melt,temp,asl)
% Get fall speed weights
% https://doi.org/10.1175/2009JAS3132.1

wr=nan(size(echo));
wi=nan(size(echo));

% No icing
icing=icingIn;
modeMelt=mode(melt,1);
icing(isnan(icing) & modeMelt<20)=inf;
icing(isnan(icing) & modeMelt>=20)=-inf;

% Below melting layer
iceMat=repmat(icing,size(asl,1),1);
wr(asl<iceMat)=1;
wi(asl<iceMat)=0;

% Stratiform
% T < 0 C
wr(temp<=0 & ~isnan(temp) & echo<35 & echo~=30)=0;
wi(temp<=0 & ~isnan(temp) & echo<35 & echo~=30)=1;

% Find ht
htMaskS=asl;
htMaskS(temp>0 | isnan(temp))=nan;
hts=min(htMaskS,[],1,'omitnan');
hts=fillmissing(hts,'linear','EndValues','nearest');

% Transition stratiform
tsInds=find(~isnan(temp) & isnan(wr) & echo<35 & echo~=30);
wrInit=(hts-asl).^2./(hts-icing).^2;
wiInit=1-wrInit;
wr(tsInds)=wrInit(tsInds);
wi(tsInds)=wiInit(tsInds);

% Convective
% T < -12 C
wr(temp<-12 & ~isnan(temp) & (echo>35 | echo==30))=0;
wi(temp<-12 & ~isnan(temp) & (echo>35 | echo==30))=1;

% Find ht
htMaskC=asl;
htMaskC(temp>=-12 | isnan(temp))=nan;
htc=min(htMaskC,[],1,'omitnan');
htc=fillmissing(htc,'linear','EndValues','nearest');

% Transition convective
tsInds=find(~isnan(temp) & isnan(wr) & (echo>35 | echo==30));
wrInit=(htc-asl).^2./(htc-icing).^2;
wiInit=1-wrInit;
wr(tsInds)=wrInit(tsInds);
wi(tsInds)=wiInit(tsInds);
end