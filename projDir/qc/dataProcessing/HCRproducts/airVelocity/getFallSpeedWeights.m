function [wr,wi]=getFallSpeedWeights(echo,icingIn,melt,temp,asl)
% Get fall speed weights
% https://doi.org/10.1175/2009JAS3132.1

wr=nan(size(echo));
wi=nan(size(echo));

% No icing
icing=icingIn;
modeMelt=mode(melt,1);
icing(isnan(icing) & modeMelt<20)=0;
icing(isnan(icing) & modeMelt>=20)=inf;

% Below melting layer
iceMat=repmat(icing,size(asl,1),1);
wr(asl<iceMat)=1;
wi(asl<iceMat)=0;

% T < -12 C
wr(temp<-12)=0;
wi(temp<-12)=1;

% Find ht
htMask=asl;
htMask(temp>=-12 | isnan(temp))=nan;
ht=min(htMask,[],1,'omitnan');

% Transition stratiform
%tsInds=find(~isnan(temp) & isnan(wr) & echo<20);
tsInds=find(~isnan(temp) & isnan(wr));
wrInit=(ht-asl).^2./(ht-icing).^2;
wiInit=1-wrInit;
wr(tsInds)=wrInit(tsInds);
wi(tsInds)=wiInit(tsInds);
end