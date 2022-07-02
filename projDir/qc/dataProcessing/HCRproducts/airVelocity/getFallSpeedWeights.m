function [wr,wi]=getFallSpeedWeights(echo,melt,icing,temp,asl)
% Get fall speed weights
% https://doi.org/10.1175/2009JAS3132.1

wr=nan(size(echo));
wi=nan(size(echo));

% Below melting layer
wr(melt<20)=1;
wi(melt<20)=0;

% T < -12 C
wr(temp<-12)=0;
wi(melt<-12)=0;

% Find ht


% Transition stratiform
tsInds=find(~isnan(temp) & isnan(wr) & echo<20);

end