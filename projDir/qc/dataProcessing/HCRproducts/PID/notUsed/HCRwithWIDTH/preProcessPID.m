function data = preProcessPID(data,convThresh,widthThresh)
% Remove fields where they are not suitable

% Censor spectrum width
data.WIDTH(data.SNR<5)=nan;

%Reverse up pointing vel
data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

% Fix melting layer
data.MELTING_LAYER(~isnan(data.MELTING_LAYER) & data.MELTING_LAYER<20)=10;
data.MELTING_LAYER(~isnan(data.MELTING_LAYER) & data.MELTING_LAYER>=20)=20;

% Sometimes there are artefacts in the first two gates of LDR
thirdLDR=find(isnan(data.LDR(20,:)));
data.LDR(19,thirdLDR)=nan;
data.LDR(18,thirdLDR)=nan;

% Sometimes there are artefacts in the first gate of WIDTH
secondWIDTH=data.WIDTH(19,:);
firstWIDTH=data.WIDTH(18,:);

outWIDTH=find((firstWIDTH>0.5 & firstWIDTH-secondWIDTH>0.4) | isnan(secondWIDTH));
data.WIDTH(18,outWIDTH)=nan;

% Remove data with too much VELTEXT below above layer
data.VEL_MASKED((data.VELTEXT>convThresh | isnan(data.VELTEXT)) & data.MELTING_LAYER==20)=nan;
data.WIDTH((data.VELTEXT>convThresh | isnan(data.VELTEXT)) & data.WIDTH>widthThresh & data.MELTING_LAYER==20)=nan;

end