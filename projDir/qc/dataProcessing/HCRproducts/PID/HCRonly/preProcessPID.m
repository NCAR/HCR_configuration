function data = preProcessPID(data,convThresh,widthThresh)
% Remove fields where they are not suitable

% Censor spectrum width
data.WIDTH(data.SNR<5)=nan;

%Reverse up pointing vel
data.VEL_MASKED(data.elevation>0)=-data.VEL_MASKED(data.elevation>0);

% Fix melting layer
data.MELTING_LAYER(~isnan(data.MELTING_LAYER) & data.MELTING_LAYER<20)=10;
data.MELTING_LAYER(~isnan(data.MELTING_LAYER) & data.MELTING_LAYER>=20)=20;

% Sometimes there are artefacts in the first two gates of LDR
thirdLDR=find(isnan(data.LDR(20,:)));
data.LDR(19,thirdLDR)=nan;
data.LDR(18,thirdLDR)=nan;

% Remove data with too much convectivity
data.VEL_MASKED(data.CONVECTIVITY>convThresh | isnan(data.CONVECTIVITY))=nan;
data.WIDTH((data.CONVECTIVITY>convThresh | isnan(data.CONVECTIVITY)) & data.WIDTH>widthThresh)=nan;

end