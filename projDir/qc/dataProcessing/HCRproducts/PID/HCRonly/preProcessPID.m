function data = preProcessPID(data,convThresh)
% Remove fields where they are not suitable

%Reverse up pointing vel
data.VEL_MASKED=-data.VEL_MASKED;

% Sometimes there are artefacts in the first two gates of LDR
thirdLDR=find(isnan(data.LDR(20,:)));
data.LDR(19,thirdLDR)=nan;
data.LDR(18,thirdLDR)=nan;

% Remove LDR below melting layer
data.LDR(data.MELTING_LAYER==9 & data.LDR<-20)=nan;

% Remove data with too much VELTEXT above melting layer
data.VEL_MASKED((data.VELTEXT>convThresh | isnan(data.VELTEXT)) & data.MELTING_LAYER==21)=nan;

% Remove temperature data above the melting layer and no LDR
data.TEMP(data.MELTING_LAYER==21 & isnan(data.LDR))=nan;

end