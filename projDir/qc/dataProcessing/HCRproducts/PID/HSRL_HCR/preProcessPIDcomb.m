function data = preProcessPIDcomb(data,convThresh)
% Remove fields where they are not suitable

%Reverse up pointing vel
data.HCR_VEL=-data.HCR_VEL;

% Sometimes there are artefacts in the first two gates of LDR
thirdLDR=find(isnan(data.HCR_LDR(20,:)));
data.HCR_LDR(19,thirdLDR)=nan;
data.HCR_LDR(18,thirdLDR)=nan;

% Remove LDR below melting layer
data.HCR_LDR(data.HCR_MELTING_LAYER==9 & data.HCR_LDR<-20)=nan;

% Remove data with too much VELTEXT above melting layer
data.HCR_VEL((data.VELTEXT>convThresh | isnan(data.VELTEXT)) & data.HCR_MELTING_LAYER==21)=nan;

% Remove temperature data above the melting layer and no LDR
data.TEMP(data.HCR_MELTING_LAYER==21 & isnan(data.HCR_LDR) & isnan(data.HSRL_Particle_Linear_Depolarization_Ratio))=nan;

% Mask HSRL data
data.HSRL_Aerosol_Backscatter_Coefficient(data.HSRL_Aerosol_Backscatter_Coefficient<9.9e-9)=nan;
data.HSRL_Particle_Linear_Depolarization_Ratio(data.HSRL_Aerosol_Backscatter_Coefficient<9.9e-9)=nan;

% Remove lidar below melting layer
data.HSRL_Aerosol_Backscatter_Coefficient(data.HCR_MELTING_LAYER==9)=nan;
data.HSRL_Particle_Linear_Depolarization_Ratio(data.HCR_MELTING_LAYER==9)=nan;
end