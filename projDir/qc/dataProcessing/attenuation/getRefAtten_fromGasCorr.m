function [sig0measAtt,surfFlag,refSig0,refFlag,sig0model,piaHydromet2]=getRefAtten_fromGasCorr(data)

%% Calculate sigma0 from model and from reflectivity

disp('Calculating sig0 ...');

% Find ocean surface gate
[linInd,maxGate,rangeToSurf]=hcrSurfInds(data);

% Measured sig0 from surface reflectivity
data.surfRefl=data.DBZcorrGas(linInd);
sig0measAtt=calc_sig0_surfGasCorr(data);

% sig0 from models
sig0modelAll=calc_sig0_model(data);
%sig0model=sig0modelAll(2,:); % Freilich Vanhoff
%sig0model=sig0modelAll(5,:); % Wu
sig0model=sig0modelAll(8,:); % Cox Munk

clear sig0modelAll
%% Create ocean surface mask
% 0 extinct
% 1 cloud
% 2 clear air

[surfFlag1,~]=makeSurfFlag(data,linInd);

%% Create field with reference sig0
% RefFlag
% 1 clear air
% 2 interpolated
% 3 model

[refSig0,surfFlag,refFlag]=makeRefSig0(sig0measAtt,sig0model,surfFlag1);

% Find surfFlag values that have previously been clear air but are now
% not
surfFlag(surfFlag1~=0 & surfFlag==0 & any(data.FLAG==1,1))=1;

%% 2 way path integrated attenuation from hydrometeors

piaHydromet2=refSig0-sig0measAtt;
piaHydromet2(surfFlag~=1)=nan;
end