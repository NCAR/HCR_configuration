function [sig0measAtt,surfFlag,refSig0,refFlag,sig0model,piaGas2,piaHydromet2]=getRefAtten(data)
% Calculate surface reference attenuation
%% One way and two way gaseous attenuation

disp('Calculating gaseous attenuation ...');

[~,gasAttCloud,~,~]=get_gas_atten(data);
piaGas2=2*gasAttCloud';

%% Calculate sigma0 from model and from reflectivity

disp('Calculating sig0 ...');

% Find ocean surface gate
[linInd,maxGate,rangeToSurf]=hcrSurfInds(data);

% Measured sig0 from surface reflectivity
data.surfRefl=data.DBZ(linInd);
sig0measured=calc_sig0_surfRefl(data);

sig0measAtt=sig0measured(linInd)+piaGas2;
sig0measAtt(data.elevation>-85)=nan;

% sig0 from models
sig0modelAll= calc_sig0_model(data);
%sig0model=sig0modelAll(2,:); % Freilich Vanhoff
%sig0model=sig0modelAll(5,:); % Wu
sig0model=sig0modelAll(8,:); % Cox Munk

clear sig0modelAll
%% Create ocean surface mask
% 0 extinct or not usable
% 1 cloud
% 2 clear air

[surfFlag1,atmFrac]=makeSurfFlag(data,linInd);

clear gasAttCloudMat
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