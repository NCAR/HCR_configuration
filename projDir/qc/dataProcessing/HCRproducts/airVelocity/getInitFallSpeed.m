function [vr,vi]=getInitFallSpeed(data)
% Calculate initial fall speed

% V rain
% https://doi-org.cuucar.idm.oclc.org/10.1175/1520-0426(1994)011<1656:EOVIPS>2.0.CO;2
a=3.5;
b=0.084;

zlin=10.^(data.DBZ_MASKED./10);

% V rain
% Eqs. 8 and 9 https://doi-org.cuucar.idm.oclc.org/10.1175/1520-0426(1994)011<1656:EOVIPS>2.0.CO;2

mu=0;
N0=6.4*10^4*exp(3.2*mu);
innerPar=(zlin./(N0.*10^6.*gamma(7+mu))).^(1/(7+mu));
outerPar=(1+6.*innerPar).^-(7+mu);
vr=9.65-10.3.*outerPar;

% V ice
% https://doi.org/10.1175/2009JAS3132.1
vi=nan(size(data.DBZ_MASKED));

vi(data.DBZ_MASKED<30)=1.4041+0.00237.*data.DBZ_MASKED(data.DBZ_MASKED<30);
vi(data.DBZ_MASKED>=30)=-3.4+0.1864.*data.DBZ_MASKED(data.DBZ_MASKED>=30);

% Adjust for air density
% Calculate air density
% (https://www.omnicalculator.com/physics/air-density)
Rd=287.05; % Gas constantfor dry air in J/(kg K)
Rv=461.5; % Gas constant for water vapor in J/(kg K)

satP=6.1094.*exp(17.625.*data.TEMP./(data.TEMP+243.04));
vaporP=satP.*data.RH./100;
pDry=data.PRESS-vaporP;

rho=pDry.*100./(Rd.*(data.TEMP+273.15))+vaporP.*100./(Rv.*(data.TEMP+273.15));

surfIndInit=data.asl;
surfIndInit(isnan(data.TEMP))=nan;

[minAsl surfInd]=min(surfIndInit,[],1,'omitnan');

matInd=sub2ind(size(data.asl),surfInd,1:length(data.time));

rhoSurf=rho(matInd);

vr=vr.*(rhoSurf./rho).^0.45;
vi=vi.*(rhoSurf./rho).^0.45;

end