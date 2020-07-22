% Calculate refractive index of sea water after Meissner and Wentz, 2004
function [n]= f_refractiveIndex_MeissnerWentz(T,S,frqIn)

% T= sea surface  temperature in C
% S= salinity in pro mille (mean for ocean is 35)

%frq= frequency in Hz

frq=frqIn/1e+9; % Convert to GHz

eps0=17.9751; % Vacuum electric permittivity in GHz m/S

% sigma

sig35=2.903602+8.607e-2*T+4.738817e-4*T.^2- ...
    2.991e-6*T.^3+4.3047e-9*T.^4;

R15=S*(37.5109+5.45216*S+1.4409e-2*S^2)/(1004.75+182.283*S+S^2);

alpha0=(6.9431+3.2841*S-9.9486e-2*S^2)/(84.850+69.024*S+S^2);
alpha1=49.843-0.2276*S+0.198e-2*S^2;
RToverR15=1+alpha0.*(T-15)./(alpha1+T);

sig=sig35.*R15.*RToverR15;

% Values at S=0

% a values
a0=5.7230e+0;
a1=2.2379e-2;
a2=-7.1237e-4;
a3=5.0478e+0;
a4=-7.0315e-2;
a5=6.0059e-4;
a6=3.6143e+0;
a7=2.8841e-2;
a8=1.3652e-1;
a9=1.4825e-3;
a10=2.4166e-4;

eps_s0=(3.70886e+4-8.2168e+1.*T)./(4.21854e+2+T);
eps1_s0=a0+a1*T+a2*T.^2;
nu1_s0=(45+T)./(a3+a4.*T+a5.*T.^2);
epsInf_s0=a6+a7*T;
nu2_s0=(45+T)./(a8+a9.*T+a10.*T.^2);

% With salinity

% b values

b0=-3.56417e-3;
b1=4.74868e-6;
b2=1.15574e-5;
b3=2.39357e-3;
b4=-3.1353e-5;
b5=2.52477e-7;
b6=-6.28908e-3;
b7=1.76032e-4;
b8=-9.22144e-5;
b9=-1.99723e-2;
b10=1.81176e-4;
b11=-2.04265e-3;
b12=1.57883e-4;

epsS=eps_s0.*exp(b0.*S+b1.*S^2+b2.*T.*S);
nu1=nu1_s0.*(1+S.*(b3+b4.*T+b5.*T.^2));
eps1=eps1_s0.*exp(b6.*S+b7.*S^2+b8.*T.*S);
nu2=nu2_s0.*(1+S*(b9+b10.*T));
epsInf=epsInf_s0.*(1+S*(b11+b12.*T));

% epsA=(eps_s0-eps1_s0)/(1+(i*frq/nu1_s0));
% epsB=(eps1_s0-epsInf_s0)/(1+(i*frq/nu2_s0));
% epsC=epsInf_s0-i*sig/(2*pi*eps0*frq);

epsAS=(epsS-eps1)./(1+(i.*frq./nu1));
epsBS=(eps1-epsInf)./(1+(i.*frq./nu2));
epsCS=epsInf-i.*sig./(2.*pi.*eps0.*frq);

%epsWater=epsA+epsB+epsC
epsSalt=epsAS+epsBS+epsCS;

n=sqrt(epsSalt);
end
