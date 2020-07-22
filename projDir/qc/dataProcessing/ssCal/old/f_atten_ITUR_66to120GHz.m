%Calculates attenuation per kilometer after ITU method. It is only valid
%for altitudes up to ~10km and for radar frequencies between 66GHz and 120
%GHz. The naming conventions follow the reference below and the equations
%have been split up into sub terms.

%Radiocommunication Sector of International Telecommunication Union. 
%Recommendation ITU-R P.676-10: Attenuation by atmospheric gases 2013.

%Author: Ulrike Romatschke romatsch@ucar.edu
%Last modified: 20170106

%Input:
%f is the radar freqyency in GHz
%ptot is the pressure in mb
%t is the temperature in C
%rh is the relative humidity in %
%es is the saturation pressure

%Output:
%gamma is the total attenuation
%gamma0 is the dry attenuation
%gammaW is he wet attenuation

function [gamma, gamma0, gammaW]= f_atten_ITUR_66to120GHz(f,ptot,t,rh,es)

%calculate absolute humidity rho from relative humidity rh
rho=(es.*rh.*2.167)./(273.15+t); %water vapour density in g/m^3 (absolute humidity)

rp = ptot./1013;
rt = 288./(273+t);

a4=-0.0112;
b4=0.0092;
c4=-0.1033;
d4=-0.0009;

a5=0.2705;
b5=-2.7192;
c5=-0.3016;
d5=-4.1033;

a6=0.2445;
b6=-5.9191;
c6=0.0422;
d6=-8.0719;

a7=-0.1833;
b7=6.5589;
c7=-0.2402;
d7=6.131;

xi4 = rp.^a4.*rt.^b4.*exp(c4.*(1-rp)+d4.*(1-rt));
xi5 = rp.^a5.*rt.^b5.*exp(c5.*(1-rp)+d5.*(1-rt));
xi6 = rp.^a6.*rt.^b6.*exp(c6.*(1-rp)+d6.*(1-rt));
xi7 = rp.^a7.*rt.^b7.*exp(c7.*(1-rp)+d7.*(1-rt));

A=3.02.*10.^(-4).*rt.^3.5;
B=(0.283.*rt.^3.8)./((f-118.75).^2+2.91.*rp.^2.*rt.^1.6);
C=(0.502.*xi6.*(1-0.0163.*xi7.*(f-66)))./((f-66).^(1.4346.*xi4)+1.15.*xi5);
D=f.^2.*rp.^2.*10.^(-3);

gamma0=(A+B+C).*D;

fi1=22;
fi2=557;
fi3=752;
fi4=1780;

g1=1+((f-fi1)./(f+fi1)).^2;
g2=1+((f-fi2)./(f+fi2)).^2;
g3=1+((f-fi3)./(f+fi3)).^2;
g4=1+((f-fi4)./(f+fi4)).^2;

eta1=0.955.*rp.*rt.^0.68+0.006.*rho;

eta2=0.735.*rp.*rt.^0.5+0.0353.*rt.^4.*rho;

E=(3.98.*eta1.*exp(2.23.*(1-rt)))./((f-22.235).^2+9.42.*eta1.^2).*g1;
F=(11.96.*eta1.*exp(0.7.*(1-rt)))./((f-183.31).^2+11.14.*eta1.^2);
G=(0.081.*eta1.*exp(6.44.*(1-rt)))./((f-321.226).^2+6.29.*eta1.^2);
H=(3.66.*eta1.*exp(1.6.*(1-rt)))./((f-325.153).^2+9.22.*eta1.^2);
I=(25.37.*eta1.*exp(1.09.*(1-rt)))./(f-380).^2;
J=(17.4.*eta1.*exp(1.46.*(1-rt)))./(f-448).^2;
K=(844.6.*eta1.*exp(0.17.*(1-rt)))./(f-557).^2.*g2;
L=(290.*eta1.*exp(0.41.*(1-rt)))./(f-752).^2.*g3;
M=(8.3328.*10.^4.*eta2.*exp(0.99.*(1-rt)))./(f-1780).^2.*g4;
N=f.^2.*rt.^2.5.*rho.*10.^(-4);

gammaW=(E+F+G+H+I+J+K+L+M).*N;

gamma=gamma0+gammaW;
