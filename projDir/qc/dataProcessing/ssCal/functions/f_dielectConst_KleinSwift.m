% Calculate dielectric constant after Klein and Swift 1977

% !!!!!!!!!!!!!!!!!!!! Not Correct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

T=20;
S=35;

eps_inf=4.9;
alpha=0;

%frq=9.4406e+10;
frq=3.4091e+10;

omega=2*pi*frq;
eps_0=8.854e-12;

epsT=87.134-1.949e-1*T-1.276e-2*T^2+2.491e-4*T^3;
aST=1+1.613e-5*S*T-3.656e-3*S+3.21e-5*S^2-4.232e-7*S^3;
eps_s=epsT*aST;

tauT=1.768e-11-6.086e-13*T+1.104e-14*T^2-8.111e-17*T^3;
bST=1+2.282e-5*S*T-7.638e-4*S-7.76e-6*S^2+1.105e-8*S^3;
tau=tauT*bST;

Delta=25-T;
beta=2.033e-2+1.266e-4*Delta+2.464e-6*Delta^2-S*(1.849e-5-2.551e-7*Delta+2.551e-8*Delta^2);
sig25S=S*(0.182521-1.46192e-3*S+2.09324e-5*S^2-1.28205e-7*S^3);
sig=sig25S*exp(-Delta*beta);

eps=eps_inf+(eps_s-eps_inf)/(1+(i*omega*tau)^(1-alpha))-i*(sig/(omega*eps_0))