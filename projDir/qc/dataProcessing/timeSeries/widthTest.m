% Width test

clear all;
close all;

sd=2;
mu=5;

x=-3:0.1:15;

y=1/(sd*sqrt(2*pi)).*exp(-0.5.*(x-mu).^2./sd^2);

plot(x,y);

% VEL
vel=sum(y.*x)/sum(y);

% WIDTH
width=(sum(y.*(x-vel).^2)./sum(y)).^0.5;