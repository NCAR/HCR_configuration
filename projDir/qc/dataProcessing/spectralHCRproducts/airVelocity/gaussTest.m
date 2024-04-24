% Width test

clear all;
close all;

sd1=1;
mu=3;
widthC=0.5;

x=-1:0.05:15;
sampleNum=length(x);
prt=1.0138e-04;
lambda=0.0032;

y1=1/(sd1*sqrt(2*pi)).*exp(-0.5.*(x-mu).^2./sd1^2);

%% IFFT

ifftY=ifft(y1,[],2);
ifftYs=ifftshift(ifftY,2);

x1=(-pi:2*pi/(sampleNum):pi);
x1=x1(1:end-1);

yC=ifftYs.*exp(1./((sampleNum.*x1).^2).*(widthC*4*pi.*prt./lambda)^2);

figure
hold on
plot(abs(ifftYs))
plot(abs(yC))


%% FFT
fftY=fft(yC,[],2);
y2=abs(fftY);

figure
plot(x,y1);
hold on
plot(x,y2);
grid on

%% WIDTH
% VEL
vel1=sum(y1.*x)/sum(y1);
% WIDTH
width1=(sum(y1.*(x-vel1).^2)./sum(y1)).^0.5;

% VEL
vel2=sum(y2.*x)/sum(y2);
% WIDTH
width2=(sum(y2.*(x-vel1).^2)./sum(y2)).^0.5;

widthSmall=sqrt(sd1^2-widthC.^2);

widthDiffShould=widthSmall-sd1;

widthDiffIs=width2-sd1;