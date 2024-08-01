% Width test

clear all;
close all;

sd1=0.8;
mu=6;
widthC=0.5;

x=-1:0.05:15;
sampleNum=length(x);
% prt=1.0138e-04;
% lambda=0.0032;

% Gauss
y1=1/(sd1*sqrt(2*pi)).*exp(-0.5.*((x-mu)/sd1).^2);
y1padded=cat(2,zeros(1,sampleNum),y1,zeros(1,sampleNum));

%% IFFT

ifftY=ifft(y1padded,[],2);
ifftYs=ifftshift(ifftY,2);

%% Correction
L=(x(end)-x(1))*3;
K=(0:length(ifftYs)-1)-ceil(length(ifftYs)/2);

expC=exp((2*pi.*K./L).^2.*(widthC)^2);
expThresh=10^13;
expC(expC>expThresh)=expThresh;

yC=ifftYs.*expC;

%% FFT back
fftY=fft(yC,[],2);
y2=abs(fftY);
y2=y2(sampleNum+1:2*sampleNum);

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

%% Plot

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'renderer','painters');
t=tiledlayout(2,1);

s1=nexttile(1);
hold on
plot(K,abs(ifftYs),'-b','LineWidth',1.5)
plot(K,abs(yC),'-r','LineWidth',1.5)
xlim([K(1),K(end)])
grid on
box on

yyaxis right
plot(K,expC,'-g','LineWidth',1.5)
ylim([0,100])

legend('Original','Corrected','Correction')
title(['Width in: ',num2str(sd1),'; Correction factor: ',num2str(widthC)]);

s2=nexttile(2);
plot(x,y1,'-b','LineWidth',1.5);
hold on
plot(x,y2,'-r','LineWidth',1.5);
xlim([x(1),x(end)])
grid on
box on

legend(['Width calc orig: ',num2str(width1)],['Width corrected: ',num2str(width2)])
title(['Width in: ',num2str(sd1),'; Corrected in: ',num2str(widthSmall)]);