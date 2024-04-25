% Width test

clear all;
close all;

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/airMotion/cases/widthCorr/';

sd1=0.54;
mu=7;
widthC=0.5;

x=-1:0.01:15;
sampleNum=length(x);

% Gauss
y1=1/(sd1*sqrt(2*pi)).*exp(-0.5.*((x-mu)/sd1).^2);
yWC=exp(-0.5.*((x-mu)/widthC).^2);

%% IFFT

ifftY=ifft(y1,[],2);
ifftYC=ifft(yWC,[],2);

ifftYC=ifftYC+10^(-7);

%% Correction

yC=ifftY./abs(ifftYC);
%yC(21:end-19)=0;

%% FFT back
fftY=fft(yC,[],2);
y2=real(fftY);
y2=y2./sqrt(sampleNum).*widthC*sqrt(2*pi./((x(end)-x(1)).*(x(2)-x(1))));

%% WIDTH
% VEL
vel1=sum(y1.*x)/sum(y1);
% WIDTH
width1=(sum(y1.*(x-vel1).^2)./sum(y1)).^0.5;

% VEL
vel2=sum(y2.*x)/sum(y2);
% WIDTH
width2=(sum(y2.*(x-vel2).^2)./sum(y2)).^0.5;

widthSmall=sqrt(sd1^2-widthC.^2);

widthDiffShould=widthSmall-sd1;

widthDiffIs=width2-sd1;

%% Plot

f1=figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'renderer','painters');
t=tiledlayout(3,1);

s1=nexttile(1);
plot(x,y1,'-b','LineWidth',1.5);
hold on
plot(x,yWC,'-r','LineWidth',1.5);
xlim([x(1),x(end)])
grid on
box on

legend('Original','Aircraft width')
title(['Width in: ',num2str(sd1),'; Correction factor: ',num2str(widthC)]);

s2=nexttile(2);
hold on
plot(abs(ifftY),'-b','LineWidth',1.5)
plot(abs(yC),'-r','LineWidth',1.5)
xlim([1,length(yC)])
grid on
box on

legend('Original','Corrected')
title(['Width in: ',num2str(sd1),'; Correction factor: ',num2str(widthC)]);

s3=nexttile(3);
plot(x,y1,'-b','LineWidth',1.5);
hold on
plot(x,y2,'-r','LineWidth',1.5);
xlim([x(1),x(end)])
grid on
box on

legend(['Width calc orig: ',num2str(width1)],['Width corrected: ',num2str(width2)])
title(['Width in: ',num2str(sd1),'; Corrected in: ',num2str(widthSmall)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'testWidthCorr_sig_',num2str(sd1,2),'_sigC_',num2str(widthC,1),'.png'],'-dpng','-r0');