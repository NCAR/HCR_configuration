% Correct spectrum for broadening due to aircraft motion
function y2=aircraftWidthCorr_firstTry(signalIn,corrFactor,xVel,noise_v,noiseThresh)
signalIn=signalIn-noiseThresh;
signalIn(isnan(signalIn))=0;
signalIn(signalIn<0)=0;

noiseLinV=10.^(noise_v./10);
sigInLin=10.^(signalIn./10)-noiseLinV;

% VEL
vel1=sum(sigInLin.*xVel)/sum(sigInLin);
% WIDTH
width1=(sum(sigInLin.*(xVel-vel1).^2)./sum(sigInLin)).^0.5;

sampleNum=length(xVel);

% Pad incoming signal
padEachSide=7;
sigInPadded=cat(2,zeros(1,padEachSide*sampleNum),signalIn,zeros(1,padEachSide*sampleNum));

% Gaussian fit of correction signal
yWC=1/(corrFactor*sqrt(2*pi)).*exp(-0.5.*((xVel-vel1)/corrFactor).^2);
yWCpadded=cat(2,zeros(1,padEachSide*sampleNum),yWC,zeros(1,padEachSide*sampleNum));
%% IFFT

ifftY=ifft(sigInPadded,[],2);
ifftYs=ifftshift(ifftY,2);

ifftYC=ifft(yWCpadded,[],2);
ifftYsC=ifftshift(ifftYC,2);
ifftYsC=ifftYsC+10^(-4);

%% Correction

yC=ifftYs./abs(ifftYsC);

%% FFT back
fftY=fft(yC,[],2);
y2=abs(fftY);
y2=y2./(2/corrFactor)*sqrt(xVel(2)-xVel(1));
y2=y2(sampleNum*padEachSide+1:padEachSide*sampleNum+sampleNum);

y2Lin=10.^(y2./10)-noiseLinV;
%% WIDTH

% VEL
vel2=sum(y2Lin.*xVel)/sum(y2Lin);
% WIDTH
width2=(sum(y2Lin.*(xVel-vel2).^2)./sum(y2Lin)).^0.5;

%% Plot
close all
f1=figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'renderer','painters');
t=tiledlayout(3,1);

s1=nexttile(1);
plot(xVel,signalIn,'-b','LineWidth',1.5);
hold on
plot(xVel,yWC,'-r','LineWidth',1.5);
xlim([xVel(1),xVel(end)])
grid on
box on

legend('Original','Aircraft width')
title(['Width in: ',num2str(999),'; Correction factor: ',num2str(corrFactor)]);

s2=nexttile(2);
hold on
plot(abs(ifftYs),'-b','LineWidth',1.5)
plot(abs(yC),'-r','LineWidth',1.5)
xlim([1,sampleNum*padEachSide*2+sampleNum])
grid on
box on

legend('Original','Corrected')
title(['Width in: ',num2str(999),'; Correction factor: ',num2str(corrFactor)]);

s3=nexttile(3);
plot(xVel,signalIn,'-b','LineWidth',1.5);
hold on
plot(xVel,y2,'-r','LineWidth',1.5);
xlim([xVel(1),xVel(end)])
grid on
box on

legend(['Width calc orig: ',num2str(width1)],['Width corrected: ',num2str(width2)])
%title(['Width in: ',num2str(999),'; Corrected in: ',num2str(widthSmall)]);

% set(gcf,'PaperPositionMode','auto')
% print(f1,[figdir,'testWidthCorr_sig_',num2str(999,2),'_sigC_',num2str(corrFactor,1),'.png'],'-dpng','-r0');
end