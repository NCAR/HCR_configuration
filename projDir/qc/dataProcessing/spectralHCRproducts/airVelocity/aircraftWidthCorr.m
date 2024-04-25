% Correct spectrum for broadening due to aircraft motion
function y2=aircraftWidthCorr(signalIn,corrFactor,xVel,noise_v,noiseThresh)
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

% % Pad incoming signal
% padEachSide=7;
% sigInPadded=cat(2,zeros(1,padEachSide*sampleNum),signalIn,zeros(1,padEachSide*sampleNum));

% Gaussian fit of correction signal
yWC=exp(-0.5.*((xVel-vel1)/corrFactor).^2);
%yWCpadded=cat(2,zeros(1,padEachSide*sampleNum),yWC,zeros(1,padEachSide*sampleNum));
%% IFFT

ifftY=ifft(signalIn,[],2);

ifftYC=ifft(yWC,[],2);
ifftYC=ifftYC+10^(-7);

%% Correction
yC=ifftY./abs(ifftYC);
yC(11:end-9)=0;

%% FFT back
fftY=fft(yC,[],2);
y2=real(fftY);
y2=y2./sqrt(sampleNum).*corrFactor*sqrt(2*pi./((xVel(end)-xVel(1)).*(xVel(2)-xVel(1))));
%y2=y2(sampleNum*padEachSide+1:padEachSide*sampleNum+sampleNum);

y2Lin=10.^(y2./10)-noiseLinV;
%% WIDTH

% VEL
vel2=sum(y2Lin.*xVel)/sum(y2Lin);
% WIDTH
width2=(sum(y2Lin.*(xVel-vel2).^2)./sum(y2Lin)).^0.5;

widthSmall=sqrt(width1^2-corrFactor.^2);

wDiff=widthSmall-width2;

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
title(['Correction factor: ',num2str(corrFactor)]);

s2=nexttile(2);
hold on
plot(abs(ifftY),'-b','LineWidth',1.5)
plot(abs(yC),'-r','LineWidth',1.5)
xlim([1,length(yC)])
grid on
box on

legend('Original','Corrected')

s3=nexttile(3);
plot(xVel,signalIn,'-b','LineWidth',1.5);
hold on
plot(xVel,y2,'-r','LineWidth',1.5);
xlim([xVel(1),xVel(end)])
grid on
box on

legend(['Width calc orig: ',num2str(width1)],['Width corrected: ',num2str(width2)])
title(['Corrected width analytic: ',num2str(widthSmall,3),'; Difference (analytic-numeric): ',num2str(wDiff,2)]);

% set(gcf,'PaperPositionMode','auto')
% print(f1,[figdir,'testWidthCorr_sig_',num2str(999,2),'_sigC_',num2str(corrFactor,1),'.png'],'-dpng','-r0');
end