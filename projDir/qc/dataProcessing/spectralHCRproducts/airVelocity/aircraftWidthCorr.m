% Correct spectrum for broadening due to aircraft motion
function y2=aircraftWidthCorr(signalIn,corrFactor,xVel,noise_v,noiseThresh,powNoNoise)
noiseLinV=10.^(noiseThresh./10);
sigInLin=10.^(signalIn./10)-noiseLinV;
sigInLinNoNoise=sigInLin;
sigInLinNoNoise(isnan(powNoNoise))=nan;

% VEL
vel1=sum(sigInLinNoNoise.*xVel,'omitmissing')/sum(sigInLinNoNoise,'omitmissing');
% WIDTH
width1=(sum(sigInLinNoNoise.*(xVel-vel1).^2,'omitmissing')./sum(sigInLinNoNoise,'omitmissing')).^0.5;

% VEL
vel1withNoise=sum(sigInLin.*xVel,'omitmissing')/sum(sigInLin,'omitmissing');
% WIDTH
width1withNoise=(sum(sigInLin.*(xVel-vel1withNoise).^2,'omitmissing')./sum(sigInLin,'omitmissing')).^0.5;

%signalIn=signalIn-noiseThresh;
%signalIn(isnan(signalIn))=noiseThresh;
%signalIn(signalIn<0)=0;

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
notZero=6;

yC=ifftY./abs(ifftYC);
yC(notZero:end-notZero+2)=0;

%% Original uncorrected but filtered

ifftYuf=ifftY;
ifftYuf(notZero:end-notZero+2)=0;

%% FFT back
fftY=fft(yC,[],2);
y2=real(fftY);
y2=y2./sqrt(sampleNum).*corrFactor*sqrt(2*pi./((xVel(end)-xVel(1)).*(xVel(2)-xVel(1))));
%y2=y2(sampleNum*padEachSide+1:padEachSide*sampleNum+sampleNum);

y3=fft(ifftYuf,[],2);

%% Linear space
y2Lin=10.^(y2./10)-noiseLinV;
y2LinNoNoise=y2Lin;
y2LinNoNoise(isnan(powNoNoise))=nan;

y3Lin=10.^(y3./10)-noiseLinV;
y3LinNoNoise=y3Lin;
y3LinNoNoise(isnan(powNoNoise))=nan;

%% WIDTH

% VEL
vel2=sum(y2LinNoNoise.*xVel,'omitmissing')/sum(y2LinNoNoise,'omitmissing');
% WIDTH
width2=(sum(y2LinNoNoise.*(xVel-vel2).^2,'omitmissing')./sum(y2LinNoNoise,'omitmissing')).^0.5;

% VEL
vel2withNoise=sum(y2Lin.*xVel)/sum(y2Lin);
% WIDTH
width2withNoise=(sum(y2Lin.*(xVel-vel2withNoise).^2)./sum(y2Lin)).^0.5;

% VEL
vel3=sum(y3LinNoNoise.*xVel,'omitmissing')/sum(y3LinNoNoise,'omitmissing');
% WIDTH
width3=(sum(y3LinNoNoise.*(xVel-vel3).^2,'omitmissing')./sum(y3LinNoNoise,'omitmissing')).^0.5;

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
plot(cat(2,abs(ifftY(1:50)),abs(ifftY(950:end))),'-b','LineWidth',1.5)
plot(cat(2,abs(yC(1:50)),abs(yC(950:end))),'-r','LineWidth',1.5)
%xlim([1,length(yC)])
grid on
box on

legend('Original','Corrected')

s3=nexttile(3);
plot(xVel,signalIn,'-b','LineWidth',1.5);
hold on
plot(xVel,y2,'-r','LineWidth',1.5);
plot(xVel,y3,'-g','LineWidth',1.5);
plot([xVel(1),xVel(end)],[noiseThresh,noiseThresh],'-c','LineWidth',1.5)
xlim([xVel(1),xVel(end)])
grid on
box on

legend(['Width calc orig: ',num2str(width1)],['Width corrected: ',num2str(width2)],['Width orig filt: ',num2str(real(width3))])
title(['Corrected width analytic: ',num2str(widthSmall,3),'; Difference (analytic-numeric): ',num2str(wDiff,2)]);

% set(gcf,'PaperPositionMode','auto')
% print(f1,[figdir,'testWidthCorr_sig_',num2str(999,2),'_sigC_',num2str(corrFactor,1),'.png'],'-dpng','-r0');
end