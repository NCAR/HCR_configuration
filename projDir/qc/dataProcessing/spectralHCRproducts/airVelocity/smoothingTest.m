% Correct spectrum for broadening due to aircraft motion
function [err]=smoothingTest(signalIn,corrFactor,xVel,noiseThresh,sampleNum,err,inds1,inds2,figdir)
sigWidthCorr=[];

%inds1=1:2:length(signalIn);
signalIn1=signalIn(inds1);

%inds2=2:2:length(signalIn);
signalIn2=signalIn(inds2);

if length(signalIn1)~=length(signalIn2)
    inds1(end)=[];
    signalIn1(end)=[];
end

% noiseLinV=10.^(noiseThresh./10);
% sigInLin=10.^(signalIn./10)-noiseLinV;

% VEL
% meanVel=sum(sigInLin.*xVel,'omitmissing')/sum(sigInLin,'omitmissing');

% Gaussian fit of correction signal
% yWC=exp(-0.5.*((xVel-meanVel)/corrFactor).^2);

%% IFFT

ifftY=ifft(signalIn,[],2);
ifftY1=ifft(signalIn1,[],2);
ifftY2=ifft(signalIn2,[],2);

% ifftYC=ifft(yWC,[],2);
% ifftYC=ifftYC+10^(-7);

%% Correction
notZero=2:250;

% yC=ifftY./abs(ifftYC);
% yC(notZero:end-notZero+2)=0;

%% Original uncorrected but filtered

ifftYuf=repmat(ifftY,length(notZero),1);
ifftYuf1=repmat(ifftY1,length(notZero),1);
ifftYuf2=repmat(ifftY2,length(notZero),1);
for ii=1:length(notZero)
    ifftYuf(ii,notZero(ii):end-notZero(ii)+2)=0;
    ifftYuf1(ii,notZero(ii):end-notZero(ii)+2)=0;
    ifftYuf2(ii,notZero(ii):end-notZero(ii)+2)=0;
end

%% FFT back
% fftY=fft(yC,[],2);
% sigWidthCorr=real(fftY);
% sigWidthCorr=sigWidthCorr./sqrt(sampleNum).*corrFactor*sqrt(2*pi./((xVel(end)-xVel(1)).*(xVel(2)-xVel(1))));

%sigFiltered=real(fft(ifftYuf,[],2));
sigFiltered1=real(fft(ifftYuf1,[],2));
sigFiltered2=real(fft(ifftYuf2,[],2));

%% Error

% Interpolate
sigFiltered1i=(interp1(inds1,sigFiltered1',inds2,'linear','extrap'))';
sigFiltered2i=(interp1(inds2,sigFiltered2',inds1,'linear','extrap'))';

err12=rmse(sigFiltered1i,signalIn2,2);
err21=rmse(sigFiltered2i,signalIn1,2);

errCat=cat(2,err12,err21);

try
    err=cat(2,err,errCat);
end

%% Plots

% close all
% figure('Position',[200 500 1000 600],'DefaultAxesFontSize',12,'renderer','painters')
% plot(notZero,err12)
% hold on
% plot(notZero,err21)
% hold off
% ylim([0,6])
% 
% grid on
% box on
% 
% filterAt=9;
% 
% %close all
% f1=figure('Position',[200 500 1000 900],'DefaultAxesFontSize',12,'renderer','painters');
% t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
% s1=nexttile(1);
% 
% grid on
% box on
% hold on
% plot(xVel,signalIn,'-b','LineWidth',2)
% plot(xVel(inds1),signalIn1,'-c','LineWidth',1)
% plot(xVel(inds2),signalIn2,'-k','LineWidth',1)
% xlim([xVel(1),xVel(end)]);
% 
% legend('Original signal','Split signal 1','Split signal 2')
% 
% s2=nexttile(2);
% 
% hold on
% plot(xVel(inds1),sigFiltered1(filterAt,:),'-c','LineWidth',1)
% plot(xVel(inds2),sigFiltered2(filterAt,:),'-k','LineWidth',1)
% xlim([xVel(1),xVel(end)]);
% 
% legend('Filtered split signal 1','Filtered split signal 2')
% 
% grid on
% box on
% 
% %plot(xVel,sigFiltered(ii,:),'-g','LineWidth',2)
% % plot(xVel,sigWidthCorr,'-r','LineWidth',2)
% %plot(xVel,ones(1,length(xVel))*noiseThresh,'-c','LineWidth',2)
% hold off
% 
% set(gcf,'PaperPositionMode','auto')
% print(f1,[figdir,'zero',num2str(filterAt),'_random'],'-dpng','-r0');

end