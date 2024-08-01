% Correct spectrum for broadening due to aircraft motion
function [err]=smoothingTestTime(signalIn,xVel,signalIn1,xVel1,signalIn2,err,figdir)

%% IFFT

ifftY=ifft(signalIn,[],2);
ifftY1=ifft(signalIn1,[],2);
ifftY2=ifft(signalIn2,[],2);

%% Correction
notZero=2:250;

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

sigFiltered1=real(fft(ifftYuf1,[],2));
sigFiltered2=real(fft(ifftYuf2,[],2));

%% Error

err12=rmse(sigFiltered1,signalIn2,2);
err21=rmse(sigFiltered2,signalIn1,2);

errCat=cat(2,err12,err21);

try
    err=cat(2,err,errCat);
end

%% Plots

close all
figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'renderer','painters')
plot(notZero,err12)
hold on
plot(notZero,err21)
hold off
ylim([0,6])

grid on
box on

filterAt=248;

%close all
f1=figure('Position',[200 500 1000 900],'DefaultAxesFontSize',12,'renderer','painters');
t = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

grid on
box on
hold on
%plot(xVel,signalIn,'-b','LineWidth',2)
plot(xVel1,signalIn1,'-c','LineWidth',1)
plot(xVel1,signalIn2,'-k','LineWidth',1)
%xlim([xVel(1),xVel(end)]);

legend('Original signal','Split signal 1','Split signal 2')

s2=nexttile(2);

hold on
%plot(real(ifftY),'-b','LineWidth',1)
%plot(imag(ifftY),'-r','LineWidth',1)
plot(real(ifftY1),'-k','LineWidth',1)
plot(imag(ifftY1),'-m','LineWidth',1)
plot(real(ifftY2),'-c','LineWidth',1)
plot(imag(ifftY2),'-y','LineWidth',1)

legend('real(orig)','imag(orig)','real(split1)','imag(split1)','real(split2)','imag(split2)');

grid on
box on

s3=nexttile(3);

hold on
plot(xVel1,sigFiltered1(filterAt,:),'-c','LineWidth',1)
plot(xVel1,sigFiltered2(filterAt,:),'-k','LineWidth',1)
%xlim([xVel(1),xVel(end)]);

legend('Filtered split signal 1','Filtered split signal 2')

grid on
box on

%plot(xVel,sigFiltered(ii,:),'-g','LineWidth',2)
% plot(xVel,sigWidthCorr,'-r','LineWidth',2)
%plot(xVel,ones(1,length(xVel))*noiseThresh,'-c','LineWidth',2)
hold off

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'zero',num2str(filterAt),'_timeDB'],'-dpng','-r0');

end