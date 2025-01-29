% Width test

clear all;
close all;

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/forReviewers/';

numIter=100000;

sd1=1;
mu=7;
widthC=0.5;

lambda=0.0032;
prt=1.013759974739514e-04;

x=-1:0.05:15;
sampleNum=length(x);

% Gauss
y1=1/(sd1*sqrt(2*pi)).*exp(-0.5.*((x-mu)/sd1).^2);

%% WIDTH SPECTRAL
% VEL
vel1=sum(y1.*x)/sum(y1);
% WIDTH
width1=(sum(y1.*(x-vel1).^2)./sum(y1)).^0.5;

%% WIDTH PulsePair

y1shift=fftshift(y1,2);
y1Lin=10.^(y1shift./10);

cIQ1=ifft(y1Lin,[],2).*sqrt(size(y1,2));

R1v1=mean(cIQ1(:,1:end-1).*conj(cIQ1(:,2:end)),2);
R2v1=mean(cIQ1(:,1:end-2).*conj(cIQ1(:,3:end)),2);

widthPP1=2/(pi.*sqrt(6))*(x(end)-x(1))/2*abs(log(abs(R1v1./R2v1))).^0.5;

widthDiffall=[];

for ii=1:numIter
    % Noise
    noiseR=randn(size(y1))*sqrt(0.01);

    y2=y1+noiseR;

    % With noise
    y2shift=fftshift(y2,2);
    y2Lin=10.^(y2shift./10);

    cIQ2=ifft(y2Lin,[],2).*sqrt(size(y2,2));

    R1v2=mean(cIQ2(:,1:end-1).*conj(cIQ2(:,2:end)),2);
    R2v2=mean(cIQ2(:,1:end-2).*conj(cIQ2(:,3:end)),2);

    widthPP2=2/(pi.*sqrt(6))*(x(end)-x(1))/2*abs(log(abs(R1v2./R2v2))).^0.5;

    widthDiffall=cat(1,widthDiffall,widthPP2-sd1);
end

meanDiff=mean(widthDiffall);

%% Plot

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'renderer','painters');
t=tiledlayout(2,1);

s1=nexttile(1);
hold on
plot(x,y1,'-b','LineWidth',1.5);
plot(x,y2,'-g','LineWidth',1);
xlim([x(1),x(end)])
grid on
box on

legend(['Width spectral: ',num2str(width1)],['Width pulse pair: ',num2str(widthPP2)])
title(['Width in: ',num2str(sd1)]);

s2=nexttile(2);
edgesIn=-1:0.1:1;
hc=histcounts(widthDiffall,edgesIn);
bar(edgesIn(1:end-1)+0.05,hc,1);

title(['Mean diff: ',num2str(meanDiff),', ',num2str(numIter),' iterations'])

print([figdir,'widthNoise.png'],'-dpng','-r0');