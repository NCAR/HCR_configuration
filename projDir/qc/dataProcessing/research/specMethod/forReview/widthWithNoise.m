% Width test

clear all;
close all;

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/forReviewers/';

numIter=10000;

sd1=[0.3;1;3];
mu=7;

noiseVar=0.05;

x=-1:0.05:15;
x=repmat(x,length(sd1),1);

% Gauss
%y1=1./(sd1.*sqrt(2*pi)).*exp(-0.5.*((x-mu)./sd1).^2);
y1=exp(-0.5.*((x-mu)./sd1).^2);

%% WIDTH SPECTRAL
% VEL
vel1=sum(y1.*x)/sum(y1);
% WIDTH
width1=(sum(y1.*(x-vel1).^2,2)./sum(y1,2)).^0.5;

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
    noiseR=randn(size(y1))*sqrt(noiseVar);

    y2=y1+noiseR;

    % With noise
    y2shift=fftshift(y2,2);
    y2Lin=10.^(y2shift./10);

    cIQ2=ifft(y2Lin,[],2).*sqrt(size(y2,2));

    R1v2=mean(cIQ2(:,1:end-1).*conj(cIQ2(:,2:end)),2);
    R2v2=mean(cIQ2(:,1:end-2).*conj(cIQ2(:,3:end)),2);

    widthPP2=2/(pi.*sqrt(6))*(x(end)-x(1))/2*abs(log(abs(R1v2./R2v2))).^0.5;

    widthDiffall=cat(2,widthDiffall,widthPP2-sd1);
end

meanDiff=mean(widthDiffall,2);

%% Plot

close all

figure('Position',[200 500 1200 800],'DefaultAxesFontSize',12,'renderer','painters');
t=tiledlayout(2,3,'TileIndexing', 'columnmajor','TileSpacing','compact','Padding','tight');

tiles=[1,2;3,4;5,6];

for ii=1:length(sd1)
    s1=nexttile(tiles(ii,1));
    hold on
    l1=plot(x(ii,:),y1(ii,:),'-b','LineWidth',1.5);
    l2=plot(x(ii,:),y2(ii,:),'-g','LineWidth',1);
    xlim([x(ii,1),x(ii,end)])
    ylim([-0.3,1.4])
    grid on
    box on

    legend([l1,l2],{['Width spectral: ',num2str(width1(ii))];['Width pulse pair: ',num2str(widthPP2(ii))]})
    title(['Width in: ',num2str(sd1(ii))]);

    s2=nexttile(tiles(ii,2));
    edgesIn=-1:0.01:1;
    hc=histcounts(widthDiffall(ii,:),edgesIn);
    bar(edgesIn(1:end-1)+0.05,hc,1);

    title(['Noisy - smooth: ',num2str(meanDiff(ii))])
end

print([figdir,'noiseTest/widthNoise_',num2str(noiseVar),'.png'],'-dpng','-r0');