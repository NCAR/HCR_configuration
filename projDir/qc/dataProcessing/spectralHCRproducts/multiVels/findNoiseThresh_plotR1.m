function [noiseThresh,meanNoise,R2]=findNoiseThresh_plotR1R2(powIn,avNum,aa)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
close all

noiseThresh=nan;
meanNoise=nan;
R2=nan;

checkInds=0:5:770;

if ~ismember(aa,checkInds)
    return
end

noiseThresh=10.^(max(powIn)./10);
powLin=10.^(powIn./10);
powLin(isnan(powLin))=[];

stepDown=0.0000005;

powR2All=[];

while noiseThresh>min(powLin,[],2,'omitmissing');

    noiseThresh=noiseThresh-stepDown;

    powLin(powLin>noiseThresh)=[];
    sampleNum=length(powLin);

    % R1
    % xCalc=1:sampleNum;
    % fn=xCalc./sampleNum;
    % sig2r=(sum(fn.^2.*powLin,'omitmissing')./sum(powLin,'omitmissing')) ...
    %     -(sum(fn.*powLin,'omitmissing')./sum(powLin,'omitmissing')).^2;
    % sigN2=sampleNum.^2/12;
    % R1=sigN2/sig2r;

    % R2
    meanNoise=sum(powLin)/sampleNum;
    Q=sum(powLin.^2/sampleNum)-meanNoise.^2;
    R2=meanNoise.^2/(Q*avNum);

    powR2All=cat(1,powR2All,[noiseThresh,R2]);
end

powR2All(any(isnan(powR2All),2),:)=[];

if isempty(powR2All)
    return
end

powR2All(:,1)=real(10*log10(powR2All(:,1)));
outInd=find(powR2All(:,2)==1);
if ~isempty(outInd)
    noiseThreshOut=powR2All(outInd,2);
else
    below1=max(find(powR2All(:,2)<1));
    above1=min(find(powR2All(:,2)>1));
    if ~isempty(below1) & ~isempty(above1)
        noiseThreshOut=interp1([powR2All(below1,2),powR2All(above1,2)],[powR2All(below1,1),powR2All(above1,1)],1);
    else
        noiseThreshOut=nan;
    end
end

powR2All(:,1)=max(powIn)-powR2All(:,1);

figure('Position',[200 500 700 1000],'DefaultAxesFontSize',12);

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
plot([0,100],[1,1],'-k');
scatter(powR2All(:,1),powR2All(:,2),'filled','o','MarkerFaceColor','b');
scatter(max(powIn)-noiseThreshOut,1,'filled','o','MarkerFaceColor','r');
xlim([0,max(powR2All(:,1))+5]);
ylim([0,3]);

xlabel('dB below peak');
ylabel('R2');

grid on
box on

s2=nexttile(2);
hold on
plot([0,length(powIn)],[noiseThreshOut,noiseThreshOut],'-k','LineWidth',2);
plot(powIn,'-b','LineWidth',2);

xlim([0,length(powIn)]);

%xlabel('dB below peak');
ylabel('Power (dB)');

grid on
box on

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/airMotion/cases/width/noiseFloor/';
set(gcf,'PaperPositionMode','auto');
print([figdir,'noiseFloor_',num2str(aa),'_smoothing',num2str(avNum)],'-dpng','-r0');
end