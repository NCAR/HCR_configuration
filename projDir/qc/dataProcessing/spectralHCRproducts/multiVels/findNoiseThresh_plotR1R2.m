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

powR1R2All=[];

while noiseThresh>min(powLin,[],2,'omitmissing');

    noiseThresh=noiseThresh-stepDown;

    powLin(powLin>noiseThresh)=[];
    sampleNum=length(powLin);

    % R1
    xCalc=1:sampleNum;
    fn=xCalc./sampleNum;
    sig2r=(sum(fn.^2.*powLin,'omitmissing')./sum(powLin,'omitmissing')) ...
        -(sum(fn.*powLin,'omitmissing')./sum(powLin,'omitmissing')).^2;
    % sigN2=sampleNum.^2/12;
    sigN2=1/12;
    R1=sigN2/sig2r;

    % R2
    meanNoise=sum(powLin)/sampleNum;
    Q=sum(powLin.^2/sampleNum)-meanNoise.^2;
    R2=meanNoise.^2/(Q*avNum);

    powR1R2All=cat(1,powR1R2All,[noiseThresh,R1,R2]);
end

powR1R2All(any(isnan(powR1R2All),2),:)=[];

if isempty(powR1R2All)
    return
end

% Find intersection point
powR1R2All(:,1)=real(10*log10(powR1R2All(:,1)));
Rdiff=powR1R2All(:,2)-powR1R2All(:,3);
outInd=find(Rdiff==0);
y=nan;
if ~isempty(outInd)
    noiseThreshOut=powR1R2All(outInd,2);
else
    beforeX=max(find(Rdiff>0));
    afterX=beforeX+1;
    if ~isempty(beforeX) & ~isempty(afterX) & afterX<=size(powR1R2All,1)
        A=[powR1R2All(beforeX,1),powR1R2All(beforeX,2)];
        B=[powR1R2All(afterX,1),powR1R2All(afterX,2)];
        C=[powR1R2All(beforeX,1),powR1R2All(beforeX,3)];
        D=[powR1R2All(afterX,1),powR1R2All(afterX,3)];
        a1=B(2)-A(2);
        b1=A(1)-B(1);
        c1=a1*(A(1))+b1*(A(2));
        a2=D(2)-C(2);
        b2=C(1)-D(1);
        c2=a2*C(1)+b2*C(2);
        det=a1*b2-a2*b1;
        noiseThreshOut=(b2*c1-b1*c2)/det;
        y=(a1*c2-a2*c1)/det;
        %noiseThreshOut=interp1([powR1R2All(beforeX,2),powR1R2All(afterX,2)],[powR1R2All(beforeX,1),powR1R2All(afterX,1)],1);
    else
        noiseThreshOut=nan;
    end
end

powR1R2All(:,1)=max(powIn)-powR1R2All(:,1);

figure('Position',[200 500 700 1000],'DefaultAxesFontSize',12);

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
scatter(powR1R2All(:,1),powR1R2All(:,2),'filled','o','MarkerFaceColor','b');
scatter(powR1R2All(:,1),powR1R2All(:,3),'filled','o','MarkerFaceColor','g');
scatter(max(powIn)-noiseThreshOut,y,'filled','o','MarkerFaceColor','r');
xlim([0,max(powR1R2All(:,1))+5]);
ylim([0,3]);

legend({'R1';'R2';['Intersection at ',num2str(y)]},'Location','northwest');

xlabel('dB below peak');
ylabel('R1 and R2');

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