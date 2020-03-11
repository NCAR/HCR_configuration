function f_plot_all(PLTall,sigXlims,sigYlims,reflXlims,reflYlims,projName,directories,uniqueCases)
%% Plot all sig0 measured vs elev angle
close all

f1=figure;
set(gcf,'Position',[200 500 800 600]);
hold on

for ii=1:size(PLTall,1)
    if min(isnan(PLTall{ii,1}.sig0measured))==0 & uniqueCases(ii,14)==1
        sig0measuredsmaller=PLTall{ii,1}.sig0measured(find(PLTall{ii,1}.rota<=180));
        elevSmaller=PLTall{ii,1}.elev(find(PLTall{ii,1}.rota<=180));
        
        sig0measuredbigger=PLTall{ii,1}.sig0measured(find(PLTall{ii,1}.rota>180));
        elevBigger=PLTall{ii,1}.elev(find(PLTall{ii,1}.rota>180));
        
        plot(elevSmaller,sig0measuredsmaller,'r');
        plot(elevBigger,sig0measuredbigger,'b');
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Measured radar cross section vs elevation angle'}],'interpreter','none');
legend('Rot angle <= 180 deg','Rot angle > 180 deg','location','southwest');
plabel = [directories.figdir projName '_sig0measured_all'];
grid on
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')


%% To smooth, put data in bins

edges=(0:0.5:25);

smElev=[];
%smRefl=[];
smSig0measured=[];

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isnan(PLT.sig0measured) | uniqueCases(ii,14)==0
        PLT.sig0measured=nan(size(PLT.elev));
    end
    elRefSig=cat(2,PLT.elev,PLT.refl,PLT.sig0measured);
    binInds=discretize(PLT.elev,edges);
    
    binData=nan(length(edges-1),3);
    
    for jj=1:length(edges-1)
        rowInds=find(binInds==jj);
        data=elRefSig(rowInds,:);
        binData(jj,:)=mean(data,1,'omitnan');
    end
    
    smElev=cat(2,smElev,binData(:,1));
    %smRefl=cat(2,smRefl,binData(:,2));
    smSig0measured=cat(2,smSig0measured,binData(:,3));
end

%% Plot smoothed sig0measured

legendGet2={};
colMap=cat(1,hsv(round(sum(uniqueCases(:,14))/1.5)),pink(round(sum(uniqueCases(:,14)/2))));

figure;
set(gcf,'Position',[200 500 800 600]);
hold on

countCol=1;
for ii=1:size(PLTall,1)
    if min(isnan(smSig0measured(:,ii)))==0 & uniqueCases(ii,14)==1
        legendGet2{end+1}=datestr(datetime(uniqueCases(ii,1:6)),'yyyymmdd_HHMMSS');
        plot(smElev(:,ii),smSig0measured(:,ii),'Color',colMap(countCol,:),'linewidth',1.5);
        countCol=countCol+1;
    end
end

set(gca,'XLim',sigXlims+[0 10]);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Measured radar cross section vs elevation angle, smoothed'}],'interpreter','none');
legend(legendGet2,'location','northeast','Interpreter','none');
grid on
set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0measured_all_smooth'],'-dpng','-r0')

%% Make sig0 measured plot with mean and std

figure;
set(gcf,'Position',[200 500 800 600]);
hold on
grid on

errorbar(mean(smElev,2,'omitnan'),mean(smSig0measured,2,'omitnan'),std(smSig0measured,1,2,'omitnan'),'Color','k', 'LineStyle', 'none');
plot(mean(smElev,2,'omitnan'),mean(smSig0measured,2,'omitnan'),'Color','r','linewidth',2);

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Mean measured radar cross section vs elevation angle'}],'interpreter','none');
set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0measured_totMeanCases'],'-dpng','-r0')

%% Make residual plot
figure;
set(gcf,'Position',[200 500 800 600]);

s1=subplot(2,1,1);
hold on
grid on

resids=[];%% To smooth, put data in bins

edges=(0:0.5:25);

smElev=[];
%smRefl=[];
smSig0measured=[];

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isnan(PLT.sig0measured) | uniqueCases(ii,14)==0
        PLT.sig0measured=nan(size(PLT.elev));
    end
    elRefSig=cat(2,PLT.elev,PLT.refl,PLT.sig0measured);
    binInds=discretize(PLT.elev,edges);
    
    binData=nan(length(edges-1),3);
    
    for jj=1:length(edges-1)
        rowInds=find(binInds==jj);
        data=elRefSig(rowInds,:);
        binData(jj,:)=mean(data,1,'omitnan');
    end
    
    smElev=cat(2,smElev,binData(:,1));
    %smRefl=cat(2,smRefl,binData(:,2));
    smSig0measured=cat(2,smSig0measured,binData(:,3));
end
countCol=1;

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model') & uniqueCases(ii,14)==1   
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(1,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(countCol,:),'linewidth',1.5);
        countCol=countCol+1;
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
title([{[projName ' FreiVan']};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
yticks(-5:1:5)
legend(legendGet2,'location','eastoutside','Interpreter','none');
drawnow
s1pos=s1.InnerPosition;

s2=subplot(2,1,2);
hold on
grid on

if ~isempty(resids)
    meanRes=mean(resids,1,'omitnan');
    stdRes=std(resids,1,'omitnan');
    
    errorbar(PLT.polyElev,meanRes,stdRes,'Color','k', 'LineStyle', 'none');
    plot(PLT.polyElev,meanRes,'-r','linewidth',1.5);
    
    meanRes(PLT.polyElev<5 | PLT.polyElev>15)=nan;
    text(sigXlims(1)+2,-2,['Max bias (5-15 deg): ',num2str(max(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-3,['Min bias (5-15 deg): ',num2str(min(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-4,['Mean bias: ',num2str(mean(reshape(meanRes,1,[]),'omitnan')),' dB'],'fontsize',12);
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
yticks(-5:1:5)

s2pos=s2.InnerPosition;
s2.InnerPosition=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0residuals_freiVan'],'-dpng','-r0')

figure;
set(gcf,'Position',[200 500 800 600]);

s1=subplot(2,1,1);
hold on
grid on

resids=[];
countCol=1;
for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model') & uniqueCases(ii,14)==1
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(2,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(countCol,:),'linewidth',1.5);
        countCol=countCol+1;
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
title([{[projName ' Wu']};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
yticks(-5:1:5)
legend(legendGet2,'location','eastoutside','Interpreter','none');
drawnow
s1pos=s1.InnerPosition;

s2=subplot(2,1,2);
hold on
grid on

if ~isempty(resids)
    meanRes=mean(resids,1,'omitnan');
    stdRes=std(resids,1,'omitnan');
    
    errorbar(PLT.polyElev,meanRes,stdRes,'Color','k', 'LineStyle', 'none');
    plot(PLT.polyElev,meanRes,'-r','linewidth',1.5);
    
    meanRes(PLT.polyElev<5 | PLT.polyElev>15)=nan;
    text(sigXlims(1)+2,-2,['Max bias (5-15 deg): ',num2str(max(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-3,['Min bias (5-15 deg): ',num2str(min(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-4,['Mean bias: ',num2str(mean(reshape(meanRes,1,[]),'omitnan')),' dB'],'fontsize',12);
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
yticks(-5:1:5)

s2pos=s2.InnerPosition;
s2.InnerPosition=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0residuals_wu'],'-dpng','-r0')

figure;
set(gcf,'Position',[200 500 800 600]);

s1=subplot(2,1,1);
hold on
grid on

resids=[];
countCol=1;
for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model') & uniqueCases(ii,14)==1
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(3,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(countCol,:),'linewidth',1.5);
        countCol=countCol+1;
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
title([{[projName ' CoxMunk']};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
yticks(-5:1:5)
legend(legendGet2,'location','eastoutside','Interpreter','none');
drawnow
s1pos=s1.InnerPosition;

s2=subplot(2,1,2);
hold on
grid on

if ~isempty(resids)
    meanRes=mean(resids,1,'omitnan');
    stdRes=std(resids,1,'omitnan');
    
    errorbar(PLT.polyElev,meanRes,stdRes,'Color','k', 'LineStyle', 'none');
    plot(PLT.polyElev,meanRes,'-r','linewidth',1.5);
    
    meanRes(PLT.polyElev<5 | PLT.polyElev>15)=nan;
    text(sigXlims(1)+2,-2,['Max bias (5-15 deg): ',num2str(max(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-3,['Min bias (5-15 deg): ',num2str(min(meanRes)),' dB'],'fontsize',12);
    text(sigXlims(1)+2,-4,['Mean bias: ',num2str(mean(reshape(meanRes,1,[]),'omitnan')),' dB'],'fontsize',12);
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',[-5 5]);
xlabel('Elevation angle [deg]');
ylabel('Bias [dB]');
yticks(-5:1:5)

s2pos=s2.InnerPosition;
s2.InnerPosition=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0residuals_coxMunk'],'-dpng','-r0')

%% Combine all data

sig0measAll=[];
sig0modAll=[];
elevAll=[];

% Cat all data together
for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if min(isnan(smSig0measured(:,ii)))==0 & uniqueCases(ii,14)==1
        sig0meas1=PLT.sig0measured;
        sig0mod1=PLT.sig0model;
        elev1=PLT.elev;
        
        % Only use data between 5 and 15 deg elevation
        sig0meas1(elev1<5 | elev1>15)=nan;
        sig0mod1((elev1<5 | elev1>15),:)=nan;
        
        sig0measAll=cat(1,sig0measAll,sig0meas1);
        sig0modAll=cat(1,sig0modAll,sig0mod1);
        elevAll=cat(1,elevAll,elev1);
    end    
end

%% Meas-Model
measMinusMod=nan(size(elevAll,1),3);
measMinusMod(:,1)=sig0measAll-sig0modAll(:,2);
measMinusMod(:,2)=sig0measAll-sig0modAll(:,5);
measMinusMod(:,3)=sig0measAll-sig0modAll(:,8);

meanBias=mean(measMinusMod,1,'omitnan');
stdBias=std(measMinusMod,[],1,'omitnan');

%% Bin data
elMeasModDiff=cat(2,elevAll,sig0measAll,sig0modAll,measMinusMod);
binIndsAll=discretize(elevAll,edges);

meanDataAll=nan(length(edges-1),size(elMeasModDiff,2));
stdDataAll=nan(length(edges-1),size(elMeasModDiff,2));

for jj=1:length(edges-1)
    rowIndsAll=find(binIndsAll==jj);
    data=elMeasModDiff(rowIndsAll,:);
    meanDataAll(jj,:)=mean(data,1,'omitnan');
    stdDataAll(jj,:)=std(data,[],1,'omitnan');
end
    
%% Plot mean and std

figure('DefaultAxesFontSize',11);
set(gcf,'Position',[200 500 800 800]);
set(gcf,'renderer','painters');

subplot(2,1,1)
hold on

errorbar(edges+0.25,meanDataAll(:,2),stdDataAll(:,2),'Color','k', 'LineStyle', 'none');
plot(edges+0.25,meanDataAll(:,2),'Color','r','linewidth',2);

yticks(-4:2:14);
set(gca,'XLim',[5 15]);
set(gca,'YLim',[-5 15]);
grid on
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Mean measured radar cross section'}],'interpreter','none');
set(gcf,'PaperPositionMode','auto')

subplot(2,1,2)
hold on

errorbar(edges+0.25,meanDataAll(:,12),stdDataAll(:,12),'Color','c', 'LineStyle', 'none');
l1=plot(edges+0.25,meanDataAll(:,12),'Color','c','linewidth',2);
errorbar(edges+0.25,meanDataAll(:,13),stdDataAll(:,13),'color',[0 0.5 0], 'LineStyle', 'none');
l2=plot(edges+0.25,meanDataAll(:,13),'color',[0 0.5 0],'linewidth',2);
errorbar(edges+0.25,meanDataAll(:,14),stdDataAll(:,14),'Color','g', 'LineStyle', 'none');
l3=plot(edges+0.25,meanDataAll(:,14),'Color','g','linewidth',2);

yticks(-8:1:8);
set(gca,'XLim',[5 15]);
set(gca,'YLim',[-4 7]);
grid on
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title(['Measured minus modeled radar cross section'],'interpreter','none');
legend([l1,l2,l3],{'FV model','Wu model','CM model'},'location','northeast','fontsize',12);
text(5.5,6.5,['FV model: mean bias ',num2str(meanBias(1),2),' dB, st. dev. ',num2str(stdBias(1),2),' dB'],'fontsize',12);
text(5.5,5.5,['Wu model: mean bias ',num2str(meanBias(2),2),' dB, st. dev. ',num2str(stdBias(2),2),' dB'],'fontsize',12);
text(5.5,4.5,['CM model: mean bias ',num2str(meanBias(3),2),' dB, st. dev. ',num2str(stdBias(3),2),' dB'],'fontsize',12);

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_bias'],'-dpng','-r0')

end

