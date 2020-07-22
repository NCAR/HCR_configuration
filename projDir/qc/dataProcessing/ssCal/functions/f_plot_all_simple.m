function f_plot_all_simple(PLTall,sigXlims,sigYlims,reflXlims,reflYlims,projName,directories,uniqueCases)
%% Plot all sig0 measured vs elev angle
close all

f1=figure;
set(gcf,'Position',[200 500 800 600]);
hold on

for ii=1:size(PLTall,1)
    if min(isnan(PLTall{ii,1}.sig0measured))==0
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
smRefl=[];
smSig0measured=[];

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isnan(PLT.sig0measured)
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
    smRefl=cat(2,smRefl,binData(:,2));
    smSig0measured=cat(2,smSig0measured,binData(:,3));
end


%% Plot smoothed reflectivity

legendGet=[];
colMap=cat(1,hsv(round(size(PLTall,1)/1.5)),pink(round(size(PLTall,1)/2)));

figure;
set(gcf,'Position',[200 500 800 600]);
hold on

for ii=1:size(PLTall,1)
    legendGet{end+1}=datestr(datetime(uniqueCases(ii,1:6)),'yyyymmdd_HHMMSS');
    plot(smElev(:,ii),smRefl(:,ii),'Color',colMap(ii,:),'linewidth',1.5);
end

set(gca,'XLim',reflXlims+[0 10]);
set(gca,'YLim',reflYlims);
xlabel('Elevation angle [deg]');
ylabel('Reflectivity [dBZ]');
title([{projName};{'Reflectivity vs elevation angle, smoothed'}],'interpreter','none');
legend(legendGet,'location','northeast','Interpreter','none');
grid on
set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_refl_all_smooth'],'-dpng','-r0')


%% Plot smoothed sig0measured

legendGet2={};

figure;
set(gcf,'Position',[200 500 800 600]);
hold on

for ii=1:size(PLTall,1)
    if min(isnan(smSig0measured(:,ii)))==0
        legendGet2{end+1}=datestr(datetime(uniqueCases(ii,1:6)),'yyyymmdd_HHMMSS');
        plot(smElev(:,ii),smSig0measured(:,ii),'Color',colMap(ii,:),'linewidth',1.5);
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
print([directories.figdir projName '_sigma0measured_totMean'],'-dpng','-r0')

%% Make poly plots
figure;
set(gcf,'Position',[200 500 800 600]);
hold on
grid on

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isfield(PLT,'polySig0Meas')
        l1=plot(PLT.polyElev,PLT.polySig0Meas,'-r','linewidth',1.5);
    end
    if isfield(PLT,'polySig0Model')
        l2=plot(PLT.polyElev,PLT.polySig0Model(1,:),'-b','linewidth',1.5);
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
legend([l1,l2],{'Sigma0 measured';'Sigma0 model FreilichVanhoff'});

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0poly_freiVan'],'-dpng','-r0')

figure;
set(gcf,'Position',[200 500 800 600]);
hold on
grid on

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isfield(PLT,'polySig0Meas')
        l1=plot(PLT.polyElev,PLT.polySig0Meas,'-r','linewidth',1.5);
    end
    if isfield(PLT,'polySig0Model')
        l2=plot(PLT.polyElev,PLT.polySig0Model(2,:),'-b','linewidth',1.5);
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
legend([l1,l2],{'Sigma0 measured';'Sigma0 model Wu'});

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0poly_wu'],'-dpng','-r0')

figure;
set(gcf,'Position',[200 500 800 600]);
hold on
grid on

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};
    if isfield(PLT,'polySig0Meas')
        l1=plot(PLT.polyElev,PLT.polySig0Meas,'-r','linewidth',1.5);
    end
    if isfield(PLT,'polySig0Model')
        l2=plot(PLT.polyElev,PLT.polySig0Model(3,:),'-b','linewidth',1.5);
    end
end

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([{projName};{'Poly fit to radar cross section vs elevation angle'}],'interpreter','none');
legend([l1,l2],{'Sigma0 measured';'Sigma0 model CoxMunk'});

set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0poly_coxMunk'],'-dpng','-r0')

%% Make residual plot
figure;
set(gcf,'Position',[200 500 800 600]);

s1=subplot(2,1,1);
hold on
grid on

resids=[];

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model')
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(1,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(ii,:),'linewidth',1.5);
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

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model')
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(2,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(ii,:),'linewidth',1.5);
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

for ii=1:size(PLTall,1)
    PLT=PLTall{ii};   
    if isfield(PLT,'polySig0Meas') & isfield(PLT,'polySig0Model')
        resids=[resids;PLT.polySig0Meas-PLT.polySig0Model(3,:)];
        l1=plot(PLT.polyElev,resids(end,:),'Color',colMap(ii,:),'linewidth',1.5);
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

end

