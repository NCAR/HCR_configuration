function plot_all(PLTall,sigXlims,sigYlims,reflXlims,reflYlims,projName,directories,uniqueCases)
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
title([projName ' measured radar cross section vs elevation angle']);
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
        binData(jj,:)=nanmean(data,1);
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
title([projName ' reflectivity vs elevation angle, smoothed']);
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
title([projName ' measured radar cross section vs elevation angle, smoothed']);
legend(legendGet2,'location','northeast','Interpreter','none');
grid on
set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0measured_all_smooth'],'-dpng','-r0')

%% Make sig0 measured plot with mean and std

figure;
set(gcf,'Position',[200 500 800 600]);
hold on
grid on

errorbar(nanmean(smElev,2),nanmean(smSig0measured,2),nanstd(smSig0measured,1,2),'Color','k', 'LineStyle', 'none');
plot(nanmean(smElev,2),nanmean(smSig0measured,2),'Color','r','linewidth',2);

set(gca,'XLim',sigXlims);
set(gca,'YLim',sigYlims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title([projName ' mean measured radar cross section vs elevation angle']);
set(gcf,'PaperPositionMode','auto')
print([directories.figdir projName '_sigma0measured_totMean'],'-dpng','-r0')
end

