%Compare attenuation calculated from sounding data vs ecmwf fake sounding
%from reanalysis data

clear all;
close all;

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/utils/');

project='aristo';

sondedir=['/scr/sci/romatsch/data/dropSonde/',project,'/'];
modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5/',project,'/'];

figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/soundECMWF/';

soundings=dir([sondedir,'*.eol']);

outMat=nan(size(soundings,1),9);

frq = 9.4406e10;

for ii=1:size(soundings,1)
    sonName=soundings(ii).name;
    disp(sonName);
    sonTime=datetime(str2num(sonName(2:5)),str2num(sonName(6:7)),str2num(sonName(8:9)),...
        str2num(sonName(11:12)),str2num(sonName(13:14)),str2num(sonName(15:16)));
    
    % sounding
    [outMat(ii,1),outMat(ii,3),outMat(ii,5),outMat(ii,7)]=f_atten_layers_sfcWnd([sondedir,soundings(ii).name],frq/1e+9);
    
    % get lon lat
    soundDataIn=importdata([sondedir,soundings(ii).name],' ',14);
    locLine=soundDataIn.textdata{5};
    quotesLoc=strfind(locLine,'''');
    commaLoc=strfind(locLine,',');
    
    if ~isempty(quotesLoc)
        lon=str2num(locLine((quotesLoc(1)+3):(commaLoc(end-1)-1)));
        lat=str2num(locLine((quotesLoc(2)+3):(commaLoc(end)-1)));
        alt=str2num(locLine((commaLoc(end)+1):end));
        
        %model
        [outMat(ii,2),outMat(ii,4),outMat(ii,6),outMat(ii,8),outMat(ii,9)]=f_atten_layers_sfcWnd_ecmwf_reanalysis(modeldir,...
            frq/1e+9,lon,lat,alt,sonTime);
    end
end

%% Plot differences

noAtt=any(isnan(outMat(:,1:4)),2);
outMat(noAtt,:)=[];

f1=figure('DefaultAxesFontSize',14);
set(f1,'Position',[200 500 1500 1200]);

subplot(5,1,1:2)
hold on;
l1=plot(outMat(:,1),'ob','markerfacecolor','b');
l2=plot(outMat(:,2),'oc','markerfacecolor','c');
l3=plot(outMat(:,3),'or','markerfacecolor','r');
l4=plot(outMat(:,4),'om','markerfacecolor','m');
legend([l1 l2 l3 l4],{'Liebe sound','Liebe model','ITU sound','ITU model'},'location','best');
title(['One way attenuation, ',project],'interpreter','none');

xlimits=[0 size(outMat,1)+1];
xlim(xlimits);
ylim([0 3]);

ylabel('Attenuation [dB]');
grid on

subplot(5,1,3:4)
hold on;
l1=plot(outMat(:,2)-outMat(:,1),'ob','markerfacecolor','b');
l2=plot(outMat(:,4)-outMat(:,3),'or','markerfacecolor','r');

xlim(xlimits);
ylimits=[-0.8 0.8];
ylim(ylimits);
plot(xlimits,[0 0],'-k');

resids=[abs(outMat(:,2)-outMat(:,1));abs(outMat(:,4)-outMat(:,3))];
text(2,ylimits(2)-0.1,['Mean offset: ',num2str(nanmean([outMat(:,2)-outMat(:,1);outMat(:,4)-outMat(:,3)]))],'fontsize',14);
text(2,ylimits(2)-0.2,['Mean diff: ',num2str(nanmean(resids))],'fontsize',14);
text(2,ylimits(2)-0.3,['Std: ',num2str(nanstd(resids))],'fontsize',14);


ylabel('Attenuation difference [dB]');

title(['Attenuation difference, model-sounding, ',project],'interpreter','none');

legend([l1 l2],{'Liebe','ITU'},'location','best');
grid on

subplot(5,1,5)
hold on;
l1=plot(outMat(:,9),'ok','markerfacecolor','k');

xlim(xlimits);
ylimits2=[0 25];
ylim(ylimits2);

ylabel('Distance [km]');

title(['Distance between model and sounding, ',project],'interpreter','none');
grid on

set(gcf,'PaperPositionMode','auto');
print(f1, [figdir,'att_modelMinusSounding_',project],'-dpng','-r0');

%% Scatter plot
f1=figure('DefaultAxesFontSize',14);
set(f1,'Position',[200 500 1500 600]);

subplot(1,2,1)
hold on
scatter(outMat(:,1),outMat(:,2),'b','filled');
scatter(outMat(:,3),outMat(:,4),'r','filled');

p = polyfit([outMat(:,1);outMat(:,3)],[outMat(:,2);outMat(:,4)],1);
x1 = 0:0.1:4;
y1 = polyval(p,x1);

xlim([0 3]);
ylim([0 3]);

plot(x1,y1,'-k','linewidth',1.5);

xlabel('One way attenuation (sounding) [dB]');
ylabel('One way attenuation (model) [dB]');
text(0.2,2.8,['y=',num2str(p(1)),'+',num2str(p(2)),'x'],'fontsize',14);

title(project)

subplot(1,2,2)
hold on
scatter(outMat(:,1),abs(outMat(:,2)-outMat(:,1)),'b','filled');
scatter(outMat(:,3),abs(outMat(:,4)-outMat(:,3)),'r','filled');

xlabel('One way attenuation (sounding) [dB]');
ylabel('Attenuation difference [dB]');

title(project)

set(gcf,'PaperPositionMode','auto');
print(f1, [figdir,'att_difference_scatt_',project],'-dpng','-r0');