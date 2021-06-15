function plotMresult(data,m,result,outname)
%Plot m and result for PID
close all

m(m==0)=nan;
result(isnan(data.DBZ))=nan;

timeMat=repmat(data.time,size(data.TEMP,1),1);

f1=figure('DefaultAxesFontSize',12,'Position',[0 300 2300 1200]);

titles={'DBZ','LDR','VEL','WIDTH','TEMP'};
ylimits=[0 3];

colormap jet

for ii=1:5
    subplot(3,2,ii);
    surf(data.time,data.asl./1000,squeeze(m(ii,:,:)),'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([0 1]);
    colorbar;
    ylabel('Altitude (km)');
    title(titles{ii});
    grid on
end

subplot(3,2,6);
surf(data.time,data.asl./1000,squeeze(result),'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0 100]);
colorbar;
ylabel('Altitude (km)');
title([outname,' weighted']);

set(gcf,'PaperPositionMode','auto')
print(f1,['/scr/sci/romatsch/HCR/pid/hcrOnly/socrates/socrates_pid_',...
datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_',outname],'-dpng','-r0')

end

