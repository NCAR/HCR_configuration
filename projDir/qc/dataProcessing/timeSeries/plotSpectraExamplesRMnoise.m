function plotSpectraExamplesRMnoise(data,momentsSpec,momentsSpecRMnoise,momentsSpecRMnoiseSmooth, ...
    specVelAdj,specPowerDBadj,powerDBsmooth,powerRMnoiseDB,powerRMnoiseDBsmooth, ...
    locsMax,locsMin,plotRangeKM,plotInd,startInd,ii,ylimUpper,figdir,project)

% Find index of specified range
rangeInd=min(find((data.range./1000)>=plotRangeKM(plotInd)));

velX=specVelAdj(rangeInd,:);

% Create figure
f0=figure('Position',[200 500 2200 1200],'DefaultAxesFontSize',12);

% Plot spectra at specified range
s1=subplot(3,3,1);
plot(velX,specPowerDBadj(rangeInd,:),'-b','LineWidth',2);
xlabel('Velocity (m s^{-1})');
ylabel('Power (dB)');
xlim([velX(1) velX(end)]);
title({[datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS')];['Power at ',num2str(plotRangeKM(plotInd)),' km range.']});

text(velX(end)+0.3,-60,{['WIDTH: ',num2str(momentsSpec.width(rangeInd,ii),2),' ms^{-1}']; ...
    ['SKEW: ',num2str(momentsSpec.skew(rangeInd,ii),2)]; ...
    ['KURT: ',num2str(momentsSpec.kurt(rangeInd,ii),2)]});

grid on
box on

% Spectra without noise smooth
s2=subplot(3,3,2);
hold on
plot(velX,specPowerDBadj(rangeInd,:),'-b','LineWidth',1);
plot(velX,powerDBsmooth(rangeInd,:),'-c','LineWidth',2);
plot(velX,powerRMnoiseDBsmooth(rangeInd,:),'-r','LineWidth',2);
scatter(velX(locsMax.(['max',num2str(ii)])),powerDBsmooth(rangeInd,locsMax.(['max',num2str(ii)])),'filled','MarkerFaceColor','g');
scatter(velX(locsMin.(['min',num2str(ii)])),powerDBsmooth(rangeInd,locsMin.(['min',num2str(ii)])),'filled','MarkerFaceColor','m');
xlabel('Velocity (m s^{-1})');
ylabel('Power (dB)');
xlim([velX(1) velX(end)]);
title({[datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS')];['Power at ',num2str(plotRangeKM(plotInd)),' km range.']});

text(velX(end)+0.3,-60,{['WIDTH: ',num2str(momentsSpecRMnoiseSmooth.width(rangeInd,ii),2),' ms^{-1}']; ...
    ['SKEW: ',num2str(momentsSpecRMnoiseSmooth.skew(rangeInd,ii),2)]; ...
    ['KURT: ',num2str(momentsSpecRMnoiseSmooth.kurt(rangeInd,ii),2)]});

s2.SortMethod='childorder';

grid on
box on

% Spectra without noise
s3=subplot(3,3,3);
hold on
plot(velX,specPowerDBadj(rangeInd,:),'-b','LineWidth',1);
plot(velX,powerRMnoiseDB(rangeInd,:),'-r','LineWidth',2);
xlabel('Velocity (m s^{-1})');
ylabel('Power (dB)');
xlim([velX(1) velX(end)]);
title({[datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS')];['Power at ',num2str(plotRangeKM(plotInd)),' km range.']});

text(velX(end)+0.3,-60,{['WIDTH: ',num2str(momentsSpecRMnoise.width(rangeInd,ii),2),' ms^{-1}']; ...
    ['SKEW: ',num2str(momentsSpecRMnoise.skew(rangeInd,ii),2)]; ...
    ['KURT: ',num2str(momentsSpecRMnoise.kurt(rangeInd,ii),2)]});

s3.SortMethod='childorder';

grid on
box on

% Waterfall plot
s4=subplot(3,3,[4,7]);
colormap('jet');
hold on
surf(1:size(specPowerDBadj,2),data.range./1000,specPowerDBadj,'EdgeColor','none');
view(2)
caxis([-60 -30]);
colorbar
%xlabel('Velocity (m s^{-1})');
ylabel('Range (km)');
xlim([1,size(specPowerDBadj,2)])
ylim([data.range(1)./1000,ylimUpper])
title('Power (dB)');

% Plot line at specified range
plot([1,size(specPowerDBadj,2)],[plotRangeKM(plotInd),plotRangeKM(plotInd)],'-w');

% Align x axis of subplots
drawnow;
s1Pos=s1.Position;
s4Pos=s4.Position;
s1.Position=[s1Pos(1),s1Pos(2),s4Pos(3),s1Pos(4)];

grid on
box on

% Waterfall plot censored smooth
s5=subplot(3,3,[5,8]);
colormap('jet');
hold on
surf(1:size(specPowerDBadj,2),data.range./1000,powerRMnoiseDBsmooth,'EdgeColor','none');
view(2)
caxis([-60 -30]);
colorbar
%xlabel('Velocity (m s^{-1})');
ylabel('Range (km)');
xlim([1,size(specPowerDBadj,2)])
ylim([data.range(1)./1000,ylimUpper])
title('Power (dB)');

% Plot line at specified range
plot([1,size(specPowerDBadj,2)],[plotRangeKM(plotInd),plotRangeKM(plotInd)],'-k');

% Align x axis of subplots
drawnow;
s2Pos=s2.Position;
s5Pos=s5.Position;
s2.Position=[s2Pos(1),s2Pos(2),s5Pos(3),s2Pos(4)];

grid on
box on

% Waterfall plot censored
s6=subplot(3,3,[6,9]);
colormap('jet');
hold on
surf(1:size(specPowerDBadj,2),data.range./1000,powerRMnoiseDB,'EdgeColor','none');
view(2)
caxis([-60 -30]);
colorbar
%xlabel('Velocity (m s^{-1})');
ylabel('Range (km)');
xlim([1,size(specPowerDBadj,2)])
ylim([data.range(1)./1000,ylimUpper])
title('Power (dB)');

% Plot line at specified range
plot([1,size(specPowerDBadj,2)],[plotRangeKM(plotInd),plotRangeKM(plotInd)],'-k');

% Align x axis of subplots
drawnow;
s3Pos=s3.Position;
s6Pos=s6.Position;
s3.Position=[s3Pos(1),s3Pos(2),s6Pos(3),s3Pos(4)];

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f0,[figdir,project,'_spectra_',datestr(data.time(startInd),'yyyymmdd_HHMMSS'),'_',num2str(plotRangeKM(plotInd)),'km.png'],'-dpng','-r0');

end