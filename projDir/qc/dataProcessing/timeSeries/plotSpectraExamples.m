function plotSpectraExamples(data,momentsSpec,specPowerDB,sampleNum,plotRangeKM,plotInd,startInd,ii,ylimUpper,figdir,project)

% X-axis
specVelVec=-pi:2*pi/(sampleNum):pi;
specVelVec=specVelVec(1:end-1);
vel=data.lambda./(4.*pi.*data.prtThis).*specVelVec;

% Find index of specified range
rangeInd=min(find((data.range./1000)>=plotRangeKM(plotInd)));

% Create figure
f0=figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12);

% Plot spectra at specified range
s1=subplot(3,1,1);
plot(vel,specPowerDB.V(rangeInd,:),'-b','LineWidth',2);
xlabel('Velocity (m s^{-1})');
ylabel('Power (dB)');
xlim([vel(1) vel(end)]);
title({[datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS')];['Power at ',num2str(plotRangeKM(plotInd)),' km range.']});

text(vel(end)+0.3,-60,{['WIDTH: ',num2str(momentsSpec.width(rangeInd,ii),2),' ms^{-1}']; ...
    ['SKEW: ',num2str(momentsSpec.skew(rangeInd,ii),2)]; ...
    ['KURT: ',num2str(momentsSpec.kurt(rangeInd,ii),2)]});

% Waterfall plot
s2=subplot(3,1,2:3);
colormap('jet');
hold on
surf(vel,data.range./1000,specPowerDB.V,'EdgeColor','none');
view(2)
caxis([-60 -30]);
colorbar
xlabel('Velocity (m s^{-1})');
ylabel('Range (km)');
xlim([vel(1) vel(end)])
ylim([data.range(1)./1000,ylimUpper])
title('Power (dB)');

% Plot line at specified range
plot([vel(1),vel(end)],[plotRangeKM(plotInd),plotRangeKM(plotInd)],'-w');

% Align x axis of subplots
drawnow;
s1Pos=s1.Position;
s2Pos=s2.Position;
s1.Position=[s1Pos(1),s1Pos(2),s2Pos(3),s1Pos(4)];

set(gcf,'PaperPositionMode','auto')
print(f0,[figdir,project,'_spectra_',datestr(data.time(startInd),'yyyymmdd_HHMMSS'),'_',num2str(plotRangeKM(plotInd)),'km.png'],'-dpng','-r0');

end