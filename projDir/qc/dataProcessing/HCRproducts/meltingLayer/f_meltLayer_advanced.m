% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function [BBfinishedOut iceLevAsl offset]= f_meltLayer_advanced(data,zeroAdjustMeters,thresholds,figdir)

debugFig=0;

%% Find zero degree altitude altitudes and indices

disp('Searching zero degree altitude ...');

oneGate=data.range(2)-data.range(1);

[layerAlts,layerInds]=zeroDegIso(data);

%% Truncate to non missing and regions with sub deg temps
gapSecs=10;
tempMinMax=[-1,7];
nonMissingInds=findNonMissingInds(data,gapSecs);

% Temperature too high
minTemp=min(data.TEMP,[],1,'omitnan');
nonMissingInds(minTemp>tempMinMax(2))=0;

dataInVars=fields(data);

dataShort=[];
for ii=1:length(dataInVars)
    dataShort.(dataInVars{ii})=data.(dataInVars{ii})(:,nonMissingInds==1);
end

%% Prepare VEL data

% Reverse sign in zenith and remove upward motion
dataShort.VEL_MASKED(:,dataShort.elevation>0)=-dataShort.VEL_MASKED(:,dataShort.elevation>0);
dataShort.VEL_MASKED(dataShort.VEL_MASKED<0)=nan;

% Median filter
dataShort.VEL_MASKED=medfilt2(dataShort.VEL_MASKED,[3,7]);

%% Velocity derivative
velDiff=diff(dataShort.VEL_MASKED,1);

%% Input fields
velDiff=cat(1,nan(size(dataShort.time)),velDiff);
velField=velDiff;
velField(dataShort.TEMP<tempMinMax(1) | velField<0.2)=nan;
ldrField=dataShort.LDR;
ldrField(dataShort.TEMP<tempMinMax(1) | ldrField<-18)=nan;

%% Mask
maskFields=(~isnan(velField) | ~isnan(ldrField));
maskFields=bwareaopen(maskFields,5);

velField(maskFields==0)=nan;
ldrField(maskFields==0)=nan;

%% Maxima
[velMax,velMaxInd]=max(velField,[],1,'omitnan');
%velMaxInd(isnan(velMax))=nan;

[ldrMax,ldrMaxInd]=max(ldrField,[],1,'omitnan');
%ldrMaxInd(isnan(ldrMax))=nan;

%% For plot
velMaxPlot=velMax(~isnan(velMax));
velMaxPlotInd=velMaxInd(~isnan(velMax));
colvel=1:length(dataShort.time);
colvel(isnan(velMax))=[];
linvel=sub2ind(size(dataShort.VEL_MASKED),velMaxPlotInd,colvel);

ldrMaxPlot=ldrMax(~isnan(ldrMax));
ldrMaxPlotInd=ldrMaxInd(~isnan(ldrMax));
colldr=1:length(dataShort.time);
colldr(isnan(ldrMax))=[];
linldr=sub2ind(size(dataShort.LDR),ldrMaxPlotInd,colldr);

maskPlot=double(maskFields);
maskPlot(maskPlot==0)=nan;

%% Plot 1

disp('Plotting ...')
ylimits=[0,5];

close all

newInds=1:round(length(dataShort.time)/2000):length(dataShort.time);

% Resample for plotting
newDBZ=dataShort.DBZ_MASKED(:,newInds);
newLDR=dataShort.LDR(:,newInds);
newVEL=dataShort.VEL_MASKED(:,newInds);
newVELdiff=velDiff(:,newInds);
newASL=dataShort.asl(:,newInds);
newTime=dataShort.time(newInds);

newMask=maskPlot(:,newInds);

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible','on');

ax1=subplot(4,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
view(2);
colMapDBZ(sub1);
scatter(dataShort.time(~isnan(velMax)),dataShort.asl(linvel)./1000,3.5,'filled','MarkerFaceColor','m');
scatter(dataShort.time(~isnan(ldrMax)),dataShort.asl(linldr)./1000,3.5,'filled','MarkerFaceColor','g');
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Reflectivity and melting layer')
grid on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.87 0.21];
ax1.SortMethod='childorder';

% LDR

ax2=subplot(4,1,2);
hold on;
surf(newTime,newASL./1000,newLDR,'edgecolor','none');
view(2);
caxis([-25 -5]);
ax2.Colormap=jet;
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('LDR')
grid on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.525 0.87 0.21];

% VEL

ax3=subplot(4,1,3);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ax3.Colormap=jet;
caxis([0 6]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('VEL')
grid on
ax3.Position=[0.06 0.287 0.87 0.21];

linkaxes([ax1 ax2 ax2 ax3],'xy');

% Diff VEL
ax4=subplot(4,1,4);
hold on;
surf(newTime,newASL./1000,newVELdiff,'edgecolor','none');
view(2);
ax4.Colormap=jet;
caxis([0 0.6]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('VEL')
grid on
ax4.Position=[0.06 0.05 0.87 0.21];

linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer1',datestr(dataShort.time(1),formatOut),'_to_',datestr(dataShort.time(end),formatOut)],'-dpng','-r0');

%% Plot 2

fig2=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible','on');

ax1=subplot(4,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newMask,'edgecolor','none');
view(2);
%scatter(dataShort.time(~isnan(velMax)),dataShort.asl(linvel)./1000,3.5,'filled','MarkerFaceColor','m');
%scatter(dataShort.time(~isnan(ldrMax)),dataShort.asl(linldr)./1000,3.5,'filled','MarkerFaceColor','g');
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
%title('Reflectivity and melting layer')
grid on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.87 0.21];
ax1.SortMethod='childorder';

% % LDR
% 
% ax2=subplot(4,1,2);
% hold on;
% surf(newTime,newASL./1000,newLDR,'edgecolor','none');
% view(2);
% caxis([-25 -5]);
% ax2.Colormap=jet;
% colorbar
% ylim(ylimits);
% ylabel('Altitude (km)');
% xlim([dataShort.time(1),dataShort.time(end)]);
% title('LDR')
% grid on
% set(gca,'xticklabel',[])
% ax2.Position=[0.06 0.525 0.87 0.21];
% 
% % VEL
% 
% ax3=subplot(4,1,3);
% hold on;
% surf(newTime,newASL./1000,newVEL,'edgecolor','none');
% view(2);
% ax3.Colormap=jet;
% caxis([0 6]);
% colorbar
% ylim(ylimits);
% ylabel('Altitude (km)');
% xlim([dataShort.time(1),dataShort.time(end)]);
% title('VEL')
% grid on
% ax3.Position=[0.06 0.287 0.87 0.21];
% 
% linkaxes([ax1 ax2 ax2 ax3],'xy');
% 
% % Diff VEL
% ax4=subplot(4,1,4);
% hold on;
% surf(newTime,newASL./1000,newVELdiff,'edgecolor','none');
% view(2);
% ax4.Colormap=jet;
% caxis([0 0.6]);
% colorbar
% ylim(ylimits);
% ylabel('Altitude (km)');
% xlim([dataShort.time(1),dataShort.time(end)]);
% title('VEL')
% grid on
% ax4.Position=[0.06 0.05 0.87 0.21];
% 
% linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer2',datestr(dataShort.time(1),formatOut),'_to_',datestr(dataShort.time(end),formatOut)],'-dpng','-r0');

end