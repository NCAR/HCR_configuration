% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

infileCSET='/scr/snow2/rsfdata/projects/cset/GV/RF03H.20150709.163000_233800.PNI.nc';
csetDT=datetime(2015,7,9);
infileSOC='/scr/snow2/rsfdata/projects/socrates/GV/RF11H.20180217.013000_061700.PNI.nc';
socDT=datetime(2018,2,17);
infileOTREC='/scr/snow1/rsfdata/projects/otrec/GV/OTRECrf01.nc';
otrecDT=datetime(2019,7,7);

gvtimeIn=ncread(infileCSET,'Time');
gvtimeCSET=csetDT+seconds(gvtimeIn);
gvaltCSET=ncread(infileCSET,'GGALT');

if min(size(gvaltCSET))~=1
    gvaltCSET=gvaltCSET(1,:)';
end

gvtimeIn=ncread(infileSOC,'Time');
gvtimeSOC=socDT+seconds(gvtimeIn);
gvaltSOC=ncread(infileSOC,'GGALT');

if min(size(gvaltSOC))~=1
    gvaltSOC=gvaltSOC(1,:)';
end

gvtimeIn=ncread(infileOTREC,'Time');
gvtimeOTREC=otrecDT+seconds(gvtimeIn);
gvaltOTREC=ncread(infileOTREC,'GGALT');

if min(size(gvaltOTREC))~=1
    gvaltOTREC=gvaltOTREC(1,:)';
end

%% Plot
close all

wi=10;
hi=10;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1=subplot(3,1,1);
hold on;
outerpos1 = ax1.Position;
ax1.Position = [outerpos1(1)-0.07 outerpos1(2)+0.02 outerpos1(3)+0.15 outerpos1(4)+0.02];
plot(gvtimeCSET,gvaltCSET./1000,'-k','linewidth',1.5);
ylim([-0.2 15]);
xlim([gvtimeCSET(1),gvtimeCSET(end)]);
ylabel('Altitude (km)')
title('(a)              CSET Research Flight 03')

ax2=subplot(3,1,2);
hold on;
outerpos1 = ax2.Position;
ax2.Position = [outerpos1(1)-0.07 outerpos1(2)-0.01 outerpos1(3)+0.15 outerpos1(4)+0.02];
plot(gvtimeSOC,gvaltSOC./1000,'-k','linewidth',1.5);
ylim([-0.2 15]);
xlim([gvtimeSOC(1),gvtimeSOC(end)]);
ylabel('Altitude (km)')
title('(b)     SOCRATES Research Flight 11')

ax3=subplot(3,1,3);
hold on;
outerpos1 = ax3.Position;
ax3.Position = [outerpos1(1)-0.07 outerpos1(2)-0.04 outerpos1(3)+0.15 outerpos1(4)+0.02];
plot(gvtimeOTREC,gvaltOTREC./1000,'-k','linewidth',1.5);
ylim([-0.2 15]);
xlim([gvtimeOTREC(1),gvtimeOTREC(end)]);
ylabel('Altitude (km)')
title('(c)            OTREC Research Flight 01')

print(fig1, [figdir,'flightPatterns.png'],'-dpng','-r0');
saveas(fig1,[figdir,'flightPatterns.pdf'])