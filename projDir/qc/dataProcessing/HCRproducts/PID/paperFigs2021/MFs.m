% Plot PID memebership function values

close all
clear all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/paperFigs/'];

memCoeffsComb

dbz1.rain=[-100,0;dbz.rain(1),0;dbz.rain(2),1;100,1];
dbz1.drizzle=[-100,0;dbz.drizzle(1),0;dbz.drizzle(2),1;dbz.drizzle(3),1;dbz.drizzle(4),0;100,0];
dbz1.cloud=[-100,1;dbz.cloud(1),1;dbz.cloud(2),0;100,0];
dbz1.mixed=[-100,1;100,1];
dbz1.lfrozen=[-100,0;dbz.lfrozen(1),0;dbz.lfrozen(2),1;100,1];
dbz1.sfrozen=[-100,1;dbz.sfrozen(1),1;dbz.sfrozen(2),0;100,0];

hcrldr1.rain=[-100,1;hcrldr.rain(1),1;hcrldr.rain(2),0;100,0];
hcrldr1.drizzle=[-100,1;hcrldr.drizzle(1),1;hcrldr.drizzle(2),0;100,0];
hcrldr1.cloud=[-100,1;hcrldr.cloud(1),1;hcrldr.cloud(2),0;100,0];
hcrldr1.mixed=[-100,0;hcrldr.mixed(1),0;hcrldr.mixed(2),1;hcrldr.mixed(3),1;hcrldr.mixed(4),0;100,0];
hcrldr1.lfrozen=[-100,0;hcrldr.lfrozen(1),0;hcrldr.lfrozen(2),1;hcrldr.lfrozen(3),1;hcrldr.lfrozen(4),0;100,0];
hcrldr1.sfrozen=[-100,0;hcrldr.sfrozen(1),0;hcrldr.sfrozen(2),1;hcrldr.sfrozen(3),1;hcrldr.sfrozen(4),0;100,0];

vel1.rain=[-100,0;vel.rain(1),0;vel.rain(2),1;100,1];
vel1.drizzle=[-100,0;vel.drizzle(1),0;vel.drizzle(2),1;vel.drizzle(3),1;vel.drizzle(4),0;100,0];
vel1.cloud=[-100,1;vel.cloud(1),1;vel.cloud(2),0;100,0];
vel1.mixed=[-100,0;vel.mixed(1),0;vel.mixed(2),1;100,1];
vel1.lfrozen=[-100,0;vel.lfrozen(1),0;vel.lfrozen(2),1;100,1];
vel1.sfrozen=[-100,1;vel.sfrozen(1),1;vel.sfrozen(2),0;100,0];

temp1.rain=[-100,0;temp.rain(1),0;temp.rain(2),1;100,1];
temp1.drizzle=[-100,0;temp.drizzle(1),0;temp.drizzle(2),1;100,1];
temp1.cloud=[-100,0;temp.cloud(1),0;temp.cloud(2),1;100,1];
temp1.mixed=[-100,0;temp.mixed(1),0;temp.mixed(2),1;temp.mixed(3),1;temp.mixed(4),0;100,0];
temp1.lfrozen=[-100,1;temp.lfrozen(1),1;temp.lfrozen(2),0;100,0];
temp1.sfrozen=[-100,1;temp.sfrozen(1),1;temp.sfrozen(2),0;100,0];

back1.rain=[1e-9,0;back.rain(1),0;back.rain(2),1;100,1];
back1.drizzle=[1e-9,0;back.drizzle(1),0;back.drizzle(2),1;100,1];
back1.cloud=[1e-9,0;back.cloud(1),0;back.cloud(2),1;100,1];
back1.mixed=[1e-9,0;100,0];
back1.lfrozen=[1e-9,1;back.lfrozen(1),1;back.lfrozen(2),0;100,0];
back1.sfrozen=[1e-9,1;back.sfrozen(1),1;back.sfrozen(2),0;100,0];

hsrlldr1.rain=[-100,1;hsrlldr.rain(1),1;hsrlldr.rain(2),0;100,0];
hsrlldr1.drizzle=[-100,1;hsrlldr.drizzle(1),1;hsrlldr.drizzle(2),0;100,0];
hsrlldr1.cloud=[-100,1;hsrlldr.cloud(1),1;hsrlldr.cloud(2),0;100,0];
hsrlldr1.mixed=[-100,1;hsrlldr.mixed(1),1;hsrlldr.mixed(2),0;100,0];
hsrlldr1.lfrozen=[-100,0;hsrlldr.lfrozen(1),0;hsrlldr.lfrozen(2),1;100,1];
hsrlldr1.sfrozen=[-100,0;hsrlldr.sfrozen(1),0;hsrlldr.sfrozen(2),1;100,1];


%% Plot

close all

fieldsIn=fields(dbz);

%colorsPlot=[1,0,0;0,1,0;0,0,1;0.5,0,0;1,1,0;0,1,1];
colorsPlot=[255,0,0; 249,163,25; 136,34,185; 17,170,51; 0,0,255; 0,255,255];
colorsPlot=colorsPlot./255;

wi=8;
hi=9;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

offset=0:0.05:0.25;

% DBZ
s1=subplot(7,1,1);
hold on
for ii=1:length(fieldsIn)
    plot(dbz1.(fieldsIn{ii})(:,1),dbz1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s1.YTick=[0,0.25,1,1.25];
s1.YTickLabel={'0','0','1','1'};

xlim([-40,40]);
xlabel('Reflectivity (dBZ)');[0.025 0.1 0.95 0.13];

s1.XTick=-40:5:40;
text(-39,0.6,'(a) HCR DBZ',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% VEL
s2=subplot(7,1,2);
hold on
for ii=1:length(fieldsIn)
    plot(vel1.(fieldsIn{ii})(:,1),vel1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s2.YTick=[0,0.25,1,1.25];
s2.YTickLabel={'0','0','1','1'};

xlim([-4,6]);
xlabel('Radial velocity (m s^{-1})');

s2.XTick=-10:0.5:10;
text(-3.88,0.6,'(b) HCR VEL',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% LDR
s3=subplot(7,1,3);
hold on
for ii=1:length(fieldsIn)
    plot(hcrldr1.(fieldsIn{ii})(:,1),hcrldr1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s3.YTick=[0,0.25,1,1.25];
s3.YTickLabel={'0','0','1','1'};

xlim([-34,0]);
xlabel('Linear depolarization ratio (dB)');

s3.XTick=-32:2:0;
text(-33.5,0.6,'(c) HCR LDR',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% PLDR
s4=subplot(7,1,4);
hold on
for ii=1:length(fieldsIn)
    plot(hsrlldr1.(fieldsIn{ii})(:,1),hsrlldr1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s4.YTick=[0,0.25,1,1.25];
s4.YTickLabel={'0','0','1','1'};

xlim([-0.5,0.5]);
xlabel('Particle linear depolarization ratio');

s4.XTick=-1:0.1:1;
text(-0.485,0.6,'(d) HSRL LLDR',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% BACKSCATTER
s5=subplot(7,1,5);
hold on
for ii=1:length(fieldsIn)
    plot(back1.(fieldsIn{ii})(:,1),back1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s5.YTick=[0,0.25,1,1.25];
s5.YTickLabel={'0','0','1','1'};
set(gca,'xscale','log')

xlim([1e-8,0.002]);
text(1.2e-8,0.83,'(e) HSRL BACKSCATTER',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
xlabel('Aerosol backscatter coefficient (m^{-1} sr^{-1})');

%s5.XTick=-10:0.5:10;
grid on
box on

% TEMP
s6=subplot(7,1,6);
hold on
for ii=1:length(fieldsIn)
    plot(temp1.(fieldsIn{ii})(:,1),temp1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s6.YTick=[0,0.25,1,1.25];
s6.YTickLabel={'0','0','1','1'};

xlim([-60,20]);
text(-58.5,0.87,'(f) TEMP',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
xlabel(['Temperature (',char(176),'C)']);

s6.XTick=-60:5:20;
grid on
box on

% Dummy for legend
s7=subplot(7,1,7);
hold on
for ii=1:length(fieldsIn)
    plot(dbz1.(fieldsIn{ii})(:,1),dbz1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end

legend('Rain','Drizzle','Cloud drops','Melting','Large frozen','Small frozen',...
    'orientation','horizontal','location','northoutside','FontSize',10);

s1.Position=[0.025 0.88 0.95 0.11];

s2.Position=[0.025 0.72 0.95 0.11];

s3.Position=[0.025 0.56 0.95 0.11];

s4.Position=[0.025 0.4 0.95 0.11];

s5.Position=[0.025 0.24 0.95 0.11];

s6.Position=[0.025 0.074 0.95 0.11];

s7.Position=[0.025 -0.14 0.95 0.13];

set(gcf,'PaperPositionMode','auto')
%print(fig1,[figdir,'MFs.png'],'-dpng','-r0')
print(fig1,[figdir,'MFs.tif'],'-dtiffn','-r0')