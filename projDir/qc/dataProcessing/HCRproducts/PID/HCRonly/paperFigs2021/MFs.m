% Plot PID memebership function values

close all
clear all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/paperFigs/'];

dbz.rain=[-100,0;3,0;5,1;100,1];
dbz.drizzle=[-100,0;-18,0;-16,1;5,1;8,0;100,0];
dbz.cloud=[-100,1;-16,1;-14,0;100,0];
dbz.mixed=[-100,0;-3,0;-1,1;20,1;25,0;100,0];
dbz.lfrozen=[-100,0;7,0;9,1;18,1;20,0;100,0];
dbz.sfrozen=[-100,0;-25,0;-20,1;9,1;11,0;100,0];

ldr.rain=[-100,1;-27,1;-22,0;100,0];
ldr.drizzle=[-100,1;-27,1;-25,0;100,0];
ldr.cloud=[-100,1;-27,1;-25,0;100,0];
ldr.mixed=[-100,0;-20,0;-17,1;-8,1;-6,0;100,0];
ldr.lfrozen=[-100,0;-22,0;-20,1;-16,1;-14,0;100,0];
ldr.sfrozen=[-100,0;-26,0;-24,1;-15,1;-12,0;100,0];

vel.rain=[-100,0;2,0;3,1;100,1];
vel.drizzle=[-100,0;0,0;0.5,1;1,1;2,0;100,0];
vel.cloud=[-100,1;1,1;2,0;100,0];
vel.mixed=[-100,0;0.5,0;1,1;3,1;4,0;100,0];
vel.lfrozen=[-100,0;0.8,0;1,1;2.5,1;3.5,0;100,0];
vel.sfrozen=[-100,0;-1,0;0,1;1,1;2,0;100,0];

width.rain=[-100,0;0.1,0;0.2,1;100,1];
width.drizzle=[-100,1;0.2,1;0.3,0;100,0];
width.cloud=[-100,1;0.1,1;0.2,0;100,0];
width.mixed=[-100,0;0.2,0;0.3,1;100,1];
width.lfrozen=[-100,0;0.2,0;0.3,1;100,1];
width.sfrozen=[-100,1;0.7,1;0.9,0;100,0];

temp.rain=[-100,0;-2,0;2,1;100,1];
temp.drizzle=[-100,0;-50,0;-39,1;100,1];
temp.cloud=[-100,1;100,1];
temp.mixed=[-100,0;-2,0;0,1;3,1;6,0;100,0];
temp.lfrozen=[-100,1;0,1;6,0;100,0];
temp.sfrozen=[-100,1;-1,1;5,0;100,0];

%% Plot

close all

fieldsIn=fields(dbz);

colorsPlot=[1,0,0;0,1,0;0,0,1;0.5,0,0;1,1,0;0,1,1];

wi=8;
hi=6.5;

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
s1=subplot(6,1,1);
hold on
for ii=1:length(fieldsIn)
    plot(dbz.(fieldsIn{ii})(:,1),dbz.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s1.YTick=[0,0.25,1,1.25];
s1.YTickLabel={'0','0','1','1'};

xlim([-40,40]);
xlabel('Reflectivity (dBZ)');[0.025 0.1 0.95 0.13];

s1.XTick=-40:5:40;
text(-39,0.6,'(a) DBZ',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% VEL
s2=subplot(6,1,2);
hold on
for ii=1:length(fieldsIn)
    plot(vel.(fieldsIn{ii})(:,1),vel.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s2.YTick=[0,0.25,1,1.25];
s2.YTickLabel={'0','0','1','1'};

xlim([-4,6]);
xlabel('Radial velocity (m s^{-1})');

s2.XTick=-10:0.5:10;
text(-3.88,0.6,'(b) VEL',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% LDR
s3=subplot(6,1,3);
hold on
for ii=1:length(fieldsIn)
    plot(ldr.(fieldsIn{ii})(:,1),ldr.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s3.YTick=[0,0.25,1,1.25];
s3.YTickLabel={'0','0','1','1'};

xlim([-32,0]);
xlabel('Linear depolarization ratio (dB)');

s3.XTick=-32:2:0;
text(-31.5,0.6,'(c) LDR',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% WIDTH
s4=subplot(6,1,4);
hold on
for ii=1:length(fieldsIn)
    plot(width.(fieldsIn{ii})(:,1),width.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s4.YTick=[0,0.25,1,1.25];
s4.YTickLabel={'0','0','1','1'};

xlim([-0.1,1]);
xlabel('Spectrum width (m s^{-1})');

s4.XTick=-1:0.1:1;
text(-0.08,0.6,'(d) WIDTH',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on

% TEMP
s5=subplot(6,1,5);
hold on
for ii=1:length(fieldsIn)
    plot(temp.(fieldsIn{ii})(:,1),temp.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s5.YTick=[0,0.25,1,1.25];
s5.YTickLabel={'0','0','1','1'};

xlim([-60,20]);
text(-58.5,0.6,'(e) TEMP',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
xlabel(['Temperature (',char(176),'C)']);

s5.XTick=-60:5:20;
grid on
box on

% Dummy for legend
s6=subplot(6,1,6);
hold on
for ii=1:length(fieldsIn)
    plot(dbz.(fieldsIn{ii})(:,1),dbz.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end

legend('Rain','Drizzle','Cloud drops','Mixed phase','Large frozen','Small frozen',...
    'orientation','horizontal','location','northoutside','FontSize',10);

s1.Position=[0.025 0.87 0.95 0.11];

s2.Position=[0.025 0.68 0.95 0.11];

s3.Position=[0.025 0.49 0.95 0.11];

s4.Position=[0.025 0.3 0.95 0.11];

s5.Position=[0.025 0.11 0.95 0.11];

s6.Position=[0.025 -0.136 0.95 0.13];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'MFs.png'],'-dpng','-r0')
