% Plot PID memebership function values

close all
clear all

project='socrates';

if strcmp(project,'otrec')
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
end

figdir=[indir(1:end-5),'pidPlots/'];

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

f1=figure('DefaultAxesFontSize',12,'Position',[0 300 1200 1200],'visible','on','renderer','painters');

offset=0:0.05:0.25;

% DBZ
s1=subplot(5,1,1);
hold on
for ii=1:length(fieldsIn)
    plot(dbz.(fieldsIn{ii})(:,1),dbz.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s1.YTick=[0,0.25,1,1.25];
s1.YTickLabel={'0','0','1','1'};

xlim([-40,40]);
title('DBZ');
xlabel('Reflectivity (dBZ)');

s1.XTick=-40:5:40;
grid on
box on

s1pos=s1.Position;

legend('Rain','Drizzle','Cloud drops','Mixed','Large frozen','Small frozen','orientation','horizontal','location','northoutside');
s1.Position=s1pos;


% LDR
s2=subplot(5,1,2);
hold on
for ii=1:length(fieldsIn)
    plot(ldr.(fieldsIn{ii})(:,1),ldr.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s2.YTick=[0,0.25,1,1.25];
s2.YTickLabel={'0','0','1','1'};

xlim([-30,0]);
title('LDR');
xlabel('Linear depolarization ratio (dB)');

s2.XTick=-30:2:0;
grid on
box on

% VEL
s3=subplot(5,1,3);
hold on
for ii=1:length(fieldsIn)
    plot(vel.(fieldsIn{ii})(:,1),vel.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s3.YTick=[0,0.25,1,1.25];
s3.YTickLabel={'0','0','1','1'};

xlim([-10,10]);
title('VEL');
xlabel('Radial velocity (m s^{-1})');

s3.XTick=-10:1:10;
grid on
box on

% WIDTH
s4=subplot(5,1,4);
hold on
for ii=1:length(fieldsIn)
    plot(width.(fieldsIn{ii})(:,1),width.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s4.YTick=[0,0.25,1,1.25];
s4.YTickLabel={'0','0','1','1'};

xlim([-1,1]);
title('WIDTH');
xlabel('Spectrum width (m s^{-1})');

s4.XTick=-1:0.1:1;
grid on
box on

% TEMP
s5=subplot(5,1,5);
hold on
for ii=1:length(fieldsIn)
    plot(temp.(fieldsIn{ii})(:,1),temp.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s5.YTick=[0,0.25,1,1.25];
s5.YTickLabel={'0','0','1','1'};

xlim([-60,20]);
title('TEMP');
xlabel('Temperature (C)');

s5.XTick=-60:5:20;
grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_pid_membershipCoeffs.png'],'-dpng','-r0')
