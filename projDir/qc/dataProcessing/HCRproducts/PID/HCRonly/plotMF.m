% Plot PID memebership function values

close all
clear all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

indir=HCRdir(project,quality,qcVersion,freqData);

% if strcmp(project,'otrec')
%     indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
% elseif strcmp(project,'socrates')
%     indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
% end

memCoeffs

figdir=[indir(1:end-5),'pidPlots/'];

dbz1.rain=[-100,0;dbz.rain(1),0;dbz.rain(2),1;100,1];
dbz1.drizzle=[-100,0;dbz.drizzle(1),0;dbz.drizzle(2),1;dbz.drizzle(3),1;dbz.drizzle(4),0;100,0];
dbz1.cloud=[-100,1;dbz.cloud(1),1;dbz.cloud(2),0;100,0];
dbz1.mixed=[-100,0;dbz.mixed(1),0;dbz.mixed(2),1;100,1];
dbz1.lfrozen=[-100,0;dbz.lfrozen(1),0;dbz.lfrozen(2),1;100,1];
dbz1.sfrozen=[-100,1;dbz.sfrozen(1),1;dbz.sfrozen(2),0;100,0];

ldr1.rain=[-100,1;ldr.rain(1),1;ldr.rain(2),0;100,0];
ldr1.drizzle=[-100,1;ldr.drizzle(1),1;ldr.drizzle(2),0;100,0];
ldr1.cloud=[-100,1;ldr.cloud(1),1;ldr.cloud(2),0;100,0];
ldr1.mixed=[-100,0;ldr.mixed(1),0;ldr.mixed(2),1;ldr.mixed(3),1;ldr.mixed(4),0;100,0];
ldr1.lfrozen=[-100,0;ldr.lfrozen(1),0;ldr.lfrozen(2),1;ldr.lfrozen(3),1;ldr.lfrozen(4),0;100,0];
ldr1.sfrozen=[-100,0;ldr.sfrozen(1),0;ldr.sfrozen(2),1;ldr.sfrozen(3),1;ldr.sfrozen(4),0;100,0];

vel1.rain=[-100,0;vel.rain(1),0;vel.rain(2),1;100,1];
vel1.drizzle=[-100,0;vel.drizzle(1),0;vel.drizzle(2),1;vel.drizzle(3),1;vel.drizzle(4),0;100,0];
vel1.cloud=[-100,1;vel.cloud(1),1;vel.cloud(2),0;100,0];
vel1.mixed=[-100,0;vel.mixed(1),0;vel.mixed(2),1;100,1];
vel1.lfrozen=[-100,0;vel.lfrozen(1),0;vel.lfrozen(2),1;100,1];
vel1.sfrozen=[-100,0;vel.sfrozen(1),0;vel.sfrozen(2),1;vel.sfrozen(3),1;vel.sfrozen(4),0;100,0];

width1.rain=[-100,0;width.rain(1),0;width.rain(2),1;100,1];
width1.drizzle=[-100,0;width.drizzle(1),0;width.drizzle(2),1;100,1];
width1.cloud=[-100,0;width.cloud(1),0;width.cloud(2),1;100,1];
width1.mixed=[-100,1;width.mixed(1),1;width.mixed(2),0;100,0];
width1.lfrozen=[-100,1;width.lfrozen(1),1;width.lfrozen(2),0;100,0];
width1.sfrozen=[-100,1;width.sfrozen(1),1;width.sfrozen(2),0;100,0];

temp1.rain=[-100,0;temp.rain(1),0;temp.rain(2),1;100,1];
temp1.drizzle=[-100,1;100,1];
temp1.cloud=[-100,1;100,1];
temp1.mixed=[-100,0;temp.mixed(1),0;temp.mixed(2),1;temp.mixed(3),1;temp.mixed(4),0;100,0];
temp1.lfrozen=[-100,1;temp.lfrozen(1),1;temp.lfrozen(2),0;100,0];
temp1.sfrozen=[-100,1;temp.sfrozen(1),1;temp.sfrozen(2),0;100,0];

%% Plot

close all

fieldsIn=fields(dbz1);

colorsPlot=[1,0,0;0,1,0;0,0,1;0.5,0,0;1,1,0;0,1,1];

f1=figure('DefaultAxesFontSize',12,'Position',[0 300 1200 1200],'visible','on','renderer','painters');

offset=0:0.05:0.25;

% DBZ
s1=subplot(5,1,1);
hold on
for ii=1:length(fieldsIn)
    plot(dbz1.(fieldsIn{ii})(:,1),dbz1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
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
    plot(ldr1.(fieldsIn{ii})(:,1),ldr1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
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
    plot(vel1.(fieldsIn{ii})(:,1),vel1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s3.YTick=[0,0.25,1,1.25];
s3.YTickLabel={'0','0','1','1'};

xlim([-4,5]);
title('VEL');
xlabel('Radial velocity (m s^{-1})');

s3.XTick=-10:0.5:10;
grid on
box on

% WIDTH
s4=subplot(5,1,4);
hold on
for ii=1:length(fieldsIn)
    plot(width1.(fieldsIn{ii})(:,1),width1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
end
ylim([-0.25,1.5]);
s4.YTick=[0,0.25,1,1.25];
s4.YTickLabel={'0','0','1','1'};

xlim([-0.5,1]);
title('WIDTH');
xlabel('Spectrum width (m s^{-1})');

s4.XTick=-1:0.1:1;
grid on
box on

% TEMP
s5=subplot(5,1,5);
hold on
for ii=1:length(fieldsIn)
    plot(temp1.(fieldsIn{ii})(:,1),temp1.(fieldsIn{ii})(:,2)+offset(ii),'-','color',colorsPlot(ii,:),'linewidth',2);
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
print(f1,[figdir,project,'_pid_membershipCoeffs_LDR.png'],'-dpng','-r0')
