% Remove stripes in HX

clear all;
close all;

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.1_full/hxPlots/'];

dataDir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/txt/';
infile='SPICULE.PsVoltage.txt';

indata=table2array(readtable([dataDir,infile]));

PsVoltage=indata(:,10);

addSecFrac=repmat([0;0.5],ceil(length(PsVoltage)/2),1);
addSecFrac=addSecFrac(1:length(PsVoltage));

time=datetime(1970,1,1)+seconds(indata(:,8));
time=time+seconds(addSecFrac);

%% Plot
close all

figure('DefaultAxesFontSize',11,'position',[1,100,1500,600],'renderer','painters');

plot(time,PsVoltage,'-k','LineWidth',1);
grid on
ylabel('PsVoltage')
xlim([time(1),time(end)])

formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
print([figdir,datestr(time(1),formatOut),'_to_',datestr(time(end),formatOut),'_PsVoltageHighRes'],...
    '-dpng','-r0');
