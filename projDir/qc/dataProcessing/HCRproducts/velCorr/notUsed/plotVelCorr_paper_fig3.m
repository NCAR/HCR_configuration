% Compare data from before and after velocity correction

clear all
close all

fig=3;
savefig=1;

addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/utils/colormaps/');

formatOut = 'yyyymmdd_HHMM';
figdir='/h/eol/romatsch/hcrCalib/velCorr/velFigs_paper/';

startTime=datetime(2015,2,2,19,15,0);
    endTime=datetime(2015,2,2,19,20,0);
xtoff=15;
ytoff=0.05;
minEdge1=-1.5;
maxEdge1=1.5;

indir='/scr/eldora1/rsfdata/noreaster/cfradial/moments/qcv1/10hz/';

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if ~isempty(fileList)
    
    data.time=[];
       
    % Get uncorrected data
    for ii=1:size(fileList,2)
        infile=fileList{ii};
        
        startTimeIn=ncread(infile,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(infile,'time')';
        timeIn=startTimeFile+seconds(timeRead);
        data.time=[data.time,timeIn];
    end
    
    timeInd=find(data.time>=startTime & data.time<=endTime);
    data.time=data.time(timeInd);

load([figdir,'Fig3_data.mat'])

%% Plot

close all

hm=2;
wm=2;
wi=10;
hi=4;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

ax1=subtightplot(1,1,1,[hm,wm]);
hold on;
outerpos1 = ax1.Position;
ax1.Position = [outerpos1(1)+0.03 outerpos1(2)+0.07 outerpos1(3)-0.04 outerpos1(4)-0.04];

plot(data.time,VEL_ground','-b','linewidth',0.8);
plot(data.time,VG_fir_tofilter','-r','linewidth',1);
plot(data.time,VG_fir_final(3,:)','-k','linewidth',1.5);
xlim([startTime,endTime]);
ylim([minEdge1 maxEdge1]);
ylabel('V_r [m s^{-1}]');

legend('v_{surf}^{meas}','Intermediate','Filtered');

if savefig
    set(gcf,'PaperPositionMode','auto')
    print(fig1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_fig',num2str(fig)],'-dpng','-r0');
end

end