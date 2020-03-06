% Compare data from before and after velocity correction

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

indir=HCRdir(project,quality,freqData);

maxEdge=10; % Upper edge for plotting
color_map=colormap(vel_default(29));

limits=-4.05:0.3:4.05;
limits=[-inf limits inf];

startTime=datetime(2019,8,16,18,18,0);
endTime=datetime(2019,8,16,18,28,0);

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data.VEL=[];
data.VEL_RAW=[];
data.VEL_CORR=[];

dataVars=fieldnames(data);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Plot vel field
close all

wi=10;
hi=10;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%% VEL_RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1=subplot(3,1,1);
hold on;
outerpos1 = ax1.Position;
ax1.Position = [outerpos1(1)-0.07 outerpos1(2)+0.02 outerpos1(3)+0.14 outerpos1(4)+0.02];

fig2=surf(data.time,data.asl./1000,data.VEL_RAW);
fig2.EdgeColor='none';
ylim([-0.2 maxEdge]);
xlim([startTime,endTime]);
view(2);

fld=fig2.CData;

col_def1 = nan(size(fld));
col_def2 = nan(size(fld));
col_def3 = nan(size(fld));

for ii=1:size(color_map,1)
    col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
    col_def1(col_ind)=color_map(ii,1);
    col_def2(col_ind)=color_map(ii,2);
    col_def3(col_ind)=color_map(ii,3);
end
if ~isequal(size(col_def1),(size(fld)))
    col_def=cat(3,col_def1',col_def2',col_def3');
else
    col_def=cat(3,col_def1,col_def2,col_def3);
end
fig2.CData=col_def;

hcb=colorbar;
set(get(hcb,'Title'),'String','m/s');
colormap(gca,color_map);
caxis([0 size(color_map,1)]);
caxis_ytick_labels=num2str((-4:1:4)');
caxis_yticks=0:3.6:32;
set(hcb,'ytick',caxis_yticks);
set(hcb,'YTickLabel',caxis_ytick_labels);
ylabel('Altitude [km]');

title('(a)       VEL_RAW','interpreter','none');

%%%%%%%%%%%%%%% Vel %%%%%%%%%%%%%%%%%%%%%%%%%%
ax2=subplot(3,1,2);
hold on;
outerpos1 = ax2.Position;
ax2.Position = [outerpos1(1)-0.07 outerpos1(2)-0.01 outerpos1(3)+0.14 outerpos1(4)+0.02];

fig2=surf(data.time,data.asl./1000,data.VEL);
fig2.EdgeColor='none';
ylim([-0.2 maxEdge]);
xlim([startTime,endTime]);
view(2);

fld=fig2.CData;

col_def1 = nan(size(fld));
col_def2 = nan(size(fld));
col_def3 = nan(size(fld));

for ii=1:size(color_map,1)
    col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
    col_def1(col_ind)=color_map(ii,1);
    col_def2(col_ind)=color_map(ii,2);
    col_def3(col_ind)=color_map(ii,3);
end
if ~isequal(size(col_def1),(size(fld)))
    col_def=cat(3,col_def1',col_def2',col_def3');
else
    col_def=cat(3,col_def1,col_def2,col_def3);
end
fig2.CData=col_def;

hcb=colorbar;
set(get(hcb,'Title'),'String','m/s');
colormap(gca,color_map);
caxis([0 size(color_map,1)]);
caxis_ytick_labels=num2str((-4:1:4)');
caxis_yticks=0:3.6:32;
set(hcb,'ytick',caxis_yticks);
set(hcb,'YTickLabel',caxis_ytick_labels);
ylabel('Altitude [km]');

title('(b)                VEL');

%%%%%%%%%%%%%%%%%%%%%% VelCorr %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3=subplot(3,1,3);
hold on;
outerpos1 = ax3.Position;
ax3.Position = [outerpos1(1)-0.07 outerpos1(2)-0.04 outerpos1(3)+0.14 outerpos1(4)+0.02];
fig2=surf(data.time,data.asl./1000,data.VEL_CORR);
fig2.EdgeColor='none';
ylim([-0.2 maxEdge]);
xlim([startTime,endTime]);
view(2);

fld=fig2.CData;

col_def1 = nan(size(fld));
col_def2 = nan(size(fld));
col_def3 = nan(size(fld));

for ii=1:size(color_map,1)
    col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
    col_def1(col_ind)=color_map(ii,1);
    col_def2(col_ind)=color_map(ii,2);
    col_def3(col_ind)=color_map(ii,3);
end
if ~isequal(size(col_def1),(size(fld)))
    col_def=cat(3,col_def1',col_def2',col_def3');
else
    col_def=cat(3,col_def1,col_def2,col_def3);
end
fig2.CData=col_def;

hcb=colorbar;
set(get(hcb,'Title'),'String','m/s');
colormap(gca,color_map);
caxis([0 size(color_map,1)]);
caxis_ytick_labels=num2str((-4:1:4)');
caxis_yticks=0:3.6:32;
set(hcb,'ytick',caxis_yticks);
set(hcb,'YTickLabel',caxis_ytick_labels);
ylabel('Altitude [km]');

title(['(c)      VEL_CORR'],'interpreter','none');

print(fig1, [figdir,'velExample.png'],'-dpng','-r0');

