% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.1';
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));
%figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];
figdir='/scr/sci/romatsch/paperFigs/meltLayer/';

ylimits=[0 5.5];

%indir=HCRdir(project,quality,freqData);
%indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/',freqData,'/'];
indir=HCRdir(project,quality,qcVersion,freqData);

%dropsondedir=['/run/media/romatsch/RSF0006/rsf/dropsondes/',project,'/'];
dropsondedir='/scr/snow2/rsfdata/projects/socrates/dropsondes/';

dropFormat='eol';

% Make N by 2 matrix of fieldname + value type
variable_names_types = [["time", "datetime"]; ...
    ["sondeAlt", "double"]; ...
    ["meltAltMeas", "double"]; ...
    ["meltAltInt", "double"]; ...
    ["meltAltEst", "double"]; ...
    ["zeroDegAlt", "double"]];
% Make table using fieldnames & value types from above
compAlts = table('Size',[0,size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));


startTime=datetime(2018,2,7,21,20,0);
endTime=datetime(2018,2,7,23,0,0);


disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);

%% Loading dropsonde data
[dropList,dropTimes]=dropSondeList(startTime,endTime,dropsondedir,dropFormat);

[dropAlt,dropT]=getDropData(dropList,dropFormat);

%% Load HCR data

disp('Loading HCR data ...');
data=[];

data.VEL_CORR=[];
data.MELTING_LAYER=[];
data.ICING_LEVEL=[];
data.FLAG=[];
data.TEMP=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

if strcmp(freqData,'combined')
    data.DBZ=data.HCR_DBZ;
    data=rmfield(data,'HCR_DBZ');
    data.LDR=data.HCR_LDR;
    data=rmfield(data,'HCR_LDR');
    data.VEL_CORR=data.HCR_VEL;
    data=rmfield(data,'HCR_VEL');
end

%% Get indices

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);

%% Get melting layer alt and type
compAltsHourIn=array2table(nan(length(dropAlt),size(variable_names_types,1)-1),...
    'VariableNames', variable_names_types(2:end,1));
compAltsHour=cat(2,array2table(dropTimes,'VariableNames',{'time'}),compAltsHourIn);

thisMeltAlt=nan;
allSondeAlts={};
for jj=1:length(dropAlt)
    
    % Melting layer altitudes
    [minval sondeTimeInd]=min(abs(etime(datevec(data.time),datevec(dropTimes(jj)))));
    altCol=data.asl(:,sondeTimeInd);
    meltCol=data.MELTING_LAYER(:,sondeTimeInd);
    iceAlt=data.ICING_LEVEL(sondeTimeInd);
    % Melting layer alt and type
    meltIndSonde=find(meltCol==12 | meltCol==13 | meltCol==14);
    meltType=meltCol(meltIndSonde);
    if ~isempty(meltType)
        if meltCol(meltType==12)
            compAltsHour.meltAltMeas(jj)=min(altCol(meltIndSonde));
            thisMeltAlt=min(altCol(meltIndSonde));
        elseif meltCol(meltType==13)
            compAltsHour.meltAltInt(jj)=min(altCol(meltIndSonde));
            thisMeltAlt=min(altCol(meltIndSonde));
        else
            compAltsHour.meltAltEst(jj)=min(altCol(meltIndSonde));
            thisMeltAlt=min(altCol(meltIndSonde));
        end
    end
    % Zero degree alt
    zeroIndSonde=find(meltCol==11 | meltCol==21);
    if ~isempty(zeroIndSonde)
        compAltsHour.zeroDegAlt(jj)=min(altCol(zeroIndSonde));
    end
    
    % Dropsonde altitude
    tempAlt=cat(2,dropT{jj},dropAlt{jj});
    tempAlt(any(isnan(tempAlt),2),:) = [];
    signChT=diff(sign(tempAlt(:,1)));
    sondeAlts=tempAlt(signChT~=0,2);
    if ~isempty(sondeAlts) & ~isnan(thisMeltAlt)
        minDiffAlts=abs(sondeAlts-thisMeltAlt);
        compAltsHour.sondeAlt(jj)=sondeAlts(minDiffAlts==min(minDiffAlts));
    end
    allSondeAlts{end+1}=sondeAlts;
end

compAlts=cat(1,compAlts,compAltsHour);

%% Plot prepare

timeMat=repmat(data.time,size(data.VEL_CORR,1),1);

close all

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end

data.VEL_CORR(data.FLAG>1)=nan;

% Resample for plotting
newVEL=data.VEL_CORR(:,newInds);
newASL=data.asl(:,newInds);
newTEMP=data.TEMP(:,newInds);
newFindMelt=data.MELTING_LAYER(:,newInds);
newTime=data.time(newInds);

%% Plot
close all;
wi=5;
hi=6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

ax1=subplot(2,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
caxis([-10 10]);
hcb1=colorbar('XTick',-10:2:10);
view(2);
colormap(jet)

scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,7,'k','filled','MarkerEdgeColor','k');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l1=scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,7,'k','filled','MarkerEdgeColor','k');
l2=scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
l3=scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
l4=scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');


% Dropsondes
for jj=1:length(dropAlt)
    timeVec=repmat(dropTimes(jj),length(dropAlt{jj}),1);
    scatter(timeVec,dropAlt{jj}./1000,20,dropT{jj},'filled');
    %set(gca,'clim',[-10 10])
    set(gca,'colormap',jet)
    if ~isempty(allSondeAlts{jj})
        scatter(repmat(dropTimes(jj),1,length(allSondeAlts{jj})),allSondeAlts{jj}/1000,20,'k','filled');
    end
end

ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title({'(a) VEL (m s^{-1}), melting layer,';['and dropsonde temperatures (',char(176),'C)']})
grid on

ax1.Position=[0.09 0.577 0.79 0.34];
hcb1.Position=[0.895 0.58 0.04 0.34];

ax2=subplot(2,1,2)
hold on;
sub1=surf(newTime,newASL./1000,newTEMP,'edgecolor','none');
caxis([-10 10]);
hcb2=colorbar('XTick',-10:2:10);
view(2);
scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,7,'k','filled','MarkerEdgeColor','k');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l1=scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,7,'k','filled','MarkerEdgeColor','k');
l2=scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
l3=scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
l4=scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

% Dropsondes
for jj=1:length(dropAlt)
    timeVec=repmat(dropTimes(jj),length(dropAlt{jj}),1);
    scatter(timeVec,dropAlt{jj}./1000,20,dropT{jj},'filled');
    set(gca,'clim',[-10 10])
    set(gca,'colormap',jet)
    if ~isempty(allSondeAlts{jj})
        scatter(repmat(dropTimes(jj),1,length(allSondeAlts{jj})),allSondeAlts{jj}/1000,20,'k','filled');
    end
end

ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title({['(b) TEMP and dropsonde temperatures (',char(176),'C),'];'and melting layer'})
grid on

ax2.Position=[0.09 0.08 0.79 0.34];
hcb2.Position=[0.9 0.08 0.04 0.34];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'temperatures'],'-dpng','-r0');
