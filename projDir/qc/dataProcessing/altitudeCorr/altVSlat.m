% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made

project='otrec'; %socrates, aristo, cset, otrec
quality='raw'; %raw, qc1, or qc2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

outName=[project,'_',quality];

directories.figdir=['/h/eol/romatsch/hcrCalib/latVSalt/',outName,'/'];

directories.planeDir=['/scr/snow2/rsfdata/projects/',project,'/GV/'];

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

if strcmp(project,'otrec')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/'; % raw data
    %elseif strcmp(quality,'qc1')
    %    directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    %elseif strcmp(quality,'qc2')
    %    directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/'; %Final
    end
    planeFiles=dir([directories.planeDir,'OTREC*.nc']);
end

if strcmp(project,'socrates')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data
    elseif strcmp(quality,'qc1')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    elseif strcmp(quality,'qc2')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/'; %Final
    end
    planeFiles=dir([directories.planeDir,'RF*.nc']);
end

if strcmp(project,'aristo')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/'; % raw data
    elseif strcmp(quality,'qc2')
        directories.dataDir=''; %Final
    end
    planeFiles=dir([directories.planeDir,'RF*.nc']);
end

if strcmp(project,'cset')
    oceanTemp=25; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/'; %raw data
    elseif strcmp(quality,'qc1')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/moments/10hz/'; %qc data
    elseif strcmp(quality,'qc2')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/moments/10hz/'; %qc data
    end
    planeFiles=dir([directories.planeDir,'RF*.nc']);
end

directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5/',project,'/'];

%File with start date, end date
infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

% Data files with plane data
directories.planeDir=['/scr/snow2/rsfdata/projects/',project,'/GV/'];

%% Run processing

allData=[];

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    data.none=[];
    
    fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    %% Get INS data from plane
    planeFileIn=[planeFiles(ii).folder,'/',planeFiles(ii).name];
    
    planeAlt=ncread(planeFileIn,'GGALT');
    
    if min(size(planeAlt))~=1
        planeAlt=planeAlt(1,:)';
    end
    
    planeTimeOrig=ncread(planeFileIn,'Time');
    
    info=ncinfo(planeFileIn);
    if strcmp(project,'cset')
        planeRefTimeIn=info.Variables(133).Attributes(3).Value;
    else
        planeRefTimeIn=info.Variables(1).Attributes(3).Value;
    end
    planeRefTime=datetime(str2num(planeRefTimeIn(15:18)),str2num(planeRefTimeIn(20:21)),str2num(planeRefTimeIn(23:24)),...
        str2num(planeRefTimeIn(26:27)),str2num(planeRefTimeIn(29:30)),str2num(planeRefTimeIn(32:33)));
    
    planeTime=planeRefTime+seconds(planeTimeOrig);
    
    TThcr=timetable(data.time',data.altitude',data.latitude');
    TThcr=rmmissing(TThcr);
    TThcr=sortrows(TThcr);
    
    TTplane=timetable(planeTime,planeAlt);
    TTplane=rmmissing(TTplane);
    TTplane=sortrows(TTplane);
    
    TTcombined=synchronize(TThcr,TTplane,'last','nearest');
    TTcombined.Properties.VariableNames{'Var1'} = 'hcrAlt';
    TTcombined.Properties.VariableNames{'Var2'} = 'lat';
    
    timeInds=find(TTcombined.Time>=data.time(1) & TTcombined.Time<=data.time(end));
    TTcombined=TTcombined(timeInds,:);
    
    TTcombined.altDiff=TTcombined.planeAlt-TTcombined.hcrAlt;
    
    allData=cat(1,allData,TTcombined);
    
    %% Plot
    close all
    f1 = figure('Position',[200 500 2000 1000],'DefaultAxesFontSize',12);
    
    subplot(2,1,1)
    hold on
    l1=plot(TTcombined.Time,TTcombined.planeAlt./1000,'-k','linewidth',2);
    l2=plot(TTcombined.Time,TTcombined.hcrAlt./1000,'-b','linewidth',2);
    ylabel('Altitude (km)');
    ylim([0 10]);
    xlim([TTcombined.Time(1),TTcombined.Time(end)]);
    grid on
    
    legend([l1 l2],{'GGALT','HCRins'},'location','northeast');
    
    title([project,' Flight ',num2str(ii),' ',datestr(data.time(1)),' to ',datestr(data.time(end))])
    
    subplot(2,1,2)
    hold on
    l1=plot(TTcombined.Time,TTcombined.altDiff,'-k','linewidth',2);
    ylabel('Altitude (m)');
    ylim([0 80]);
    xlim([TTcombined.Time(1),TTcombined.Time(end)]);
    grid on
    
    yyaxis right
    l2=plot(TTcombined.Time,TTcombined.lat,'-b','linewidth',2);
    ylabel('Latitude (deg)');
    %ylim([0 10]);
    xlim([TTcombined.Time(1),TTcombined.Time(end)]);
    grid on
    
    legend([l1 l2],{'GGALT-HCRins','Latitude'},'location','northeast');
    
    set(gcf,'PaperPositionMode','auto')
    print([directories.figdir,project,'_Flight',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    clear data
end
%% Scatter plot
% Remove outliers
TF=isoutlier(allData.altDiff);
allData(TF,:)=[];

if strcmp(project,'otrec')
    TF=isoutlier(allData.altDiff);
    allData(TF,:)=[];
end

%Remove empty data
nanInds=find(isnan(allData.altDiff));
allData(nanInds,:)=[];
nanInds=find(isnan(allData.lat));
allData(nanInds,:)=[];

f2 = figure('Position',[200 500 800 800],'DefaultAxesFontSize',12);

scatter(allData.lat,allData.altDiff);

hold on
xlimits=xlim;
ylimits=ylim;

fitOrth=gmregress(allData.lat,allData.altDiff,1);
fitAll=[fitOrth(2) fitOrth(1)];
xFit = xlimits(1):0.1:xlimits(2);
yFit = polyval(fitAll, xFit);

plot(xFit, yFit,'-k','linewidth',3);

xlim(xlimits);
ylim(ylimits);

xlabel('Latitude [deg]');
ylabel('GGALT-HCRins [m]');
title([project,' latitude vs altitude difference']);

textbp(['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
%ylim([0 55])

set(f2,'PaperPositionMode','auto')
print([directories.figdir,project,'_latVSaltdiff'],'-dpng','-r0')

save([directories.figdir,project,'_table'],'allData');