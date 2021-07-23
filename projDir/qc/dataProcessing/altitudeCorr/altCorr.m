% Correct HCR altitude
clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz';

figdir=['/h/eol/romatsch/hcrCalib/latVSalt/altCorr/',project,'/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];
%infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'_testFlights.txt'];

if strcmp(project,'socrates')
    gvfiles=dir('/scr/snow2/rsfdata/projects/socrates/GV/RF*.nc');
elseif strcmp(project,'otrec')
    gvfiles=dir('/scr/snow1/rsfdata/projects/otrec/GV/OTRECrf*.nc');
    
    %gvfiles1=dir('/scr/snow2/rsfdata/projects/otrec/GV/OTRECtf*.nc');
    %gvfiles2=dir('/scr/snow2/rsfdata/projects/otrec/GV/OTRECff*.nc');
    %gvfiles=cat(1,gvfiles1,gvfiles2);
elseif strcmp(project,'cset')
    gvfiles=dir('/scr/snow2/rsfdata/projects/cset/GV/RF*.nc');
end

%% Alt Corr file
altFile='/scr/sci/romatsch/data/topo/EGM_2008_WGS84_2.5minx2.5min.nc';

altLat=ncread(altFile,'y0');
altLon=ncread(altFile,'x0');

altHt=ncread(altFile,'GeoidHt');

altLatMat=repmat(altLat',size(altHt,1),1);
altLonMat=repmat(altLon,1,size(altHt,2));
%% HCR

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freqData);

for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    data.dummy=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isempty(data.time)
        disp('No data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars)
        if ~isfield(data,dataVars{kk})
            dataVars{kk}=[];
        end
    end
    
    %% Corr
    HtInt = interpn(altLonMat,altLatMat,altHt,data.longitude,data.latitude);
    
    %% GV
    
    gvfile=[gvfiles(ii).folder,'/',gvfiles(ii).name];
    
    gvtimeIn=ncread(gvfile,'Time');
    refTimegv=datetime(year(data.time(1)),month(data.time(1)),day(data.time(1)),0,0,0);
    gvtime=refTimegv+seconds(gvtimeIn);
    gvalt=ncread(gvfile,'GGALT');
    
    if min(size(gvalt))~=1
        gvalt=gvalt(1,:)';
    end
    
    %% Synchornize HCR and GV
    TThcr=timetable(data.time',data.altitude',HtInt');
    TTgv=timetable(gvtime,gvalt);
    
    TT=synchronize(TThcr,TTgv,'first','linear');
    
    %% Plot
    close all
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 2000 600]);
    
    hold on
    plot(TT.Time,TT.Var1-TT.gvalt,'linewidth',2);
    plot(TT.Time,(TT.Var1-TT.Var2)-TT.gvalt,'linewidth',2);
    xlim([data.time(1) data.time(end)]);
    ylim([-50 10]);
    ylabel('Altitude diff (m)');
    title([project,' flight ',num2str(ii)])
    grid on
    legend('HCR-GV','HCR-Corr-GV')
    
    print([figdir,project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_altCorr'],'-dpng','-r0')
end