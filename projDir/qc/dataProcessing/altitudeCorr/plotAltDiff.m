% Plot altitude difference. Run for field data and then for corrected.
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='noreaster'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz';
qcVersion='moments';

formatOut = 'yyyymmdd_HHMM';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

if strcmp(project,'socrates')
    gvfiles=dir('/scr/snow2/rsfdata/projects/socrates/GV/RF*.nc');
elseif strcmp(project,'otrec')
    gvfiles=dir('/scr/sleet2/rsfdata/projects/otrec/GV/LRT/RF*.nc');
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/moments/10hz/';
elseif strcmp(project,'cset')
    gvfiles=dir('/scr/snow2/rsfdata/projects/cset/GV/RF*.nc');
elseif strcmp(project,'spicule')
    gvfiles=dir('/scr/sleet2/rsfdata/projects/spicule/GV/*rf*.nc');
elseif strcmp(project,'noreaster')
    gvfiles=dir('/scr/snow2/rsfdata/projects/noreaster/gv/*rf*.nc');
    indir='/scr/snow2/rsfdata/projects/noreaster/hcr/qc2/cfradial/moments/10hz/';
    %indir='/scr/snow2/rsfdata/projects/noreaster/hcr/cfradial/moments/orig/10hz/';
end

%% HCR

caseList = table2array(readtable(infile));

%indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'checkAltCorr/'];

for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    disp('Loading HCR data ...');
    
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
    
    %% GV
    
    disp('Loading GV data ...');
    
    gvfile=[gvfiles(ii).folder,'/',gvfiles(ii).name];
    
    gvtimeIn=ncread(gvfile,'Time');
    refTimegv=datetime(year(data.time(1)),month(data.time(1)),day(data.time(1)),0,0,0);
    gvtime=refTimegv+seconds(gvtimeIn);
    gvalt=ncread(gvfile,'GGALT');
    
    if min(size(gvalt))~=1
        gvalt=gvalt(1,:)';
    end
    
    %% Synchornize HCR and GV
    TThcr=timetable(data.time',data.altitude');
    TTgv=timetable(gvtime,gvalt);
    
    TT=synchronize(TThcr,TTgv,'first','linear');
    
    %% Plot
    close all
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 2000 600]);
    
    hold on
    plot(TT.Time,TT.Var1-TT.gvalt,'linewidth',2);
    xlim([data.time(1) data.time(end)]);
    ylim([-50 10]);
    ylabel('Altitude diff (m)');
    title([project,' flight ',num2str(ii)])
    grid on
    legend('HCR-GV')
    
    print([figdir,project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_altDiff'],'-dpng','-r0')
end