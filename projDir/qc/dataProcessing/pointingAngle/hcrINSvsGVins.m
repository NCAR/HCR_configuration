% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

upwind_limit=2;
crosswind_limit=2;

% If attitude correction was sometimes turned off we need to filter that
% data
filterAttitudeCorrOff=1;

project='otrec'; % socrates, cset, aristo, otrec
quality='field'; % field, qc1, qc2
freq='10hz';

%figdir=['/h/eol/romatsch/hcrCalib/pointAng/',project,'/wholeFlights/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

figdir=['/h/eol/romatsch/hcrCalib/pointAng/',project,'/insCompare/'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freq);

lineCols=lines;

gvfiles=dir('/scr/snow2/rsfdata/projects/otrec/GV/OTRECrf*.nc');

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
                    
        % Desired variables. The variable name comes after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        data.roll=[];
        data.pitch=[];
        %data.rotation=[];
        %data.tilt=[];
        data.drift=[];
        %data.northward_velocity=[];
        %data.eastward_velocity=[];
        %data.vertical_velocity=[];
        %data.altitude_agl=[];
        
        %data.VEL_RAW=[];
        %data.DBZ=[];
        
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
       
           
        %% Drift from GV
        close all
        
        gvfile=[gvfiles(ii).folder,'/',gvfiles(ii).name];
        
        gvtimeIn=ncread(gvfile,'Time');
        refTimegv=datetime(year(data.time(1)),month(data.time(1)),day(data.time(1)),0,0,0);
        gvdrift=ncread(gvfile,'DRFTA');
        gvroll=ncread(gvfile,'ROLL');
        gvpitch=ncread(gvfile,'PITCH');
        
        gvtime=refTimegv+seconds(gvtimeIn);
        
        hcrTT=timetable(data.time',data.pitch',data.roll',data.drift');
        gvTT=timetable(gvtime,gvpitch,gvroll,gvdrift);
        
        ttsync=synchronize(hcrTT,gvTT,'first','linear');
        
        f2=figure('DefaultAxesFontSize',12);
        set(f2,'Position',[200 500 2000 800]);
        
        hold on
        plot(data.time,ttsync.Var2-ttsync.gvroll,'linewidth',2);
         plot(data.time,ttsync.Var1-ttsync.gvpitch,'linewidth',2);
        plot(data.time,ttsync.Var3-ttsync.gvdrift,'linewidth',2);
        
        legend('HCR roll - GV roll','HCR pitch - GV pitch','HCR drift - GV drift')
        
        grid on
        ylabel('Angles (deg)');
        ylim([-3 3]);
        xlim([data.time(1),data.time(end)]);
        
        title([project,' RF ',num2str(ii)]);
        
        set(gcf,'PaperPositionMode','auto')
        print(f2,[figdir,project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
end