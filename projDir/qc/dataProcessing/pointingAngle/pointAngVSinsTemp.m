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
quality='qc0'; % field, qc1, qc2
freq='10hz';

figdir=['/h/eol/romatsch/hcrCalib/pointAng/',project,'/tiltVSinsTemp/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freq);

lineCols=lines;

gvfiles=dir('/scr/snow2/rsfdata/projects/otrec/GV/OTRECrf*.nc');

%% Run processing
tableAll=[];

% Go through flights
for ii=4:9
    disp(['Flight ',num2str(ii),' of 9.']);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    data.roll=[];
    data.pitch=[];
    data.rotation=[];
    data.tilt=[];
    data.drift=[];
    data.northward_velocity=[];
    data.eastward_velocity=[];
    data.vertical_velocity=[];
    data.altitude_agl=[];
    
    data.VEL_RAW=[];
    data.DBZ=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load insTemp
    insT=[];
    timeTemp=[];
    
    for kk=1:size(fileList,2)
        tempFile=fileList{kk};
        insT=[insT varFromCfRadialString(tempFile,'InsTemp')];
       
        startTimeIn=ncread(tempFile,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeTemp=[timeTemp startTimeFile];
    end
    
    % Load HCR data
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
    
    % Remove zenith pointing
    zenith=find(data.elevation>0);
    
    for kk=1:length(dataVars)
        dataTemp=data.(dataVars{kk});
        dataTemp(:,zenith)=nan;
        data.(dataVars{kk})=dataTemp;
    end
    
    [VEL_ground ground_speed]=calcVelCorr(data.drift,data.pitch,data.roll,data.rotation,data.tilt,...
        data.eastward_velocity,data.northward_velocity,data.vertical_velocity,...
        data.VEL_RAW,data.DBZ);
    
    %% Calculate pointing angle errors
    
    % Tilt error for all data
    tilt_error_updown=asind(VEL_ground./ground_speed);
    
    % Rotation error for all data
    ground_speed_cross_all=ground_speed.*sin(deg2rad(data.drift)); % Cross track component
    ratio_all=VEL_ground./ground_speed_cross_all;
    ratio_all(abs(ratio_all)>1)=nan;
    rotation_error_cross_all=-asind(ratio_all);
    
    tilt_error_good=tilt_error_updown;
    
    vel_ground_mean=movmean(VEL_ground,200);
    
    % Sort out bad data
    i_level=find(abs(data.roll)>5 | abs(data.drift)>upwind_limit | data.pitch>3.5 | data.pitch<1 ...
        | abs(data.rotation-180)>5 | abs(VEL_ground)>7 | abs(vel_ground_mean)>5);
    tilt_error_good(i_level)=nan;
    
    % Exclude data where attitude correction was turned off
    if filterAttitudeCorrOff
        attInd=find(abs(data.tilt)<0.2);
        tilt_error_good(attInd)=nan;
    end
   
    abs_tilt_error=tilt_error_good;
    posIndsTilt=find(data.tilt>0);
    abs_tilt_error(posIndsTilt)=-abs_tilt_error(posIndsTilt);
           
    Ttilt=timetable(data.time',abs_tilt_error');
    Ttemp=timetable(timeTemp',insT');
    
    TT = synchronize(Ttemp,Ttilt,'first','mean');
    
    tableAll=cat(1,tableAll,TT);
    
    close all
    figure
    scatter(TT.Var1_Ttilt,TT.Var1_Ttemp)
    title([project,' RF ',num2str(ii)]);
    xlabel('Tilt error (deg)');
    ylabel('INS temperature (C)');
    xlim([-0.3 0.1]);
    print([figdir,project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_tiltVSinsTemp'],'-dpng','-r0')
end
%% 
edges=-40:2:30;

meanTilts=[];
meanTemps=[];

for ii=1:length(edges)-1
    tempInds=find(tableAll.Var1_Ttemp>edges(ii) & tableAll.Var1_Ttemp<=edges(ii+1));
    if ~isempty(tempInds)
        tiltsInds=tableAll.Var1_Ttilt(tempInds);
        meanTilts=[meanTilts nanmean(tiltsInds)];
        meanTemps=[meanTemps edges(ii)+1];
    end
end
    

figure
hold on
scatter(tableAll.Var1_Ttilt,tableAll.Var1_Ttemp)
plot(meanTilts,meanTemps,'linewidth',2)
title('OTREC RF04 to RF09');
xlabel('Tilt error (deg)');
ylabel('INS temperature (C)');
xlim([-0.3 0.1]);
legend('All data','2 deg mean','location','northwest')
print([figdir,project,'_RF04to09_tiltVSinsTemp'],'-dpng','-r0')