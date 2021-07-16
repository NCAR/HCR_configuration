% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

upwind_limit=2;
crosswind_limit=2;

% If attitude correction was sometimes turned off we need to filter that
% data
filterAttitudeCorrOff=1;

project='spicule'; % socrates, cset, aristo, otrec
quality='field'; % field, qc0, qc1, qc2
qcVersion='v0.1';
freq='10hz';

figdir=['/h/eol/romatsch/hcrCalib/pointAng/',project];
formatOut = 'yyyymmdd_HHMM';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freq);

lineCols=lines;

gvfiles=dir('/scr/sleet2/rsfdata/projects/spicule/GV/SPICULErf*.nc');

%% Run processing

tiltSum=0;
tiltPix=0;
tiltErrorFlight=[];

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
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
    
    tiltSum=tiltSum+nansum(abs_tilt_error);
    tiltPix=tiltPix+length(find(~isnan(abs_tilt_error)==1));
    
    tiltErrorFlight=[tiltErrorFlight nansum(abs_tilt_error)/length(find(~isnan(abs_tilt_error)==1))];
    
    disp(['Tilt error for flight ',num2str(ii),' is ',num2str(tiltErrorFlight(end)),' deg.']);
    
    %% Find tilt angle events
    nanInds=find(isnan(abs_tilt_error));
    timeTemp=data.time';
    tiltTemp=abs_tilt_error;
    
    timeTemp(nanInds)=[];
    tiltTemp(nanInds)=[];
    
    timeDiff=etime(datevec(timeTemp(2:end)),datevec(timeTemp(1:end-1)));
    endGapInd=find(timeDiff > 60);
    startGapInd=endGapInd+1;
    endGapInd=cat(1,endGapInd,length(timeTemp));
    startGapInd=cat(1,1,startGapInd);
    
    tiltEvent=[];
    timeEvent=[];
    
    for jj=1:length(endGapInd)
        tiltDataE=tiltTemp(startGapInd(jj):endGapInd(jj));
        timeDataE=timeTemp(startGapInd(jj):endGapInd(jj));
        if length(tiltDataE)>600
            tiltEvent=[tiltEvent,nanmedian(tiltDataE)];
            timeEvent=[timeEvent,timeDataE(round(length(tiltDataE)/2))];
        end
    end
    
    % If there are not tilt events then create one with 0 error
    if isempty(tiltEvent)
        tiltEvent=0;
        timeEvent=startTime;
    end
    
    %% Find rotation error events
    rotOnes=ones(size(data.rotation));
    
    % Sort out bad data
    i_level2=find(abs(data.roll)>5 | abs(data.drift)<crosswind_limit | data.pitch>3.5 | data.pitch<1 ...
        | abs(data.rotation-180)>5 | abs(VEL_ground)>7 | abs(vel_ground_mean)>5 | data.elevation>0);
    rotOnes(i_level2)=nan;
    
    % Exclude data where attitude correction was turned off
    if filterAttitudeCorrOff
        rotOnes(attInd)=nan;
    end
    
    nanIndsR=find(isnan(rotOnes));
    timeTempR=data.time';
    indsTemp=1:1:length(data.time);
    
    timeTempR(nanIndsR)=[];
    indsTemp(nanIndsR)=[];
    
    timeDiffR=etime(datevec(timeTempR(2:end)),datevec(timeTempR(1:end-1)));
    endGapIndR=find(timeDiffR > 60);
    startGapIndR=endGapIndR+1;
    endGapIndR=cat(1,endGapIndR,length(timeTempR));
    startGapIndR=cat(1,1,startGapIndR);
    
    rotEvent=[];
    timeEventR=[];
    
    goodRot=nan(size(data.rotation));
   
    % Go through events and calculate rotatio error
    for jj=1:length(endGapIndR)-1
        timeDataE=timeTempR(startGapIndR(jj):endGapIndR(jj));
        indsE=indsTemp(startGapIndR(jj):endGapIndR(jj));
        if length(timeDataE)>600
            timeEventR=[timeEventR,timeDataE(round(length(timeDataE)/2))];
            % Find closest tilt
            timeDiff2tilt=abs(etime(datevec(timeEventR(end)),datevec(timeEvent)));
            minDiff=find(timeDiff2tilt==min(timeDiff2tilt));
            tiltRotEvent=tiltEvent(minDiff);
            [VEL_groundE ground_speedE]=calcVelCorr(data.drift(indsE),data.pitch(indsE),...
                data.roll(indsE),data.rotation(indsE),data.tilt(indsE)-tiltRotEvent,...
                data.eastward_velocity(indsE),data.northward_velocity(indsE),data.vertical_velocity(indsE),...
                data.VEL_RAW(:,indsE),data.DBZ(:,indsE));
            
            % Rotation error
            ground_speed_cross=ground_speedE.*sin(deg2rad(data.drift(indsE))); % Cross track component
            ratio=VEL_groundE./ground_speed_cross;
            ratio(abs(ratio)>1)=nan;
            rotation_error_cross=-asind(ratio);
   
            goodRot(indsE)=rotation_error_cross;
            rotEvent=[rotEvent,nanmedian(rotation_error_cross)];
        end
    end
    
    
    %%
    close all
    
    dataInds=find(~isnan(vel_ground_mean));
    xlims=[data.time(min(dataInds)),data.time(max(dataInds))];
    
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 2000 1300]);
    
    s1=subplot(4,1,1);
    hold on
    yyaxis right
    l1=plot(data.time,data.roll,'-','linewidth',2,'color',lineCols(1,:));
    l2=plot(data.time,data.pitch,'-','linewidth',2,'color',lineCols(2,:));
    l3=plot(data.time,data.drift,'-','linewidth',2,'color',lineCols(3,:));
           
    ylabel('Angles (deg)');
    
    ylim([-14 14]);
    yticks(-20:4:20);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=[0 0 0];
    grid on
    
    yyaxis left
    l4=plot(data.time,movmean(data.elevation+90,100),'-','linewidth',2,'color',lineCols(4,:));
    ylim([-0.175 0.175]);
    yticks(-5:0.05:5);
    ylabel('Elev+90 (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(4,:);
    
    xlim(xlims);
       
    title([project,' RF ',num2str(ii)]);
    
    plot(xlims,[0 0],'-k');
    s1.SortMethod = 'childorder';
    legend([l1 l2 l3 l4],{'Elev+90','Roll','Pitch','Drift'});
        
    s2=subplot(4,1,2);
    hold on
    l1=plot(data.time,data.vertical_velocity,'linewidth',1);
    l2=plot(data.time,VEL_ground,'linewidth',1);
    l3=plot(data.time,vel_ground_mean,'linewidth',2);
    
    ylabel('Velocity (m/s)');
    
    xlim(xlims);
    ylim([-2 2]);
    yticks(-10:0.5:10);
    
    t1=text(xlims(1)+minutes(5),-1.5,['Mean ground velocity error: ',num2str(mean(VEL_ground,'omitnan')),' m s^{-1}'],'fontsize',11);
    t1.BackgroundColor=[1,1,1];
    grid on
    plot(xlims,[0 0],'-k');
    s2.SortMethod = 'childorder';
    legend([l1 l2 l3],{'Vel vert','Vel ground','Vel ground mean'});
    
    s3=subplot(4,1,3);
    hold on
      
    yyaxis left
    %plot(data.time,rotation_error_cross_all,'linewidth',1,'color',lineCols(1,:));
    l1=plot(data.time,movmean(rotation_error_cross_all,100),'-','linewidth',2,'color',lineCols(6,:));
    ylim([-2 2]);
    yticks(-5:0.5:5);
    ylabel('Rotation angle (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(6,:);
    
    yyaxis right
    %plot(data.time,tilt_error_updown,'linewidth',1,'color',lineCols(4,:));
    l2=plot(data.time,movmean(tilt_error_updown,100),'-','linewidth',2,'color',lineCols(2,:));
    ylim([-0.4 0.4]);
    yticks(-5:0.1:5);
    ylabel('Tilt angle (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(2,:);
    
    %legend('Tilt error all','Tilt error good','Rotation error all','Rotation error good');
        
    xlim(xlims);
        
    grid on
    
    plot(xlims,[0 0],'-k');
    s3.SortMethod = 'childorder';
    legend([l1,l2],{'Max rotation error','Max tilt error'});
    
    s4=subplot(4,1,4);
    hold on
    
    yyaxis left
    plot(data.time,goodRot,'linewidth',1,'color',lineCols(1,:));
    l1=plot(data.time,movmean(goodRot,100),'-','linewidth',2,'color',lineCols(6,:));
    l3=scatter(timeEventR,rotEvent,70,'o','filled','markerfacecolor',lineCols(5,:),'markeredgecolor','w');
    ylim([-2 2]);
    yticks(-5:0.5:5);
    ylabel('Rotation angle (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(6,:);
    
    yyaxis right
    plot(data.time,abs_tilt_error,'linewidth',1,'color',lineCols(4,:));
    l2=plot(data.time,movmean(abs_tilt_error,100),'-','linewidth',2,'color',lineCols(2,:));    
    l4=scatter(timeEvent,tiltEvent,70,'o','filled','markerfacecolor',lineCols(3,:),'markeredgecolor','w');
    ylim([-0.4 0.4]);
    yticks(-5:0.1:5);
    ylabel('Tilt angle (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(2,:);
            
    xlim(xlims);
        
    t1=text(xlims(1)+minutes(5),-0.3,['Mean tilt error: ',num2str(mean(abs_tilt_error,'omitnan')),...
        ' deg, mean rotation error: ',num2str(mean(goodRot,'omitnan')),' deg'],'fontsize',11);
    t1.BackgroundColor=[1,1,1];
    grid on
    plot(xlims,[0 0],'-k');
    s4.SortMethod = 'childorder';
    legend([l1,l2,l3,l4],{'Rotation error','Tilt error','Rot error event','Tilt error event'});
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'/wholeFlights/',project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    %% Drift from GV
  
    gvfile=[gvfiles(ii).folder,'/',gvfiles(ii).name];
    
    gvtimeIn=ncread(gvfile,'Time');
    refTimegv=datetime(year(data.time(1)),month(data.time(1)),day(data.time(1)),0,0,0);
    gvdrift=ncread(gvfile,'DRFTA');
    gvroll=ncread(gvfile,'ROLL');
    gvpitch=ncread(gvfile,'PITCH');
    
    if min(size(gvdrift))~=1
        gvdrift=gvdrift(1,:)';
        gvroll=gvroll(1,:)';
        gvpitch=gvpitch(1,:)';
    end
    
    gvtime=refTimegv+seconds(gvtimeIn);
    
    hcrTT=timetable(data.time',data.pitch',data.roll',data.drift');
    gvTT=timetable(gvtime,gvpitch,gvroll,gvdrift);
    
    ttsync=synchronize(hcrTT,gvTT,'first','linear');
    
    %% Plot GV data
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 2000 1100]);
    
    s1=subplot(3,1,1);
    
    hold on
    
    yyaxis left
    l1=plot(data.time,data.roll,'-','linewidth',2,'color',lineCols(1,:));
    l2=plot(gvtime,gvroll,'-','linewidth',2,'color',lineCols(6,:));
    ylabel('Roll (deg)');
    ylim([-20 20]);
    yticks(-20:5:20);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(6,:);
    
    yyaxis right
    l3=plot(data.time,data.pitch,'-','linewidth',2,'color',lineCols(4,:));
    l4=plot(gvtime,gvpitch,'-','linewidth',2,'color',lineCols(2,:));
    l5=plot(data.time,data.drift,'-','linewidth',2,'color',lineCols(5,:));
    l6=plot(gvtime,gvdrift,'-','linewidth',2,'color',lineCols(3,:));
    ylabel('Pitch and drift (deg)');
    ylim([-8 8]);
    yticks(-8:2:8);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(2,:);
    
    plot(xlims,[0 0],'-k');
    s1.SortMethod = 'childorder';
    
    legend([l1 l2 l3 l4 l5 l6],{'HCR roll','GV roll','HCR pitch','GV pitch','HCR drift','GV drift'});
    
    grid on
    xlim(xlims);
    
    title([project,' RF ',num2str(ii)]);
    
    s2=subplot(3,1,2);
    hold on
    l1=plot(data.time,ttsync.Var2-ttsync.gvroll,'linewidth',2);
    l2=plot(data.time,ttsync.Var1-ttsync.gvpitch,'linewidth',2);
    l3=plot(data.time,ttsync.Var3-ttsync.gvdrift,'linewidth',2);
    
    plot(xlims,[0 0],'-k');
    s2.SortMethod = 'childorder';
    
    legend([l1 l2 l3],{'HCR roll - GV roll','HCR pitch - GV pitch','HCR drift - GV drift'})
    
    grid on
    ylabel('Angles (deg)');
    ylim([-3 3]);
    xlim(xlims);
    
    s3=subplot(3,1,3);
    hold on
    
    yyaxis left
    l1=plot(data.time,ttsync.Var1-ttsync.gvpitch,'linewidth',2,'color',lineCols(2,:));
    ylim([-1 1]);
    yticks(-1:0.2:1);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(2,:);
    ylabel('Pitch (deg)');
    
    yyaxis right
    l2=plot(data.time,movmean(abs_tilt_error,100),'-','linewidth',2,'color',lineCols(4,:));
    ylim([-0.2 0.2]);
    yticks(-1:0.04:1);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(4,:);
    
    plot(xlims,[0 0],'-k');
    s3.SortMethod = 'childorder';
    
    legend([l1 l2],{'HCR pitch - GV pitch','Tilt error'})
    
    grid on
    ylabel('Tilt (deg)');
    xlim(xlims);
    
    set(gcf,'PaperPositionMode','auto')
    print(f2,[figdir,'/insCompare/',project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
 
    figure
    scatter(movmean(abs_tilt_error,100),ttsync.Var1-ttsync.gvpitch)
    title([project,' RF ',num2str(ii)]);
    xlabel('Tilt error (deg)');
    ylabel('HCR pitch - GV pitch (deg)');
    print([figdir,'/insCompare/',project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_scatter'],'-dpng','-r0')
end
%% 
close all
f3=figure('DefaultAxesFontSize',12);
set(f3,'renderer','painters');
set(f3,'Position',[200 500 1000 600]);

plot(tiltErrorFlight,'linewidth',2)
xlim([1 length(tiltErrorFlight)])
ylim([-0.15 0.15])
xticks(1:length(tiltErrorFlight))
xlabel('Flight')
ylabel('Mean tilt error (deg)')
grid on

title([project,' mean tilt errors'])

set(gcf,'PaperPositionMode','auto')
print(f3,[figdir,'/wholeFlights/',project,'_meanTiltErrors'],'-dpng','-r0')

disp(['Mean tilt error is ',num2str(tiltSum/tiltPix),' deg']);