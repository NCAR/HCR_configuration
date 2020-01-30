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

if strcmp(project,'otrec')
    infile='/h/eol/romatsch/hcrCalib/pointAng/flights_otrec.txt';
else
    infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];
end

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freq);

lineCols=lines;

gvfile='/scr/snow2/rsfdata/projects/otrec/GV/OTRECtf02.nc';

%% Run processing

% Go through flights
for ii=4:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
                    
        % Desired variables. The variable name comes after the . and must be spelled exactly
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
        
        %% Calculate Vel_corr
        
        % Convert to radian
        drift_rad=deg2rad(data.drift);
        pitch_rad=deg2rad(data.pitch);
        roll_rad=deg2rad(data.roll);
        theta_suba_rad=deg2rad(data.rotation);
        tau_suba_rad=deg2rad(data.tilt);
        
        % Compute y_t following equation 9 Lee et al (1994)
        y_subt=-cos(theta_suba_rad+roll_rad).*cos(drift_rad).*cos(tau_suba_rad).*sin(pitch_rad)...
            +sin(drift_rad).*sin(theta_suba_rad+roll_rad).*cos(tau_suba_rad)...
            +cos(pitch_rad).*cos(drift_rad).*sin(tau_suba_rad);
        
        % Compute z following equation 9 Lee et al (1994)
        z=cos(pitch_rad).*cos(tau_suba_rad).*cos(theta_suba_rad+roll_rad)+sin(pitch_rad).*sin(tau_suba_rad);
        
        % compute tau_t following equation 11 Lee et al (1994)
        tau_subt=asin(y_subt);
        
        % Compute phi following equation 17 Lee et al (1994)
        phi=asin(z);
        
        % Compute platform motion based on Eq 27 from Lee et al (1994)
        ground_speed=sqrt(data.eastward_velocity.^2 + data.northward_velocity.^2);
        vr_platform=-ground_speed.*sin(tau_subt)-data.vertical_velocity.*sin(phi);
        
        VEL_corr=data.VEL_RAW-vr_platform;
        
        % Compute the velocity of the ground
        DBZ_tmp=data.DBZ;
        DBZ_tmp(1:15,:)=0;
        [bla ground_index]=max(DBZ_tmp,[],1);
        
        ground_index(ground_index==1)=nan;
        
        VEL_ground_raw=nan(1,size(data.VEL_RAW,2));
        VEL_ground=nan(1,size(data.VEL_RAW,2));
        
        for jj=1:size(data.VEL_RAW,2)
            if ~isnan(ground_index(jj))
                if data.DBZ(ground_index(jj),jj)>20
                    VEL_ground_raw(jj)=data.VEL_RAW(ground_index(jj),jj);
                    VEL_ground(jj)=VEL_corr(ground_index(jj),jj);
                end
            end
        end
        
        %% Calculate pointing angle errors
        
        % Tilt error
        tilt_error_updown=asind(VEL_ground./ground_speed);
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
        
        negTilt=find(data.tilt<0);
        tilt_error_abs=tilt_error_good;
        tilt_error_abs(negTilt)=-tilt_error_good(negTilt);
        
                
        %%
        close all
        
        f1=figure('DefaultAxesFontSize',12);
        set(f1,'Position',[200 500 2000 1300]);
        
        subplot(3,1,1)
        hold on
        plot(data.time,data.roll,'linewidth',2);
        plot(data.time,data.pitch,'linewidth',2);
        plot(data.time,data.drift,'linewidth',2);
        plot(data.time,data.rotation-180,'linewidth',2);
        plot(data.time,data.tilt,'linewidth',2);
        
        legend('Roll','Pitch','Drift','Rotation-180','Tilt');
        
        ylabel('Angles (deg)');
        
        xlim([data.time(1) data.time(end)]);
        ylim([-14 14]);
        yticks(-20:2:20);
        
        title([project,' flight ',num2str(ii)]);
        
        grid on
        
        subplot(3,1,2)
        hold on
        plot(data.time,VEL_ground_raw,'linewidth',2);
        plot(data.time,VEL_ground,'linewidth',2);
        plot(data.time,vel_ground_mean,'linewidth',2);
        plot(data.time,vr_platform,'linewidth',2);
        plot(data.time,data.vertical_velocity,'linewidth',2);
                
        legend('Vel ground raw','Vel ground','Vel ground mean','Vel plattform','Vel vert');
        
        ylabel('Velocity (m/s)');
        
        xlim([data.time(1) data.time(end)]);
        ylim([-20 20]);
        
        grid on
        
        subplot(3,1,3)
        hold on
        %plot(data.time,tilt_error_updown,'linewidth',1,'color',[0.7 0.7 0.7]);
        plot(data.time,tilt_error_good,'linewidth',2);
                        
        %legend('Tilt error all','Tilt error good','Rotation error all','Rotation error good');
        legend('Tilt error');
        
        ylabel('Angles (deg)');
        
        xlim([data.time(1) data.time(end)]);
        ylim([-1 1]);
        yticks(-5:0.2:5);
        
        grid on
        
        set(gcf,'PaperPositionMode','auto')
        %        print(f1,[figdir,project,'_Flight',num2str(ii),'_tilt_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
           
        %% Drift from GV
        close all
        
        gvtimeIn=ncread(gvfile,'Time');
        refTimegv=datetime(2019,7,29,0,0,0);
        dr1=ncread(gvfile,'DRFTA_IRS2');
        dr2=ncread(gvfile,'DRFTA_IRS3');
        
        gvtime=refTimegv+seconds(gvtimeIn);
        
        f2=figure('DefaultAxesFontSize',12);
        set(f2,'Position',[200 500 2000 800]);
        
        hold on
        plot(data.time,data.drift,'linewidth',2);
        plot(gvtime,dr1,'linewidth',2);
        plot(gvtime,dr2,'linewidth',2);
        
        legend('HCR INS','GV INS2','GV INS3')
        
        grid on
        ylabel('Drift (deg)');
        xlim([datetime(2019,7,29,15,0,0),datetime(2019,7,29,17,30,0)]);
end