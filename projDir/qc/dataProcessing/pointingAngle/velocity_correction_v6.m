clear all

tilt_error=0.0%-0.0664
rotation_error=0.0%-.85%-0.72%-0.98%
drift_error=0.0
upwind_limit=2;

reuse='n';
% if exist('files')
%     reuse=input('Reuse same files, y or n? ','s');
% end

if reuse=='n'
    %clearvars -except tilt_error rotation_error drift_error upwind_limit
    [files, pathname, filterindex] = uigetfile('/scr/snow2/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/*.nc', 'Pick a file','MultiSelect','on');
    if iscell(files)==1
        numfiles=length(files);
    else
        numfiles=1;
    end
    
    time_long=0;
    index_time=1;
    for i=1:numfiles
        if iscell(files)==1
            nctest2([pathname,files{i}])
        else
            nctest2([pathname,files(i,:)])
        end
        if i==1
            %starttime=time_coverage_start(12:20);
            starttime=datenum(time_coverage_start(12:19),'HH:MM:SS');
            starttime_char=time_coverage_start(12:19);
            title_long=time_coverage_start(1:10);
        elseif i==numfiles
            endtime_char=time_coverage_end(12:20);
        end
        if exist('DBZ','var')==1
            [m n]=size(DBZ);
            VEL_RAW_long(index_time:index_time+m-1,:)=VEL_RAW;
            VEL_long(index_time:index_time+m-1,:)=VEL;
%             VEL_CORR_long(index_time:index_time+m-1,:)=VEL_CORR;
            VEL_CORR_long(index_time:index_time+m-1,:)=VEL;
            DBMHX_long(index_time:index_time+m-1,:)=DBMHX;
            DBMVC_long(index_time:index_time+m-1,:)=DBMVC;
            DBZ_long(index_time:index_time+m-1,:)=DBZ;
            DBZVC_long(index_time:index_time+m-1,:)=DBZVC;
            LDRV_long(index_time:index_time+m-1,:)=LDRV;
            NCP_long(index_time:index_time+m-1,:)=NCP;
            SNRHX_long(index_time:index_time+m-1,:)=SNRHX;
            SNRVC_long(index_time:index_time+m-1,:)=SNRVC;
            WIDTH_long(index_time:index_time+m-1,:)=WIDTH;

            altitude_long(index_time:index_time+m-1)=altitude;
            elevation_long(index_time:index_time+m-1)=elevation;
            azimuth_long(index_time:index_time+m-1)=azimuth;
            time_long(index_time:index_time+m-1)=time+max(time_long);

            tilt_long(index_time:index_time+m-1)=tilt;
            rotation_long(index_time:index_time+m-1)=rotation;
            
            tilt_norm_long(index_time:index_time+m-1)=tilt-nanmean(tilt);

            pitch_long(index_time:index_time+m-1)=pitch;
            heading_long(index_time:index_time+m-1)=heading;
            roll_long(index_time:index_time+m-1)=roll;
            vertical_velocity_long(index_time:index_time+m-1)=vertical_velocity;
            drift_long(index_time:index_time+m-1)=drift;
            eastward_velocity_long(index_time:index_time+m-1)=eastward_velocity;
            northward_velocity_long(index_time:index_time+m-1)=northward_velocity;
            nyquist_velocity_long(index_time:index_time+m-1)=nyquist_velocity;
                        
            latitude_long(index_time:index_time+m-1)=latitude;
            longitude_long(index_time:index_time+m-1)=longitude;


            %time_long(index_time:index_time+m-1)=max(time_long)+time/86400 + datenum(time_coverage_start(12:19),'HH:MM:SS');

            [index_time]=length(time_long)+1;
        end

        i
        numfiles
        %clearvars -except *long range index_time files numfiles nyquist_velocity...
        %                  starttime starttime_char endtime endtime_char title_long tilt_error...
        %                  rotation_error drift_error upwind_limit
    end

    time=time_long;
    ti_test=time/86400 + starttime;
    ti_test=ti_test-floor(ti_test);
    altitude=altitude_long;

    % DBZ_long(DBMVC_long<-103.50)=nan;
    % VEL_long(DBMVC_long<-103.50)=nan;
    % LDRV_long(DBMHX_long<-103.0)=nan;
    % DBZ_long(SNRVC_long<-6.50)=nan;
    % VEL_long(SNRVC_long<-6.50)=nan;
    % LDRV_long(SNRHX_long<0.0)=nan;

    rad=pi/180;
    noise_floor=-107.75;

%     northward_velocity_long=abs(northward_velocity_long)+10;
    ground_speed=sqrt(eastward_velocity_long.^2 + northward_velocity_long.^2);

    [m,n]=size(DBZ_long);

end
% ground_speed=ground_speed+5;

VEL_RAW_long(VEL_RAW_long<-900)=nan;
VEL_corr=VEL_RAW_long;
%VEL_corr(DBMVC<noise_floor)=nan;

% HCR-TEST November 2015
% tau_suba=tilt_long - 0.62;% tilt angle wrt aircraft - for 20141114
tau_suba=tilt_long - tilt_error;% tilt angle wrt aircraft - after first correction
theta_suba=rotation_long - rotation_error; % rotation angle wrt aircraft;
drift_long=drift_long - drift_error;

% IDEAS IV
% tau_suba=tilt_long -0.12;% tilt angle wrt aircraft
% theta_suba=rotation_long - .0; % rotation angle wrt aircraft;

% convert to radians
drift_rad=drift_long*rad;
pitch_rad=pitch_long*rad;
roll_rad=roll_long*rad;
theta_suba_rad=theta_suba*rad;
tau_suba_rad=tau_suba*rad;

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
vr_platform=-ground_speed.*sin(tau_subt)-vertical_velocity_long.*sin(phi);
%vr_platform=-ground_speed.*sin(tau_subt)-vertical_velocity_long.*sind(elevation_long);
%vr_platform=-ground_speed.*sin(tau_subt); %-vertical_velocity_long.*sin(phi);
%VEL_corr(VEL_corr<-900)=nan;
%vr_platform_check=vr_platform;
for i=1:m
for j=1:n
    VEL_corr(i,j)=VEL_corr(i,j)-vr_platform(i);
%     while VEL_corr(i,j)<-nyquist_velocity_long(i)
%         VEL_corr(i,j)=VEL_corr(i,j)+nyquist_velocity_long(i);
%         %vr_platform_check(i)=vr_platform_check(i)+nyquist_velocity_long(i);
%     end
%     while VEL_corr(i,j)>nyquist_velocity_long(i)
%         VEL_corr(i,j)=VEL_corr(i,j)-nyquist_velocity_long(i);
%         %vr_platform_check(i)=vr_platform_check(i)-nyquist_velocity_long(i);
%     end
%     
end
end

% Compute the velocity of the ground
DBZ_tmp=DBZ_long;
DBZ_tmp(:,1:15)=0;
[bla ground_index]=max(DBZ_tmp,[],2);clear DBZ_tmp

% VEL_corr_ground=VEL_corr;
% VEL_corr_ground2=VEL_RAW_long;
for i=1:m
    VEL_ground_raw(i)=VEL_RAW_long(i,ground_index(i));
    VEL_ground(i)=VEL_corr(i,ground_index(i));
    VEL_ground_ucorr(i)=VEL_long(i,ground_index(i));
    VEL_ground_dixon(i)=VEL_CORR_long(i,ground_index(i));
    Width_ground(i)=WIDTH_long(i,ground_index(i));
    DBZ_ground(i)=DBZ_long(i,ground_index(i));
end

figure;plot(ti_test,VEL_ground)
%hold on;plot(ti_test,VEL_ground_ucorr-vertical_velocity_long,'r')
datetick('x',15,'keeplimits');
xlabel('Time, UTC');ylabel('Corrected Ground Vr, m s^-^1')

figure;plot(ti_test,VEL_ground_raw,'r') %-vertical_velocity_long,'r')
hold on;plot(ti_test,vr_platform,'k')
legend('Uncorrected Vr ground','Vr platform');
datetick('x',15,'keeplimits');
xlabel('Time, UTC');ylabel('m s^-^1')

i_level=find(abs(roll_long<5) & elevation_long<-80 & abs(drift_long)<upwind_limit);
Vr_updown_mean=nanmean(VEL_ground(i_level))

if isempty(i_level)==0
    tilt_error_updown=asin(VEL_ground./ground_speed);
    tilt_error_deg_updown=tilt_error_updown*180/pi;
    tilt_error_updown_mean=nanmean(tilt_error_deg_updown(i_level))

    figure;hist(VEL_ground(i_level),50)
    xlabel('Ground Vr (m/s), upwind-downwind')
    title(strcat('Mean ground Vr = ',num2str(Vr_updown_mean),' Mean tilt error = ',num2str(tilt_error_updown_mean)))
end

    i_level2=find(abs(roll_long<5) & elevation_long<-80 & abs(drift_long)>upwind_limit);
    Vr_cross_mean=nanmean(VEL_ground(i_level2))
    
if isempty(i_level2)==0
    ratio=VEL_ground./(ground_speed.*sin(drift_rad));
    ratio(abs(ratio)>1)=nan;
    
    rotation_error_cross=-asind(ratio);
    rotation_error_cross_mean=nanmean(rotation_error_cross(i_level2))
    
    figure;hist(VEL_ground(i_level2),50)
    xlabel('Ground Vr (m/s), cross wind')
    title(strcat('Mean ground Vr = ',num2str(Vr_cross_mean),' Mean rotation error = ',num2str(rotation_error_cross_mean)))   
end
