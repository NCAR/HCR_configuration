% Compare data from before and after velocity correction

clear all
close all

savefig=0;

startTime=datetime(2015,8,5,0,0,0);
endTime=datetime(2015,8,12,0,0,0);

addpath('/h/eol/romatsch/git/private/utils/');

indirHCR='/scr/eldora1/rsfdata/cset/hcr/qc/cfradial/moments/10hz/';
indirCloudSat='/scr/sci/romatsch/data/cloudSat/cset/';

figdir='/h/eol/romatsch/hcrCalib/cloudSat/figs/';
formatOut = 'yyyymmdd_HHMM';

%% Load HCR data

fileListHCR=makeFileList(indirHCR,startTime,endTime,1);

if ~isempty(fileListHCR)
    
    dataHCR.lat=[];
    dataHCR.lon=[];
    dataHCR.time=[];
    %dataHCR.alt=[];
    %dataHCR.elev=[];
    %dataHCR.range=[];
    %dataHCR.dbz=[];    
    
    % Get uncorrected data
    for ii=1:size(fileListHCR,2)
        if mod(ii,10)==0
            disp([num2str(ii),' of ',num2str(size(fileListHCR,2))]);
        end
        infile=fileListHCR{ii};
        
        startTimeIn=ncread(infile,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(infile,'time')';
        timeIn=startTimeFile+seconds(timeRead);
        latIn=ncread(infile,'latitude');
        lonIn=ncread(infile,'longitude');
        
        if length(timeIn)==length(latIn) & length(timeIn)==length(lonIn)
            dataHCR.time=[dataHCR.time,timeIn];
            
            dataHCR.lat=[dataHCR.lat,ncread(infile,'latitude')'];
            dataHCR.lon=[dataHCR.lon,ncread(infile,'longitude')'];
        end
        
%         rangeIn=ncread(infile,'range');
%         rangeMat=repmat(rangeIn,1,length(timeIn));
%         dataHCR.range=[dataHCR.range,rangeMat];
%         
%         dataHCR.elev=[dataHCR.elev,ncread(infile,'elevation')'];
%         
%         dataHCR.alt=[dataHCR.alt,ncread(infile,'altitude')'];
%         
%         dataHCR.dbz=[dataHCR.dbz,ncread(infile,'DBZ')];
        
    end
end

%% Load CloudSat data

[fileListCS,granuleStart]=makeFileList_cloudSat(indirCloudSat,startTime,endTime);

if ~isempty(fileListCS)
    
    dataCS.lat=[];
    dataCS.lon=[];
    dataCS.time=[];
    %dataCS.alt=[];
    %dataCS.elev=[];
    %dataCS.range=[];
    %dataCS.dbz=[];    
    
    % Get uncorrected data
    for ii=1:size(fileListCS,2)
        if mod(ii,10)==0
            disp([num2str(ii),' of ',num2str(size(fileListCS,2))]);
        end
        FILE_NAME=fileListCS{ii};
        
        file_id = hdfsw('open', FILE_NAME, 'rdonly');
        SWATH_NAME = '2B-GEOPROF';
        swath_id = hdfsw('attach', file_id, SWATH_NAME);
        
        if ii==1
             % Read attributes.
             [long_name, status] = hdfsw('readattr', swath_id,'Radar_Reflectivity.long_name');
             [units, status] = hdfsw('readattr', swath_id,'Radar_Reflectivity.units');
             
             [units_h, status] = hdfsw('readattr', swath_id,'Height.units');
             
             [units_t, status] = hdfsw('readattr', swath_id,'Profile_time.units');
             [long_name_t, status] = hdfsw('readattr', swath_id,'Profile_time.long_name');
         end
        
        % Read data.
        
        DATAFIELD_NAME = 'Radar_Reflectivity';
        [data, status] = hdfsw('readfield', swath_id, DATAFIELD_NAME, [],[],[]);
        
        % Read lat/lon/height/time data.
        [lon, status] = hdfsw('readfield', swath_id, 'Longitude', [], [], []);
        [lat, status] = hdfsw('readfield', swath_id, 'Latitude', [], [], []);
        [height, status] = hdfsw('readfield', swath_id, 'Height', [], [], []);
        [timeIn, status] = hdfsw('readfield', swath_id, 'Profile_time', [], [], []);
        
        % Make type double for plotting.
         if length(timeIn)==length(lat) & length(timeIn)==length(lon)
            dataCS.lat=[dataCS.lat;double(lat)];
            dataCS.lon=[dataCS.lon;double(lon)];
            timeIn=double(timeIn);
            data=double(data);
            
            [scale_factor, status] = hdfsw('readattr', swath_id,'Radar_Reflectivity.factor');
            scale_factor = double(scale_factor);
            
            [valid_range, status] = hdfsw('readattr', swath_id,'Radar_Reflectivity.valid_range');
            
            hdfsw('detach', swath_id);
            hdfsw('close', file_id);
            
            %Convert time
            dataCS.time=[dataCS.time;granuleStart(ii)+seconds(timeIn)];
            
            % Process valid_range. Fill value and missing value will be handled by this
            % since they are outside of range values.
            data((data < valid_range(1)) | (data > valid_range(2))) = NaN;
            
            % Apply scale factor according to [1].
            data = data / scale_factor;
        end
            
        
    end
end

%% Find overlap
roundHCRtime=dataHCR.time;
roundHCRtime.Second=round(roundHCRtime.Second);
roundCStime=dataCS.time;
roundCStime.Second=round(roundCStime.Second);

HCRdata=table(datenum(roundHCRtime)',round(dataHCR.lon,4)',round(dataHCR.lat,4)');
CSdata=table(datenum(roundCStime),round(dataCS.lon,4),round(dataCS.lat,4));

[C,iHCR,iCS] = intersect(HCRdata,CSdata,'rows');