function [data,fileFirst,fileLast]=readHCRts(fileList,data,startTime,endTime,cutTime)

fileFirst=nan(length(fileList),1);
fileLast=nan(length(fileList),1);

if ~isempty(data)
    vars=fieldnames(data);
else
    vars=[];
end

file=fileList{1};

% Read HCR time series data
baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset_vc');

fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=fileStartTime+seconds(timeOffset)';

data.range=ncread(file,'range');
data.latitude=ncread(file,'latitude')';
data.longitude=ncread(file,'longitude')';
data.altitude=ncread(file,'altitude')';
data.elevation=ncread(file,'elevation_vc')';
data.pulse_width=ncread(file,'pulse_width_vc')';
data.prt=ncread(file,'prt_vc')';
data.lambda=ncreadatt(file,'/','radar_wavelength_cm')/100;
data.dbz1km_v=ncreadatt(file,'/','cal_base_dbz_1km_vc');
data.dbz1km_h=ncreadatt(file,'/','cal_base_dbz_1km_hc');
data.noise_v=ncreadatt(file,'/','cal_noise_dbm_vc');
data.noise_h=ncreadatt(file,'/','cal_noise_dbm_hc');
data.rx_gain_v=ncreadatt(file,'/','cal_receiver_gain_db_vc');
data.rx_gain_h=ncreadatt(file,'/','cal_receiver_gain_db_hc');
data.beamwidth_v=ncreadatt(file,'/','radar_beamwidth_deg_v');
data.beamwidth_h=ncreadatt(file,'/','radar_beamwidth_deg_h');

for ii=1:length(vars)
    data.(vars{ii})=ncread(file,(vars{ii}));
    if ~strcmp(vars{ii},'range') & size(data.(vars{ii}),2)==1
        data.(vars{ii})=data.(vars{ii})';
    end
end

fileFirst(1)=1;
fileLast(1)=length(data.time);

for jj=2:length(fileList)
    file=fileList{jj};

    baseTime=ncread(file,'base_time');
    timeOffset=ncread(file,'time_offset_vc');

    fileStartTime=datetime(1970,1,1)+seconds(baseTime);
    data.time=cat(2,data.time,fileStartTime+seconds(timeOffset)');

    data.elevation=cat(2,data.elevation,ncread(file,'elevation_vc')');
    data.latitude=cat(2,data.latitude,ncread(file,'latitude')');
    data.longitude=cat(2,data.longitude,ncread(file,'longitude')');
    data.altitude=cat(2,data.altitude,ncread(file,'altitude')');
    data.pulse_width=cat(2,data.pulse_width,ncread(file,'pulse_width_vc')');
    data.prt=cat(2,data.prt,ncread(file,'prt_vc')');

    for ii=1:length(vars)
        readThis=ncread(file,(vars{ii}));
        if ~strcmp(vars{ii},'range') & size(readThis,2)==1
            readThis=readThis';
        end
        data.(vars{ii})=cat(2,data.(vars{ii}),single(readThis));
    end
    fileFirst(jj)=fileLast(jj-1)+1;
    fileLast(jj)=length(data.time);
end

% Trimm times
if cutTime
    allVars=fieldnames(data);
    noTimeInds=find(data.time<startTime | data.time>endTime);
    for ii=1:size(allVars,1)
        if ~strcmp(allVars{ii},'range') & ~strcmp(allVars{ii},'time') & max(size(data.(allVars{ii})))~=1
            data.(allVars{ii})(:,noTimeInds)=[];
        end
    end
    data.time(:,noTimeInds)=[];
end
end