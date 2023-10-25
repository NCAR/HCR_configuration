function data=readHCRts(fileList,data,startTime,endTime)

vars=fieldnames(data);

file=fileList{1};

% Read HCR time series data
baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset_vc');

fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=fileStartTime+seconds(timeOffset)';

data.range=ncread(file,'range');
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

data.lambda=ncreadatt(file,'/','radar_wavelength_cm')/100;
data.dbz1km_v=ncreadatt(file,'/','cal_base_dbz_1km_vc');
data.dbz1km_h=ncreadatt(file,'/','cal_base_dbz_1km_hc');
data.noise_v=ncreadatt(file,'/','cal_noise_dbm_vc');
data.noise_h=ncreadatt(file,'/','cal_noise_dbm_hc');
data.rx_gain_vc=ncreadatt(file,'/','cal_receiver_gain_db_vc');
data.rx_gain_hc=ncreadatt(file,'/','cal_receiver_gain_db_hc');

for ii=1:length(vars)
    data.(vars{ii})=ncread(file,(vars{ii}));
end

for jj=2:length(fileList)
    file=fileList{jj};

    baseTime=ncread(file,'base_time');
    timeOffset=ncread(file,'time_offset_vc');

    fileStartTime=datetime(1970,1,1)+seconds(baseTime);
    data.time=cat(2,data.time,fileStartTime+seconds(timeOffset)');

    data.elevation=cat(2,data.elevation,ncread(file,'elevation_vc')');
    data.pulse_width=cat(2,data.pulse_width,ncread(file,'pulse_width_vc')');
    data.prt=cat(2,data.prt,ncread(file,'prt_vc')');
    
    for ii=1:length(vars)
        data.(vars{ii})=cat(2,data.(vars{ii}),ncread(file,(vars{ii})));
    end
end

% Trimm times
allVars=fieldnames(data);
timeInds=find(data.time>=startTime & data.time<=endTime);
for ii=1:size(allVars,1)
    if ~strcmp(allVars{ii},'range') & max(size(data.(allVars{ii})))~=1
        if min(size(data.(allVars{ii})))~=1
            data.(allVars{ii})=single(data.(allVars{ii})(:,timeInds));
        else
            data.(allVars{ii})=data.(allVars{ii})(:,timeInds);
        end
    end
end
%data.asl=HCRrange2asl(data.range,data.elevation,data.altitude);
%data.asl=single(data.asl);
end