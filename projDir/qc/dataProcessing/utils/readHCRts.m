function data=readHCRts(fileList,data,startTime,endTime)

vars=fieldnames(data);

file=fileList{1};

% Read HCR time series data
baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset');

fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=fileStartTime+seconds(timeOffset)';

data.range=ncread(file,'range');
data.elevation=ncread(file,'elevation')';
data.pulse_width=ncread(file,'pulse_width')';
data.prt=ncread(file,'prt')';
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
    timeOffset=ncread(file,'time_offset');

    fileStartTime=datetime(1970,1,1)+seconds(baseTime);
    data.time=cat(2,data.time,fileStartTime+seconds(timeOffset)');

    data.elevation=cat(2,data.elevation,ncread(file,'elevation')');
    data.pulse_width=cat(2,data.pulse_width,ncread(file,'pulse_width')');
    data.prt=cat(2,data.prt,ncread(file,'prt')');
    
    for ii=1:length(vars)
        data.(vars{ii})=cat(2,data.(vars{ii}),ncread(file,(vars{ii})));
    end
end
end