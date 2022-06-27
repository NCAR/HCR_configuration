function data=readHCRts(data,file)

vars=fieldnames(data);

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
data.noise_v=ncreadatt(file,'/','cal_noise_dbm_vc');
data.rx_gain_v=ncreadatt(file,'/','cal_receiver_gain_db_vc');

for ii=1:length(vars)
    data.(vars{ii})=ncread(file,(vars{ii}));
end

end