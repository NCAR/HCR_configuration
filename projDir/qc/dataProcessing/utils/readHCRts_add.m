function data=readHCRts_add(file,data,vars)

% Read HCR time series data

baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset_vc');

fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=cat(2,data.time,fileStartTime+seconds(timeOffset)');

data.range=ncread(file,'range');
data.latitude=cat(2,data.latitude,ncread(file,'latitude')');
data.longitude=cat(2,data.longitude,ncread(file,'longitude')');
data.altitude=cat(2,data.altitude,ncread(file,'altitude')');
data.elevation=cat(2,data.elevation,ncread(file,'elevation_vc')');
data.pulse_width=cat(2,data.pulse_width,ncread(file,'pulse_width_vc')');
data.prt=cat(2,data.prt,ncread(file,'prt_vc')');
data.lambda=ncreadatt(file,'/','radar_wavelength_cm')/100;
data.dbz1km_v=ncreadatt(file,'/','cal_base_dbz_1km_vc');
data.dbz1km_h=ncreadatt(file,'/','cal_base_dbz_1km_hc');
data.noise_v=ncreadatt(file,'/','cal_noise_dbm_vc');
data.noise_h=ncreadatt(file,'/','cal_noise_dbm_hc');
data.rx_gain_v=ncreadatt(file,'/','cal_receiver_gain_db_vc');
data.rx_gain_h=ncreadatt(file,'/','cal_receiver_gain_db_hc');
data.beamwidth_v=ncreadatt(file,'/','radar_beamwidth_deg_v');
data.beamwidth_h=ncreadatt(file,'/','radar_beamwidth_deg_h');

data.lambda=ncreadatt(file,'/','radar_wavelength_cm')/100;
data.dbz1km_v=ncreadatt(file,'/','cal_base_dbz_1km_vc');
data.dbz1km_h=ncreadatt(file,'/','cal_base_dbz_1km_hc');
data.noise_v=ncreadatt(file,'/','cal_noise_dbm_vc');
data.noise_h=ncreadatt(file,'/','cal_noise_dbm_hc');
data.rx_gain_vc=ncreadatt(file,'/','cal_receiver_gain_db_vc');
data.rx_gain_hc=ncreadatt(file,'/','cal_receiver_gain_db_hc');

for ii=1:length(vars)
    thisVar=ncread(file,(vars{ii}));
    if ~strcmp(vars{ii},'range') & size(thisVar,2)==1
        thisVar=thisVar';
    end
    data.(vars{ii})=cat(2,data.(vars{ii}),thisVar);
end
end