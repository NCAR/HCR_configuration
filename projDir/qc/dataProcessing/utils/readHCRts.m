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

for ii=1:length(vars)
    data.(vars{ii})=ncread(file,(vars{ii}));
end

end