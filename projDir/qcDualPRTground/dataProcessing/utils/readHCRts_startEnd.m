function [data,fileFirst,fileLast]=readHCRts_startEnd(fileList,data,readStart,sampleNum)

if ~isempty(data)
    vars=fieldnames(data);
else
    vars=[];
end

% Handle more than one file
if length(fileList)>1
    readStart=cat(1,readStart,ones(length(fileList)-1,1));
    sampleNumNew=nan(length(fileList),1);
    tinfo=ncinfo(fileList{1},'time_offset_vc');
    dim1=tinfo.Dimensions.Length;
    sampleNumNew(1)=dim1-readStart(1);
    for hh=2:length(fileList)-1
        tinfo=ncinfo(fileList{hh},'time_offset_vc');
        dim1=tinfo.Dimensions.Length;
        sampleNumNew(hh)=dim1;
    end
    allSamples=sum(sampleNumNew,'omitmissing');
    sampleNumNew(end)=sampleNum-allSamples;
    sampleNum=sampleNumNew;
end

file=fileList{1};

% Read HCR time series data
baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset_vc',readStart(1),sampleNum(1));

fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=fileStartTime+seconds(timeOffset)';

data.range=ncread(file,'range');
data.latitude=ncread(file,'latitude',readStart(1),sampleNum(1))';
data.longitude=ncread(file,'longitude',readStart(1),sampleNum(1))';
data.altitude=ncread(file,'altitude',readStart(1),sampleNum(1))';
data.elevation=ncread(file,'elevation_vc',readStart(1),sampleNum(1))';
data.pulse_width=ncread(file,'pulse_width_vc',readStart(1),sampleNum(1))';
data.prt=ncread(file,'prt_vc',readStart(1),sampleNum(1))';
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

dim2D=zeros(size(vars));

for ii=1:length(vars)
    vinfo=ncinfo(file,vars{ii});
    if max(size(vinfo.Dimensions))==2
        dim2D(ii)=1;
        data.(vars{ii})=ncread(file,(vars{ii}),[1,readStart(1)],[length(data.range),sampleNum(1)]);
    else
        data.(vars{ii})=ncread(file,(vars{ii}),readStart(1),sampleNum(1));
    end
    if ~strcmp(vars{ii},'range') & size(data.(vars{ii}),2)==1
        data.(vars{ii})=data.(vars{ii})';
    end
end

for jj=2:length(fileList)
    file=fileList{jj};

    baseTime=ncread(file,'base_time');
    timeOffset=ncread(file,'time_offset_vc',readStart(jj),sampleNum(jj));

    fileStartTime=datetime(1970,1,1)+seconds(baseTime);
    data.time=cat(2,data.time,fileStartTime+seconds(timeOffset)');

    data.elevation=cat(2,data.elevation,ncread(file,'elevation_vc',readStart(jj),sampleNum(jj))');
    data.latitude=cat(2,data.latitude,ncread(file,'latitude',readStart(jj),sampleNum(jj))');
    data.longitude=cat(2,data.longitude,ncread(file,'longitude',readStart(jj),sampleNum(jj))');
    data.altitude=cat(2,data.altitude,ncread(file,'altitude',readStart(jj),sampleNum(jj))');
    data.pulse_width=cat(2,data.pulse_width,ncread(file,'pulse_width_vc',readStart(jj),sampleNum(jj))');
    data.prt=cat(2,data.prt,ncread(file,'prt_vc',readStart(jj),sampleNum(jj))');

    for ii=1:length(vars)
        if dim2D(ii)
            readThis=ncread(file,(vars{ii}),[1,readStart(jj)],[length(data.range),sampleNum(jj)]);
        else
            readThis=ncread(file,(vars{ii}),readStart(jj),sampleNum(jj));
        end
        if ~strcmp(vars{ii},'range') & size(readThis,2)==1
            readThis=readThis';
        end
        data.(vars{ii})=cat(2,data.(vars{ii}),single(readThis));
    end
end

end