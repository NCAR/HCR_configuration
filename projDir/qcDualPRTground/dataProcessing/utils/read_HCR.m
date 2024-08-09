function [data] = read_HCR(fileList,indata,startTime,endTime)

vars = fieldnames(indata);

% Read first file and figure out format
infile=fileList{1};

startTimeIn=ncread(infile,'time_coverage_start')';
startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
    str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
timeRead=ncread(infile,'time')';
timeIn=startTimeFile+seconds(timeRead);
indata.time=timeIn;
indata.altitude=ncread(infile,'altitude');
indata.latitude=ncread(infile,'latitude');
indata.longitude=ncread(infile,'longitude');
indata.elevation=ncread(infile,'elevation');

for ii=1:size(vars,1)
    try
        indata.(vars{ii})=ncread(infile,vars{ii});
    catch
        disp(['Variable ' vars{ii} ' does not exist in CfRadial file ',infile]);
        indata=rmfield(indata,vars{ii});
    end
end

allVars=fieldnames(indata);
timeDim=size(indata.time);
if max(timeDim)==1
    disp(['Time dimension is one in file ',infile]);
end

if timeDim(1)>timeDim(2)
    flip.time=1;
    indata.time=indata.time';
else
    flip.time=0;
end

for ii=1:size(allVars,1)
    if ~strcmp((allVars{ii}),'time')
        dim=size(indata.(allVars{ii}));
        if dim(2)~=size(indata.time,2)
            indata.(allVars{ii})=indata.(allVars{ii})';
            flip.(allVars{ii})=1;
        else
            flip.(allVars{ii})=0;
        end
    end
end

indata.range=ncread(infile,'range');

if size(indata.range,1)==1
    indata.range=indata.range';
    flip.range=1;
else
    flip.range=0;
end

indata.range=repmat(indata.range,1,size(indata.time,2));

allVars=fieldnames(indata);

for jj=2:length(fileList)
    infile=fileList{jj};
    
    startTimeIn=ncread(infile,'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeIn=ncread(infile,'time')';
    timeIn=startTimeFile+seconds(timeIn);
    tempTime=timeIn;
    
    if flip.time
        tempTime=tempTime';
    end
    
    indata.time=cat(2,indata.time,tempTime);
    timeLength=size(tempTime,2);
    
    for ii=1:size(allVars,1)
        if ~strcmp((allVars{ii}),'time')
            temp=ncread(infile,allVars{ii});
            if flip.(allVars{ii})
                temp=temp';
            end
            if strcmp((allVars{ii}),'range')
                temp=repmat(temp,1,timeLength);
            end
            try
                indata.(allVars{ii})=cat(2,indata.(allVars{ii}),temp);
            catch
                stop1=1;
            end
        end
    end
    if size(indata.time,2)~=size(indata.(allVars{4}),2)
        disp(['Time and data length do not match up at ',datestr(indata.time(end),'yyyy-mm-dd HH:MM:SS')]);
        disp(['File ',infile]);
        error('Stopping.')
    end
end

timeInds=find(indata.time>=startTime & indata.time<=endTime);
for ii=1:size(allVars,1)
    if isfield(indata,allVars{ii})
        if min(size(indata.(allVars{ii})))~=1
            data.(allVars{ii})=single(indata.(allVars{ii})(:,timeInds));
        else
            data.(allVars{ii})=indata.(allVars{ii})(:,timeInds);
        end
    end
end
try
    data.asl=HCRrange2asl(data.range,data.elevation,data.altitude);
    data.asl=single(data.asl);
end
end

