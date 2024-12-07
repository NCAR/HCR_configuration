function [model] = read_model_longShort(model,modelDir,startTime,endTime)
% Read in model data
% Input variables:
% model: structure with empty variables to read
% Possible variables are time, p, t, rh, u, v, sst, asl
% modelDir: directory where the model data is located
% startTime, endTime: requested time frame (in datetime format) 
model.time_long=[];
varNamesIn=fields(model);

varNames={};
for ii=1:length(varNamesIn)
    if strcmp(varNamesIn{ii},'p') | strcmp(varNamesIn{ii},'rh') ...
            | strcmp(varNamesIn{ii},'temp') | strcmp(varNamesIn{ii},'u') | strcmp(varNamesIn{ii},'v')
        varNames{end+1}=[varNamesIn{ii} 'HCR'];
    else
        varNames{end+1}=varNamesIn{ii};
    end
end

fileList={};

for ii=1:length(varNames)
    allFiles=dir([modelDir,'*',varNames{ii},'.*.mat']);
    dot=strfind(allFiles(1).name,'.');
    
    fileStart=[];
    fileEnd=[];
    for jj=1:size(allFiles,1)
        fileStart=cat(1,fileStart,datetime(str2num(allFiles(jj).name(dot(2)+1:dot(2)+4)),...
            str2num(allFiles(jj).name(dot(2)+5:dot(2)+6)),...
            str2num(allFiles(jj).name(dot(2)+7:dot(2)+8)),...
            str2num(allFiles(jj).name(dot(2)+10:dot(2)+11)),...
            str2num(allFiles(jj).name(dot(2)+12:dot(2)+13)),...
            str2num(allFiles(jj).name(dot(2)+14:dot(2)+15)))-minutes(10));
        fileEnd=cat(1,fileEnd,datetime(str2num(allFiles(1).name(dot(3)-15:dot(3)-12)),...
            str2num(allFiles(jj).name(dot(3)-11:dot(3)-10)),...
            str2num(allFiles(jj).name(dot(3)-9:dot(3)-8)),...
            str2num(allFiles(jj).name(dot(3)-6:dot(3)-5)),...
            str2num(allFiles(jj).name(dot(3)-4:dot(3)-3)),...
            str2num(allFiles(jj).name(dot(3)-2:dot(3)-1)))+minutes(10));
    end
    fileInd=max(find(fileStart<startTime+hours(1)));
    if endTime<fileEnd(fileInd)+seconds(1) | abs(fileEnd(fileInd)-endTime)<hours(1)
        fileOut=[allFiles(fileInd).folder,'/',allFiles(fileInd).name];
        fileList{end+1}=fileOut;
    else
        disp(['No ',varNames{ii},' file found.']);
    end
end

if length(fileList)~=length(varNamesIn)
    disp('Some model variables were not found.');
    return
end

% Load data
for ii=1:length(varNames)
    modelTemp.(varNames{ii})=load(fileList{ii});
end

% Get right times
timeInds=find(modelTemp.time_long.timeHCR>=startTime & modelTemp.time_long.timeHCR<=endTime);

for ii=1:length(varNames)
    nameIn=fields(modelTemp.(varNames{ii}));
    dataIn=modelTemp.(varNames{ii}).(nameIn{:});
    if size(dataIn,2)==1
        dataIn=dataIn';
    end
    model.(varNamesIn{ii})=dataIn(:,timeInds);
end
model.time=model.time_long;
model=rmfield(model,'time_long');
end

