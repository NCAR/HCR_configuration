function [model] = read_model_longShort(model,modelDir,startTime,endTime)
% Read in model data
% Input variables:
% model: structure with empty variables to read
% Possible variables are time, p, t, rh, u, v, sst, asl
% modelDir: directory where the model data is located
% startTime, endTime: requested time frame (in datetime format) 
model.time=[];
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
    allFiles=dir([modelDir,'*',varNames{ii},'*.mat']);
    underSc=strfind(allFiles(1).name,'_');
    
    fileStart=[];
    fileEnd=[];
    for jj=1:size(allFiles,1)
        fileStart=cat(1,fileStart,datetime(str2num(allFiles(jj).name(underSc(2)-8:underSc(2)-5)),...
            str2num(allFiles(jj).name(underSc(2)-4:underSc(2)-3)),...
            str2num(allFiles(jj).name(underSc(2)-2:underSc(2)-1)),...
            str2num(allFiles(jj).name(underSc(2)+1:underSc(2)+2)),...
            str2num(allFiles(jj).name(underSc(2)+3:underSc(2)+4)),...
            str2num(allFiles(jj).name(underSc(2)+5:underSc(2)+6)))-minutes(10));
        fileEnd=cat(1,fileEnd,datetime(str2num(allFiles(1).name(underSc(5)-8:underSc(5)-5)),...
            str2num(allFiles(jj).name(underSc(5)-4:underSc(5)-3)),...
            str2num(allFiles(jj).name(underSc(5)-2:underSc(5)-1)),...
            str2num(allFiles(jj).name(underSc(5)+1:underSc(5)+2)),...
            str2num(allFiles(jj).name(underSc(5)+3:underSc(5)+4)),...
            str2num(allFiles(jj).name(underSc(5)+5:underSc(5)+6)))+minutes(10));
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
timeInds=find(modelTemp.time.timeHCR>=startTime & modelTemp.time.timeHCR<=endTime);

for ii=1:length(varNames)
    nameIn=fields(modelTemp.(varNames{ii}));
    dataIn=modelTemp.(varNames{ii}).(nameIn{:});
    if size(dataIn,2)==1
        dataIn=dataIn';
    end
    model.(varNamesIn{ii})=dataIn(:,timeInds);
end
end

