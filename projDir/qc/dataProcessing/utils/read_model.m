function [model] = read_model(model,modelDir,startTime,endTime)
% Read in model data
% Input variables:
% model: structure with empty variables to read
% Possible variables are time, p, t, rh, u, v, sst, asl
% modelDir: directory where the model data is located
% startTime, endTime: requested time frame (in datetime format) 

varNamesIn=fields(model);
if ~strcmp(varNamesIn{1},'momentsSpecParams')
    model.time=[];
end

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
        fileStart=cat(1,fileStart,datetime(str2num(allFiles(jj).name(underSc(1)-8:underSc(1)-5)),...
            str2num(allFiles(jj).name(underSc(1)-4:underSc(1)-3)),...
            str2num(allFiles(jj).name(underSc(1)-2:underSc(1)-1)),...
            str2num(allFiles(jj).name(underSc(1)+1:underSc(1)+2)),...
            str2num(allFiles(jj).name(underSc(1)+3:underSc(1)+4)),...
            str2num(allFiles(jj).name(underSc(1)+5:underSc(1)+6)))-minutes(10));
        fileEnd=cat(1,fileEnd,datetime(str2num(allFiles(1).name(underSc(4)-8:underSc(4)-5)),...
            str2num(allFiles(jj).name(underSc(4)-4:underSc(4)-3)),...
            str2num(allFiles(jj).name(underSc(4)-2:underSc(4)-1)),...
            str2num(allFiles(jj).name(underSc(4)+1:underSc(4)+2)),...
            str2num(allFiles(jj).name(underSc(4)+3:underSc(4)+4)),...
            str2num(allFiles(jj).name(underSc(4)+5:underSc(4)+6)))+minutes(10));
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
if ~strcmp(varNamesIn{1},'momentsSpecParams')
    timeInds=find(modelTemp.time.timeHCR>=startTime & modelTemp.time.timeHCR<=endTime);

    for ii=1:length(varNames)
        nameIn=fields(modelTemp.(varNames{ii}));
        dataIn=modelTemp.(varNames{ii}).(nameIn{:});
        if size(dataIn,2)==1
            dataIn=dataIn';
        end
        model.(varNamesIn{ii})=dataIn(:,timeInds);
    end
else
    timeInds=find(modelTemp.momentsSpecParams.momentsSpecParams.time>=startTime & ...
        modelTemp.momentsSpecParams.momentsSpecParams.time<=endTime);

    varNames=fields(modelTemp.momentsSpecParams.momentsSpecParams);

    for ii=1:length(varNames)
        dataIn=modelTemp.momentsSpecParams.momentsSpecParams.(varNames{ii});
        % if size(dataIn,2)==1
        %     dataIn=dataIn';
        % end
        model.(varNames{ii})=dataIn(:,timeInds);
    end
end
end

