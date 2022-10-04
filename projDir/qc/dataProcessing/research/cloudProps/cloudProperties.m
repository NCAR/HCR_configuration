% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';
whichModel='era5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

cfDir=HCRdir(project,quality,qcVersion,freqData);

[~,matDir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[cfDir(1:end-5),'cloudPuzzleEchoType/wholeFlights/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

matFilesClass=dir([matDir,'*cloudClass*']);
matFilesPuzzle=dir([matDir,'*cloudPuzzle*']);

classTypes={'CloudLow','CloudMid','CloudHigh',...
        'StratShallow','StratMid','StratDeep',...
        'ConvYoungShallow','ConvYoungMid','ConvYongDeep',...
        'ConvMatureShallow','ConvMatureMid','ConvMatureDeep'};

%% Initiate output
for ii=1:length(classTypes)
    maxReflAll.(classTypes{ii})=[];
    meanReflAll.(classTypes{ii})=[];
    maxConvAll.(classTypes{ii})=[];
    meanConvAll.(classTypes{ii})=[];
    cloudDepthAll.(classTypes{ii})=[];
end

%% Loop through flights

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading cloud class data ...')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    cloudClass=load([matDir,matFilesClass(aa).name]);
    cloudClass=cloudClass.cloudClass;

    cloudPuzzle=load([matDir,matFilesPuzzle(aa).name]);
    cloudPuzzle=cloudPuzzle.cloudPuzzle;

    %% Get cfradial data

    disp("Loading cfradial data ...");

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    fileList=makeFileList(cfDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ_MASKED=[];
    data.CONVECTIVITY=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    % Check time
    if ~isequal(size(cloudClass),size(data.DBZ_MASKED))
        error('Times do not match up.')
    end

    %% Loop through clouds

    cloudClass(cloudClass==11)=4;
    cloudClass(cloudClass==12)=5;
    cloudClass(cloudClass==13)=6;
    cloudClass(cloudClass==21)=7;
    cloudClass(cloudClass==22)=8;
    cloudClass(cloudClass==23)=9;
    cloudClass(cloudClass==31)=10;
    cloudClass(cloudClass==32)=11;
    cloudClass(cloudClass==33)=12;

    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));

    for ii=1:length(uClouds)
        disp(['Processing cloud ',num2str(ii),' of ',num2str(length(uClouds)),' ...'])

        cloudInds=find(cloudPuzzle==uClouds(ii));

        % Max and average refl
        maxRefl=max(data.DBZ_MASKED(cloudInds),[],'omitnan');
        meanRefl=mean(data.DBZ_MASKED(cloudInds),'omitnan');

        % Max and average convectivity
        maxConv=max(data.CONVECTIVITY(cloudInds),[],'omitnan');
        meanConv=mean(data.CONVECTIVITY(cloudInds),'omitnan');

        % Cloud depth
        cloudAsl=data.asl(cloudInds);
        cloudDepth=max(cloudAsl,[],'omitnan')-min(cloudAsl,[],'omitnan');

        % Process mat
        [clR,clC]=ind2sub(size(data.DBZ_MASKED),cloudInds);

        maskBig=zeros(size(data.DBZ_MASKED));
        maskBig(cloudInds)=1;

        typeMap=maskBig(min(clR):max(clR),min(clC):max(clC));

        % Add output
        cloudType=unique(cloudClass(cloudInds));

        if cloudType==0
            continue
        end

        maxReflAll.(classTypes{cloudType})=cat(1,maxReflAll.(classTypes{cloudType}),maxRefl);
        meanReflAll.(classTypes{cloudType})=cat(1,meanReflAll.(classTypes{cloudType}),meanRefl);
        maxConvAll.(classTypes{cloudType})=cat(1,maxConvAll.(classTypes{cloudType}),maxConv);
        meanConvAll.(classTypes{cloudType})=cat(1,meanConvAll.(classTypes{cloudType}),meanConv);
        cloudDepthAll.(classTypes{cloudType})=cat(1,cloudDepthAll.(classTypes{cloudType}),cloudDepth);
    end
end

%% Plot
% fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900],'visible',showPlot);
%
%
% formatOut = 'yyyymmdd_HHMM';
% set(gcf,'PaperPositionMode','auto')
% print([figdir,project,'_.png'],'-dpng','-r0');
