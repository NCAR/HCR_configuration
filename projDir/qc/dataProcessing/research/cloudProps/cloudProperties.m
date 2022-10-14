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

figdir=[cfDir(1:end-5),'cloudProps/'];

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
    maxTempAll.(classTypes{ii})=[];
    minTempAll.(classTypes{ii})=[];
    meanTempAll.(classTypes{ii})=[];
    maxPressAll.(classTypes{ii})=[];
    minPressAll.(classTypes{ii})=[];
    meanPressAll.(classTypes{ii})=[];
    iceLevAll.(classTypes{ii})=[];
    if ~strcmp(project,'spicule')
        sstAll.(classTypes{ii})=[];
    end
    lonAll.(classTypes{ii})=[];
    latAll.(classTypes{ii})=[];
    upFracAll.(classTypes{ii})=[];
    upMeanStrengthAll.(classTypes{ii})=[];
    downMeanStrengthAll.(classTypes{ii})=[];
    upMaxStrengthAll.(classTypes{ii})=[];
    downMaxStrengthAll.(classTypes{ii})=[];
    upRegsAll.(classTypes{ii})=[];
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
    data.VEL_MASKED=[];
    data.CONVECTIVITY=[];
    data.TEMP=[];
    data.PRESS=[];
    data.ICING_LEVEL=[];
    if ~strcmp(project,'spicule')
        data.SST=[];
    end
    data.eastward_velocity=[];
    data.northward_velocity=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.CONVECTIVITY(data.CONVECTIVITY>1)=1;
    data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

    groundSpeed=sqrt(data.eastward_velocity.^2+data.northward_velocity.^2);
    groundDist=groundSpeed.*etime(datevec(data.time(2)),datevec(data.time(1)));

    % Check time
    if ~isequal(size(cloudClass),size(data.DBZ_MASKED))
        error('Times do not match up.')
    end

    %% Prepare cloud Class

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

    %% Loop through clouds

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
        cloudAsl=data.asl(cloudInds)./1000;
        cloudDepth=max(cloudAsl,[],'omitnan')-min(cloudAsl,[],'omitnan');

        % Temperature
        maxTemp=max(data.TEMP(cloudInds),[],'omitnan');
        minTemp=min(data.TEMP(cloudInds),[],'omitnan');
        meanTemp=mean(data.TEMP(cloudInds),'omitnan');

        % Pressure
        maxPress=max(data.PRESS(cloudInds),[],'omitnan');
        minPress=min(data.PRESS(cloudInds),[],'omitnan');
        meanPress=mean(data.PRESS(cloudInds),'omitnan');

        % 1D
        [clR,clC]=ind2sub(size(data.DBZ_MASKED),cloudInds);

        % Icing level
        iceLev=mean(data.ICING_LEVEL(clC)./1000,'omitnan');

        % Sea surface temperature
        if ~strcmp(project,'spicule')
            sst=mean(data.SST(clC),'omitnan');
        end

        % Longitude and latitude
        lon=mean(data.longitude(clC),'omitnan');
        lat=mean(data.latitude(clC),'omitnan');

        % Velocity
        velBig=nan(size(data.VEL_MASKED));
        velBig(cloudInds)=data.VEL_MASKED(cloudInds);

        aslBig=nan(size(data.asl));
        aslBig(cloudInds)=data.asl(cloudInds);

        velMap=velBig(min(clR):max(clR),min(clC):max(clC));
        aslMap=aslBig(min(clR):max(clR),min(clC):max(clC));
        
        [upRegs,upFrac,upMaxStrength,downMaxStrength,upMeanStrength,downMeanStrength]=upDownDraft(velMap,aslMap,data.range(2)-data.range(1),mean(groundDist(clC)));
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
        maxTempAll.(classTypes{cloudType})=cat(1,maxTempAll.(classTypes{cloudType}),maxTemp);
        minTempAll.(classTypes{cloudType})=cat(1,minTempAll.(classTypes{cloudType}),minTemp);
        meanTempAll.(classTypes{cloudType})=cat(1,meanTempAll.(classTypes{cloudType}),meanTemp);
        maxPressAll.(classTypes{cloudType})=cat(1,maxPressAll.(classTypes{cloudType}),maxPress);
        minPressAll.(classTypes{cloudType})=cat(1,minPressAll.(classTypes{cloudType}),minPress);
        meanPressAll.(classTypes{cloudType})=cat(1,meanPressAll.(classTypes{cloudType}),meanPress);
        iceLevAll.(classTypes{cloudType})=cat(1,iceLevAll.(classTypes{cloudType}),iceLev);
        if ~strcmp(project,'spicule')
            sstAll.(classTypes{cloudType})=cat(1,sstAll.(classTypes{cloudType}),sst);
        end
        lonAll.(classTypes{cloudType})=cat(1,lonAll.(classTypes{cloudType}),lon);
        latAll.(classTypes{cloudType})=cat(1,latAll.(classTypes{cloudType}),lat);
        upFracAll.(classTypes{cloudType})=cat(1,upFracAll.(classTypes{cloudType}),upFrac);
        upMeanStrengthAll.(classTypes{cloudType})=cat(1,upMeanStrengthAll.(classTypes{cloudType}),upMeanStrength);
        downMeanStrengthAll.(classTypes{cloudType})=cat(1,downMeanStrengthAll.(classTypes{cloudType}),downMeanStrength);
        upMaxStrengthAll.(classTypes{cloudType})=cat(1,upMaxStrengthAll.(classTypes{cloudType}),upMaxStrength);
        downMaxStrengthAll.(classTypes{cloudType})=cat(1,downMaxStrengthAll.(classTypes{cloudType}),downMaxStrength);
        upRegsAll.(classTypes{cloudType})=cat(1,upRegsAll.(classTypes{cloudType}),upRegs);
    end
end

%% Save properties

disp('Saving output ...');

if ~strcmp(project,'spicule')
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','sstAll','upFracAll','upRegsAll', ...
        'upMeanStrengthAll','downMeanStrengthAll','upMaxStrengthAll','downMaxStrengthAll','latAll','lonAll');
else
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','upRegNumAll','upFracAll','upRegsAll', ...
        'upMeanStrengthAll','downMeanStrengthAll','upMaxStrengthAll','downMaxStrengthAll','latAll','lonAll');
end
