% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.1';
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
    upNumAll.(classTypes{ii})=[];
    upFracAll.(classTypes{ii})=[];
    upMaxWidthAll.(classTypes{ii})=[];
    upMaxDepthAll.(classTypes{ii})=[];
    upMaxStrengthAll.(classTypes{ii})=[];
    downMaxStrengthAll.(classTypes{ii})=[];
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

        velMap=velBig(min(clR):max(clR),min(clC):max(clC));
        [upNum,upFrac,upMaxWidth,upMaxDepth,upMaxStrength,downMaxStrength]=upDownDraft(velMap,data.range(2)-data.range(1),mean(groundDist(clC)));

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
        upNumAll.(classTypes{cloudType})=cat(1,upNumAll.(classTypes{cloudType}),upNum);
        upFracAll.(classTypes{cloudType})=cat(1,upFracAll.(classTypes{cloudType}),upFrac);
        upMaxWidthAll.(classTypes{cloudType})=cat(1,upMaxWidthAll.(classTypes{cloudType}),upMaxWidth);
        upMaxDepthAll.(classTypes{cloudType})=cat(1,upMaxDepthAll.(classTypes{cloudType}),upMaxDepth);
        upMaxStrengthAll.(classTypes{cloudType})=cat(1,upMaxStrengthAll.(classTypes{cloudType}),upMaxStrength);
        downMaxStrengthAll.(classTypes{cloudType})=cat(1,downMaxStrengthAll.(classTypes{cloudType}),downMaxStrength);
    end
end

%% Save properties

disp('Saving output ...');

if ~strcmp(project,'spicule')
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','sstAll','upNumAll','upFracAll','upMaxWidthAll','upMaxDepthAll','upMaxStrengthAll','downMaxStrengthAll');
else
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','upNumAll','upFracAll','upMaxWidthAll','upMaxDepthAll','upMaxStrengthAll','downMaxStrengthAll');
end
%% Plot

disp('Plotting ...')

colmapCC=[204,255,204;
    153,204,0;
    0,128,0;
    0,204,255;
    51,102,255;
    0,0,180;
    255,204,0;
    255,102,0;
    220,0,0;
    255,153,220;
    204,153,255;
    128,0,128];

colmapCC=colmapCC./255;

close all

%% Max refl

close all

edges=-50:30;
xlab='Maximum reflectivity (dBZ)';
figname=[figdir,project,'_maxRefl.png'];

plotStats(maxReflAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean refl

close all

edges=-50:30;
xlab='Mean reflectivity (dBZ)';
figname=[figdir,project,'_meanRefl.png'];

plotStats(meanReflAll,edges,xlab,figname,classTypes,colmapCC);

%% Max convectivity

close all

edges=0:0.02:1;
xlab='Max convectivity';
figname=[figdir,project,'_maxConv.png'];

plotStats(maxConvAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean convectivity

close all

edges=0:0.02:1;
xlab='Mean convectivity';
figname=[figdir,project,'_meanConv.png'];

plotStats(meanConvAll,edges,xlab,figname,classTypes,colmapCC);

%% Cloud depth

close all

edges=0:0.2:15;
xlab='Cloud depth (km)';
figname=[figdir,project,'_cloudDepth.png'];

plotStats(cloudDepthAll,edges,xlab,figname,classTypes,colmapCC);

%% Max temperature

close all

edges=-80:2:30;
xlab='Max temperature (C)';
figname=[figdir,project,'_maxTemp.png'];

plotStats(maxTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Min temperature

close all

edges=-80:2:30;
xlab='Min temperature (C)';
figname=[figdir,project,'_minTemp.png'];

plotStats(minTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean temperature

close all

edges=-80:2:30;
xlab='Mean temperature (C)';
figname=[figdir,project,'_meanTemp.png'];

plotStats(meanTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Max pressure

close all

edges=100:10:1100;
xlab='Max pressure (hPa)';
figname=[figdir,project,'_maxPress.png'];

plotStats(maxPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Min pressure

close all

edges=100:10:1100;
xlab='Min pressure (hPa)';
figname=[figdir,project,'_minPress.png'];

plotStats(minPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean pressure

close all

edges=100:10:1100;
xlab='Mean pressure (hPa)';
figname=[figdir,project,'_meanPress.png'];

plotStats(meanPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Icing level

close all

edges=0:0.1:8;
xlab='Icing level (km)';
figname=[figdir,project,'_iceLev.png'];

plotStats(iceLevAll,edges,xlab,figname,classTypes,colmapCC);

%% SST

if ~strcmp(project,'spicule')
    close all

    edges=-35:1:35;
    xlab='SST (C)';
    figname=[figdir,project,'_sst.png'];

    plotStats(sstAll,edges,xlab,figname,classTypes,colmapCC);
end

%% Longitude

close all

if strcmp(project,'cset')
    edges=-160:0.5:-120;
elseif strcmp(project,'socrates')
    edges=130:0.5:180;
elseif strcmp(project,'otrec')
    edges=-95:0.5:-65;
elseif strcmp(project,'spicule')
    edges=-140:0.5:-60;
end
xlab='Longitude (deg)';
figname=[figdir,project,'_lon.png'];

plotStats(lonAll,edges,xlab,figname,classTypes,colmapCC);

%% Latitude

close all

if strcmp(project,'cset')
    edges=15:0.5:50;
elseif strcmp(project,'socrates')
    edges=-70:0.5:-30;
elseif strcmp(project,'otrec')
    edges=-5:0.5:15;
elseif strcmp(project,'spicule')
    edges=20:0.5:55;
end
xlab='Latitude (deg)';
figname=[figdir,project,'_lat.png'];

plotStats(latAll,edges,xlab,figname,classTypes,colmapCC);

%% Updraft number

close all

edges=1:1:30;
xlab='Number of updrafts';
figname=[figdir,project,'_upNum.png'];

plotStats(upNumAll,edges,xlab,figname,classTypes,colmapCC);

%% Updraft fraction

close all

edges=0:0.02:1;
xlab='Updraft fraction';
figname=[figdir,project,'_upFrac.png'];

plotStats(upFracAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum width of updrafts

close all

edges=0:2:100;
xlab='Max width of updrafts (km)';
figname=[figdir,project,'_upMaxWidth.png'];

plotStats(upMaxWidthAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum depth of updrafts

close all

edges=0:0.1:10;
xlab='Max depth of updrafts (km)';
figname=[figdir,project,'_upMaxDepth.png'];

plotStats(upMaxDepthAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum strength of updraft

close all

edges=0:0.5:20;
xlab='Max up velocity (m s^{-1})';
figname=[figdir,project,'_upMaxStrength.png'];

plotStats(upMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum strength of downdraft

close all

edges=0:0.5:20;
xlab='Max down velocity (m s^{-1})';
figname=[figdir,project,'_downMaxStrength.png'];

plotStats(downMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);
