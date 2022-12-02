% Calculate cloud properties

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';
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
        'ConvYoungShallow','ConvYoungMid','ConvYoungDeep',...
        'ConvMatureShallow','ConvMatureMid','ConvMatureDeep'};

%% Initiate output
for ii=1:length(classTypes)
    maxReflAll.(classTypes{ii})=[];
    meanReflAll.(classTypes{ii})=[];
    maxConvAll.(classTypes{ii})=[];
    meanConvAll.(classTypes{ii})=[];
    cloudDepthAll.(classTypes{ii})=[];
    cloudLengthAll.(classTypes{ii})=[];
    cloudTopAll.(classTypes{ii})=[];
    cloudBaseAll.(classTypes{ii})=[];
    cloudLayersAll.(classTypes{ii})=[];
    maxTempAll.(classTypes{ii})=[];
    minTempAll.(classTypes{ii})=[];
    meanTempAll.(classTypes{ii})=[];
    maxPressAll.(classTypes{ii})=[];
    minPressAll.(classTypes{ii})=[];
    meanPressAll.(classTypes{ii})=[];
    iceLevAll.(classTypes{ii})=[];
    divLevAll.(classTypes{ii})=[];
    meltDetAll.(classTypes{ii})=[];
    coldFracAll.(classTypes{ii})=[];
    if ~strcmp(project,'spicule')
        sstAll.(classTypes{ii})=[];
    end
    lonAll.(classTypes{ii})=[];
    latAll.(classTypes{ii})=[];
    upFracAll.(classTypes{ii})=[];
    upNumAll.(classTypes{ii})=[];
    upMeanStrengthAll.(classTypes{ii})=[];
    downMeanStrengthAll.(classTypes{ii})=[];
    upMaxStrengthAll.(classTypes{ii})=[];
    downMaxStrengthAll.(classTypes{ii})=[];
    upRegsAll.(classTypes{ii})=[];
    precShaftsAll.(classTypes{ii})=[];
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
    data.MELTING_LAYER=[];
    data.FLAG=[];
    if ~strcmp(project,'spicule')
        data.SST=[];
    end
    data.eastward_velocity=[];
    data.northward_velocity=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.CONVECTIVITY(data.CONVECTIVITY>1)=1;
    data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

    % Shrink velocity to remove outliers at the edges
    VELmask=zeros(size(data.VEL_MASKED));
    VELmask(~isnan(data.VEL_MASKED))=1;

    VELsmall=imerode(VELmask,strel('disk',5));

    data.VEL_MASKED(VELmask==1 & VELsmall==0)=nan;

    % Get distance traveled
    groundSpeed=sqrt(data.eastward_velocity.^2+data.northward_velocity.^2);
    groundDist=groundSpeed.*etime(datevec(data.time(2)),datevec(data.time(1)));

    % Prepare precip shafts
    shaftMask=~isnan(data.DBZ_MASKED);
    shaftMask(data.FLAG==3)=1;

    % Count cloud layers
    cloudMask=~isnan(cloudPuzzle);
    numLayers=zeros(1,size(cloudMask,2));
    for kk=1:size(cloudMask,2)
        maskRay=cloudMask(:,kk);
        maskRay=bwareaopen(maskRay,5);
        cloudsRay=bwconncomp(maskRay);
        numLayers(kk)=cloudsRay.NumObjects;
    end
    numLayers(data.altitude<5000 & data.elevation<0)=nan;
    numLayers(data.altitude>500 & data.elevation>0)=nan;

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
        
        % Check cloud type and size
        cloudInds=find(cloudPuzzle==uClouds(ii));
        cloudType=unique(cloudClass(cloudInds));

        if cloudType==0 | length(cloudInds)<5000
            continue
        end

        disp(['Processing cloud ',num2str(ii),' of ',num2str(length(uClouds)),' ...'])

        % Max and average refl
        maxRefl=max(data.DBZ_MASKED(cloudInds),[],'omitnan');
        meanRefl=mean(data.DBZ_MASKED(cloudInds),'omitnan');

        % Max and average convectivity
        maxConv=max(data.CONVECTIVITY(cloudInds),[],'omitnan');
        meanConv=mean(data.CONVECTIVITY(cloudInds),'omitnan');

        % Cloud depth
        cloudAsl=data.asl(cloudInds)./1000;
        cloudDepth=max(cloudAsl,[],'omitnan')-min(cloudAsl,[],'omitnan');

        % Cloud top and base
        cloudTop=max(cloudAsl,[],'omitnan');
        cloudBase=nan;
        if cloudType<=3
            cloudBase=min(cloudAsl,[],'omitnan');
        end

        % Melting layer detection
        meltLayerCloud=data.MELTING_LAYER(cloudInds);
        meltDet=double(max(ismember([12,22],meltLayerCloud)));

        % Cold fraction
        coldFrac=length(find(meltLayerCloud>=20))./length(meltLayerCloud);

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

        % Cloud length
        cloudLength=(mean(groundDist(clC))*(max(clC)-min(clC)))/1000;

        % Icing level
        iceLev=mean(data.ICING_LEVEL(clC)./1000,'omitnan');

        % Divergence level
        tempCols=data.TEMP(:,clC);
        aslCols=data.asl(:,clC);

        aslCols(tempCols>-24.8 | tempCols<-25.2 | isnan(tempCols))=nan;
        divLev=mean(aslCols(:),'omitnan')./1000;

        % Sea surface temperature
        if ~strcmp(project,'spicule')
            sst=mean(data.SST(clC),'omitnan');
        end

        % Longitude and latitude
        lon=mean(data.longitude(clC),'omitnan');
        lat=mean(data.latitude(clC),'omitnan');

        % Number of cloud layers
        cloudLayers=median(numLayers(clC));

        % Cloud mask area
        bigMask=zeros(size(cloudPuzzle));
        bigMask(cloudInds)=1;
        smallMask=bigMask(min(clR):max(clR),min(clC):max(clC));
        smallInds=find(smallMask==1);

        % Velocity
        velMap=nan(size(smallMask));
        velMap(smallInds)=data.VEL_MASKED(cloudInds);

        aslMap=nan(size(smallMask));
        aslMap(smallInds)=data.asl(cloudInds);
        
        [upRegs,upFrac,upNum,upMaxStrength,downMaxStrength,upMeanStrength,downMeanStrength]=upDownDraft(velMap,aslMap,data.range(2)-data.range(1),mean(groundDist(clC)),data.longitude(min(clC):max(clC)),data.latitude(min(clC):max(clC)));
        
        % Precip shafts
        shaftMap=zeros(size(smallMask));
        shaftMap(smallInds)=shaftMask(cloudInds);

        dbzMap=nan(size(smallMask));
        dbzMap(smallInds)=data.DBZ_MASKED(cloudInds);

        precShafts=table(nan,nan,nan,nan,nan,nan, ...
            'VariableNames',{'shaftKM','frac','meanRef','maxRefl','meanVel','maxVel'});
        if cloudType>3
            precShafts=precipShafts(shaftMap,dbzMap,velMap,aslMap,mean(groundDist(clC)));
        end
        
        % Add output
        
        maxReflAll.(classTypes{cloudType})=cat(1,maxReflAll.(classTypes{cloudType}),maxRefl);
        meanReflAll.(classTypes{cloudType})=cat(1,meanReflAll.(classTypes{cloudType}),meanRefl);
        maxConvAll.(classTypes{cloudType})=cat(1,maxConvAll.(classTypes{cloudType}),maxConv);
        meanConvAll.(classTypes{cloudType})=cat(1,meanConvAll.(classTypes{cloudType}),meanConv);
        cloudDepthAll.(classTypes{cloudType})=cat(1,cloudDepthAll.(classTypes{cloudType}),cloudDepth);
        cloudLengthAll.(classTypes{cloudType})=cat(1,cloudLengthAll.(classTypes{cloudType}),cloudLength);
        cloudTopAll.(classTypes{cloudType})=cat(1,cloudTopAll.(classTypes{cloudType}),cloudTop);
        cloudBaseAll.(classTypes{cloudType})=cat(1,cloudBaseAll.(classTypes{cloudType}),cloudBase);
        cloudLayersAll.(classTypes{cloudType})=cat(1,cloudLayersAll.(classTypes{cloudType}),cloudLayers);
        maxTempAll.(classTypes{cloudType})=cat(1,maxTempAll.(classTypes{cloudType}),maxTemp);
        minTempAll.(classTypes{cloudType})=cat(1,minTempAll.(classTypes{cloudType}),minTemp);
        meanTempAll.(classTypes{cloudType})=cat(1,meanTempAll.(classTypes{cloudType}),meanTemp);
        maxPressAll.(classTypes{cloudType})=cat(1,maxPressAll.(classTypes{cloudType}),maxPress);
        minPressAll.(classTypes{cloudType})=cat(1,minPressAll.(classTypes{cloudType}),minPress);
        meanPressAll.(classTypes{cloudType})=cat(1,meanPressAll.(classTypes{cloudType}),meanPress);
        iceLevAll.(classTypes{cloudType})=cat(1,iceLevAll.(classTypes{cloudType}),iceLev);
        divLevAll.(classTypes{cloudType})=cat(1,divLevAll.(classTypes{cloudType}),divLev);
        meltDetAll.(classTypes{cloudType})=cat(1,meltDetAll.(classTypes{cloudType}),meltDet);
        coldFracAll.(classTypes{cloudType})=cat(1,coldFracAll.(classTypes{cloudType}),coldFrac);
        if ~strcmp(project,'spicule')
            sstAll.(classTypes{cloudType})=cat(1,sstAll.(classTypes{cloudType}),sst);
        end
        lonAll.(classTypes{cloudType})=cat(1,lonAll.(classTypes{cloudType}),lon);
        latAll.(classTypes{cloudType})=cat(1,latAll.(classTypes{cloudType}),lat);
        upFracAll.(classTypes{cloudType})=cat(1,upFracAll.(classTypes{cloudType}),upFrac);
        upNumAll.(classTypes{cloudType})=cat(1,upNumAll.(classTypes{cloudType}),upNum);
        upMeanStrengthAll.(classTypes{cloudType})=cat(1,upMeanStrengthAll.(classTypes{cloudType}),upMeanStrength);
        downMeanStrengthAll.(classTypes{cloudType})=cat(1,downMeanStrengthAll.(classTypes{cloudType}),downMeanStrength);
        upMaxStrengthAll.(classTypes{cloudType})=cat(1,upMaxStrengthAll.(classTypes{cloudType}),upMaxStrength);
        downMaxStrengthAll.(classTypes{cloudType})=cat(1,downMaxStrengthAll.(classTypes{cloudType}),downMaxStrength);
        upRegsAll.(classTypes{cloudType})=cat(1,upRegsAll.(classTypes{cloudType}),upRegs);
        precShaftsAll.(classTypes{cloudType})=cat(1,precShaftsAll.(classTypes{cloudType}),precShafts);
    end
end

%% Save properties

disp('Saving output ...');

if ~strcmp(project,'spicule')
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','cloudLengthAll','cloudTopAll','cloudBaseAll','cloudLayersAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','divLevAll','meltDetAll','coldFracAll','sstAll','upFracAll','upNumAll','upRegsAll','precShaftsAll', ...
        'upMeanStrengthAll','downMeanStrengthAll','upMaxStrengthAll','downMaxStrengthAll','latAll','lonAll');
else
    save([figdir,project,'_cloudProps.mat'],'maxReflAll','meanReflAll','maxConvAll','meanConvAll', ...
        'cloudDepthAll','cloudLengthAll','cloudTopAll','cloudBaseAll','cloudLayersAll','maxTempAll','minTempAll','meanTempAll','maxPressAll','minPressAll','meanPressAll', ...
        'iceLevAll','divLevAll','meltDetAll','coldFracAll','upRegNumAll','upFracAll','upNumAll','upRegsAll','precShaftsAll', ...
        'upMeanStrengthAll','downMeanStrengthAll','upMaxStrengthAll','downMaxStrengthAll','latAll','lonAll');
end
