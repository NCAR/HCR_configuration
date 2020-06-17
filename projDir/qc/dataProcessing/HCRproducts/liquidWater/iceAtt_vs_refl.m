% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made

project='otrec'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/liquidWaterHCR/iceAttVSrefl/'];

dataDir=HCRdir(project,quality,dataFreq);

% startTime=datetime(2019,8,7,17,5,0);
% endTime=datetime(2019,8,7,17,16,0);

ylimUpper=15;

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%% Load data

reflAttAll=[];

for mm=1:size(caseList,1)
    disp(['Flight ',num2str(mm)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(mm,1:6));
    endTime=datetime(caseList(mm,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.TEMP=[];
    data.TOPO=[];
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    data.freq=ncread(fileList{1},'frequency');
    
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
    
    surfMask=nan(size(data.time));
    
    %sort out non nadir pointing
    surfMask(data.elevation>-85)=0;
    
    %sort out land
    surfMask(data.TOPO>0)=0;
    
    % sort out data from below 2500m altitude
    surfMask(data.altitude<2500)=0;
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Calculate reflectivity sum inside and outside ocean surface to
    % distinguish clear air and cloud
    reflTemp=data.DBZ;
    
    % Remove bang
    reflTemp(data.FLAG==6)=nan;
    
    reflLin=10.^(reflTemp./10);
    reflOceanLin=nan(size(data.time));
    reflNoOceanLin=nan(size(data.time));
    
    for ii=1:length(data.time)
        if (~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-5)) & ~isnan(maxGate(ii))
            reflRay=reflLin(:,ii);
            reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
            reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
        end
    end
    
    % Remove data where reflectivity outside of ocean swath is more than
    % 0.8
    clearAir=find(reflNoOceanLin<=0.8);
    surfMask(clearAir)=2;
    surfMask(isnan(reflOceanLin))=0;
    
    % Find cloud data
    dbzMasked=data.DBZ;
    dbzMasked(data.FLAG>1)=nan;
    
    surfMask(find(any(~isnan(dbzMasked),1) & surfMask~=0 & surfMask~=2))=1;
    
    % Remove noise source cal, ant trans, and missing
    surfMask(find(any(data.FLAG>9,1) & surfMask~=0))=0;
    
    % Remove extinct
    surfMask(find(any(data.FLAG==3,1) & surfMask~=0))=0;
    
    %% Find melting layer and separate warm and cold precip
    meltInd=nan(size(data.time));
    warmRefl=nan(size(data.DBZ));
    coldRefl=nan(size(data.DBZ));
    
    for ii=1:length(meltInd)
        if surfMask(ii)~=0
            tempRay=data.TEMP(:,ii);
            if ~isempty(find(tempRay<=0))
                meltInd(ii)=max(find(tempRay<=0));
                
                warmRefl(meltInd(ii):end,ii)=dbzMasked(meltInd(ii):end,ii);
                coldRefl(1:meltInd(ii)-1,ii)=dbzMasked(1:meltInd(ii)-1,ii);
            end
        end
    end
    
    %% Calculate two way ice attenuation
    coldReflLin=10.^(coldRefl./10);
    coldReflLinSum=sum(coldReflLin,1,'omitnan');
    iceSpecAtt=0.0325.*coldReflLin;
    
    iceAttAll=iceSpecAtt.*(data.range(2)-data.range(1))./1000;
    iceAtt=sum(iceAttAll,1,'omitnan');
    
    %% Ice only inds
    surfMask(find(any(~isnan(warmRefl),1) & surfMask~=0))=0;
    
    if ~max(surfMask)==0
        %% Calculate clear air and cloudy ocean reflectivity
        
        reflSurf=data.DBZ(linInd);
        
        % Remove data with clouds
        dbzClear=nan(size(data.time));
        dbzCloud=nan(size(data.time));
        
        dbzClear(surfMask==2)=reflSurf(surfMask==2);
        dbzCloud(surfMask==1)=reflSurf(surfMask==1);
        
        % Add ice attenuation back in
        dbzCloudUsed=dbzCloud+iceAtt;
        
        % Clear air ocean refl
        clearShort=dbzClear;
        clearShort(isnan(clearShort))=[];
        
        meanClearShort=movmedian(clearShort,100,'omitnan');
        meanClear=nan(size(data.time));
        meanClear(~isnan(dbzClear))=meanClearShort;
        
        meanClear(1)=meanClear(min(find(~isnan(meanClear))));
        
        for jj=2:length(meanClear)
            if isnan(meanClear(jj))
                meanClear(jj)=meanClear(jj-1);
            end
        end
        
        %% Calculate liquid attenuation
        
        attLiq=nan(size(data.time));
        
        cloudInds=find(surfMask==1);
        
        for ii=1:length(cloudInds)
            dbzRay=dbzMasked(:,cloudInds(ii));
            cloudIndsRay=find(~isnan(dbzRay));
            
            if length(cloudIndsRay)>2
                
                attDiff=meanClear(cloudInds(ii))-dbzCloudUsed(cloudInds(ii));
                if attDiff>=0
                    attLiq(cloudInds(ii))=attDiff;
                else
                    attLiq(cloudInds(ii))=0;
                end
            end
        end
        
        coldReflLinSum(surfMask~=1)=nan;
        attLiq(surfMask~=1)=nan;
        iceAtt(surfMask~=1)=nan;
        
        %% Clean up
        reflAtt=cat(1,coldReflLinSum,attLiq);
        reflAtt=reflAtt';
        
        reflAtt(any(isnan(reflAtt),2),:) = [];
        
        reflAtt(find(reflAtt(:,2)==0),:)=[];
        reflAtt(find(reflAtt(:,1)<200),:)=[];
        
        reflAttAll=cat(1,reflAttAll,reflAtt);
        %% Scatter plot
        close all
        
        if size(reflAtt,1)>30
            f=fit(reflAtt(:,1),reflAtt(:,2),'poly2');
            
            f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12,'renderer','painters');
            
            plot(f,reflAtt(:,1),reflAtt(:,2));
            
            xlabel('Linear reflectivity sum (mm^6 m^{-3})')
            ylabel('2-way ice attenuation (dB)')
            
            ylimits=ylim;
            
            fstring = evalc('f');
            text(250,ylimits(2)-1,fstring)
            
            title([project,' RF',num2str(mm),', Ice Attenuation vs Reflectivity Sum'])
            
            print([figdir,project,'_iceAttVSrefl_RF',num2str(mm),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),...
                '_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
    end
end

%% Scatter plot
close all

f=fit(reflAttAll(:,1),reflAttAll(:,2),'poly2');

f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12,'renderer','painters');

plot(f,reflAttAll(:,1),reflAttAll(:,2));

xlabel('Linear reflectivity sum (mm^6 m^{-3})')
ylabel('2-way ice attenuation (dB)')

ylimits=ylim;

fstring = evalc('f');
text(250,ylimits(2)-1,fstring)

title([project,', Ice Attenuation vs Reflectivity Sum'])

print([figdir,project,'_iceAttVSrefl_RF',num2str(mm)],'-dpng','-r0')