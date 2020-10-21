% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz';
whichModel='era5';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%indir=HCRdir(project,quality,freqData);
indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socrates/combined/';

%[~,modeldir]=modelDir(project,whichModel,freqData);
modeldir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socratesMat/';

%% Run processing

% Go through flights
for ii=2:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    % Get model data
    disp('Loading 10hz data ...');
    
    model.meltLayer=[];
    
    model=read_model(model,modeldir,startTime,endTime);
    timeModelNum=datenum(model.time);
    
    %% Resample model data
    disp('Resampling model data ...');
    
    meltLong=cat(2,nan(size(model.meltLayer,1),20),model.meltLayer,nan(size(model.meltLayer,1),20));
    
    % Add extra columns for means etc.
    extraTime1=datetime(1899,1,1);
    extraTime=repmat(extraTime1,1,20);
    timeLong=cat(2,extraTime,model.time,extraTime);
    
    % Read first 2hz file to get time right
    timeLowResIn=ncread(fileList{1},'time');
    startTimeIn=ncread(fileList{1},'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeLowRes=startTimeFile+seconds(timeLowResIn);
    
    resolSecs=median(diff(timeLowRes));
    newTime=timeLowRes(1):resolSecs:model.time(end);
    
    newMelt=nan(size(model.meltLayer,1),length(newTime));
    
    % Create matrices with only zero flags and all others
    only0=meltLong;
    only0(only0~=0)=nan;
    no0=meltLong;
    no0(no0==0)=nan;
    no0(1:3,:)=nan;
    no0(end-3:end,:)=nan;
    
    % Find the right range pixels and write zero degree data
    
    disp('Searching range pixels and writing zero degree data ...');
    only0hor=movmean(only0,5,2,'omitnan');
    no0hor=movmean(no0,5,2,'omitnan');
    
    %maxLinesNo0=max(sum(~isnan(no0),1));
    rangePixNo0=zeros(1,length(newTime));
    
    timeIndColAll=nan(size(newTime));
    
    for kk=1:length(newTime)
        timeIndCol=find(timeLong==newTime(kk));
        if ~isempty(timeIndCol)
            timeIndColAll(kk)=timeIndCol;
            
            oneColOnly0=only0hor(:,timeIndCol);
            oneColOnly0(~isnan(oneColOnly0))=1;
            oneColOnly0(isnan(oneColOnly0))=0;
            
            connPixOnly0=bwconncomp(oneColOnly0);
            
            for ll=1:connPixOnly0.NumObjects
                rangePix=round(mean(connPixOnly0.PixelIdxList{ll}));
                newMelt(rangePix,kk)=0;
            end
            
            oneColNo0=no0hor(:,timeIndCol);
            oneColNo0(~isnan(oneColNo0))=1;
            oneColNo0(isnan(oneColNo0))=0;
            
            connPixNo0=bwconncomp(oneColNo0);
            
            for ll=1:connPixNo0.NumObjects
                if ll>5
                    stop1=1;
                end
                rangePixNo0(ll,kk)=round(mean(connPixNo0.PixelIdxList{ll}));
            end
        end
    end
    
    %% Find the right non zero data
    disp('Filling in melting layer data ...');
    
    for kk=1:length(newTime)
        pixCol=rangePixNo0(:,kk);
        pixCol(pixCol==0)=[];
        if ~isempty(pixCol) & ~isnan(timeIndColAll(kk))
            for ll=1:length(pixCol);
                subArea=no0(pixCol(ll)-2:pixCol(ll)+2,timeIndColAll(kk)-2:timeIndColAll(kk)+2);
                newMelt(pixCol(ll),kk)=min(min(subArea));
            end
        end
    end
        
    %% Prepare for plot
    
    disp('Plotting ...');
    
    zeroIndsH=find(model.meltLayer==0);
    oneIndsH=find(model.meltLayer==1);
    twoIndsH=find(model.meltLayer==2);
    threeIndsH=find(model.meltLayer==3);
    
    zeroIndsL=find(newMelt==0);
    oneIndsL=find(newMelt==1);
    twoIndsL=find(newMelt==2);
    threeIndsL=find(newMelt==3);
    
    timeMatHigh=repmat(model.time,size(model.meltLayer,1),1);
    timeMatLow=repmat(newTime,size(newMelt,1),1);
    
    pixVec=1:size(model.meltLayer,1);
    pixMatHigh=repmat(pixVec',1,size(model.meltLayer,2));
    pixMatLow=repmat(pixVec',1,size(newMelt,2));
    
    %% Plot
    
    close all
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1500,1000]);
    
    s1=subplot(2,1,1);
    hold on;
    scatter(timeMatHigh(zeroIndsH),pixMatHigh(zeroIndsH),10,'k','filled');
    scatter(timeMatHigh(oneIndsH),pixMatHigh(oneIndsH),10,'b','filled');
    scatter(timeMatHigh(twoIndsH),pixMatHigh(twoIndsH),10,'c','filled');
    scatter(timeMatHigh(threeIndsH),pixMatHigh(threeIndsH),10,'g','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylabel('Altitude (km)');
    xlim([newTime(1),newTime(end)]);
    ylim([0 770]);
    title({['Flight ',num2str(ii),', ',project,', ',...
        datestr(newTime(1),'HH:MM'),' to ',datestr(newTime(end),'HH:MM')];['10hz freezing level']});
    grid on
    
    s1=subplot(2,1,2);
    hold on;
    scatter(timeMatLow(zeroIndsL),pixMatLow(zeroIndsL),10,'k','filled');
    scatter(timeMatLow(oneIndsL),pixMatLow(oneIndsL),10,'b','filled');
    scatter(timeMatLow(twoIndsL),pixMatLow(twoIndsL),10,'c','filled');
    scatter(timeMatLow(threeIndsL),pixMatLow(threeIndsL),10,'g','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylabel('Altitude (km)');
    xlim([newTime(1),newTime(end)]);
    ylim([0 770]);
    title({['Flight ',num2str(ii),', ',project,', ',...
        datestr(newTime(1),'HH:MM'),' to ',datestr(newTime(end),'HH:MM')];['2hz freezing level']});
    grid on
        
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'meltRefl_Flight',num2str(aa)],'-dpng','-r0');
    
    %% Loop through HCR data files
    for jj=1:length(fileList)
        infile=fileList{jj};
        
        disp(infile);
        
        % Find times that are equal
        startTimeIn=ncread(infile,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(infile,'time')';
        timeHCR=startTimeFile+seconds(timeRead);
        
        timeHcrNum=datenum(timeHCR);
        
        [C,ia,ib] = intersect(timeHcrNum,timeModelNum);
        
        if length(timeHCR)~=length(ib)
            warning('Times do not match up. Skipping file.')
            continue
        end
        
        % Write output
        fillVal=-9999;
        
        modVars=fields(model);
        
        for kk=1:length(modVars)
            if ~strcmp((modVars{kk}),'time')
                modOut.(modVars{kk})=model.(modVars{kk})(:,ib);
                modOut.(modVars{kk})(isnan(modOut.(modVars{kk})))=fillVal;
                modOut.(modVars{kk})=modOut.(modVars{kk});
            end
        end
        
        % Open file
        ncid = netcdf.open(infile,'WRITE');
        netcdf.setFill(ncid,'FILL');
        
        % Get dimensions
        dimtime = netcdf.inqDimID(ncid,'time');
        dimrange = netcdf.inqDimID(ncid,'range');
        
        % Define variables
        netcdf.reDef(ncid);
        varidML = netcdf.defVar(ncid,'FREEZING_LEVEL','NC_SHORT',[dimrange dimtime]);
        netcdf.defVarFill(ncid,varidML,false,fillVal);
        netcdf.endDef(ncid);
        
        % Write variables
        netcdf.putVar(ncid,varidML,modOut.meltLayer);
        
        netcdf.close(ncid);
        
        % Write attributes
        ncwriteatt(infile,'FREEZING_LEVEL','long_name','freezing_level_and_zero_degree_level');
        ncwriteatt(infile,'FREEZING_LEVEL','standard_name','freezing_level_and_zero_degree_level');
        ncwriteatt(infile,'FREEZING_LEVEL','units','');
        ncwriteatt(infile,'FREEZING_LEVEL','flag_values',[0, 1, 2, 3]);
        ncwriteatt(infile,'FREEZING_LEVEL','flag_meanings','ERA5_zero_degree_level freezing_level_detected freezing_level_interpolated freezing_level_estimated');
        ncwriteatt(infile,'FREEZING_LEVEL','grid_mapping','grid_mapping');
        ncwriteatt(infile,'FREEZING_LEVEL','coordinates','time range');
        
    end
end