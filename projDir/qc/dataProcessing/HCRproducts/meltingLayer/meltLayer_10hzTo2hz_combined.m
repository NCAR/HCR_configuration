% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2, qc3
freqInData='10hz';
freqOutData='2hz';
qcVersion='v3.0';
whichModel='era5';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

outdir2hz=HCRdir(project,quality,qcVersion,freqOutData);

[~,indirMat10hz]=modelDir(project,whichModel,quality,qcVersion,freqInData);

figdir=[outdir2hz(1:end-4),'meltLayerPlots/2hzProcess/'];

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    clearvars -except caseList figdir formatOut freqData ii indirMat10hz infile ...
        outdir2hz project quality whichModel
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(outdir2hz,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    % Get model data
    disp('Loading 10hz data ...');
    
    model.meltLayer=[];
    model.iceLevel=[];
    
    model=read_model(model,indirMat10hz,startTime,endTime);
        
    %% Resample model data
    disp('Resampling model data ...');
    
    meltLong=cat(2,nan(size(model.meltLayer,1),20),model.meltLayer,nan(size(model.meltLayer,1),20));
    
    % Add extra columns for means etc.
    extraTime1=datetime(1899,1,1);
    extraTime=repmat(extraTime1,1,20);
    timeLong=cat(2,extraTime,model.time,extraTime);
    
    % Read first 2hz file to get time right
    timeLowRes=datetime(1899,1,1);
    aa=1;
    
    while timeLowRes(1)<startTime
        timeLowResIn=ncread(fileList{aa},'time');
        startTimeIn=ncread(fileList{aa},'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeLowRes=startTimeFile+seconds(timeLowResIn);
        aa=aa+1;
    end
    
    resolSecs=median(diff(timeLowRes));
    newTime=timeLowRes(1):resolSecs:model.time(end);
    
    newMelt=nan(size(model.meltLayer,1),length(newTime));
    
    % Create matrices with only zero flags and all others
    onlyLines=meltLong;
    onlyLines(onlyLines==10 | onlyLines==20)=nan;
    
    only0=onlyLines;
    only0(only0~=11 & only0~=21)=nan;
    no0=onlyLines;
    no0(no0==11 | no0==21)=nan;
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
                %newMelt(rangePix,kk)=0;
                newMelt(rangePix,kk)=min(only0hor(connPixOnly0.PixelIdxList{ll},timeIndCol));
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
    
    clear only0 no0 only0hor onlyLines no0hor
    
    sumMelt=sum(~isnan(newMelt),1);
    
    %% Remove data next to gaps
        
    sumMask=zeros(size(sumMelt));
    sumMask(sumMelt==0)=1;
    
    diffMask=diff(sumMask);
    
    checkMask=zeros(size(sumMelt));
    checkMask(diffMask==1)=1;
    checkMask(find(diffMask==-1)+1)=1;
    checkInds=find(checkMask==1);
    
    goThrough=-3:3;
    
    for bb=1:length(checkInds)
        for cc=1:length(goThrough)
            if checkInds(bb)+goThrough(cc)>0 & checkInds(bb)+goThrough(cc)<=length(sumMelt)
                newMelt(:,checkInds(bb)+goThrough(cc))=nan;
            end
        end
    end
    
    newMelt(:,find(sumMelt>10))=nan;
    
    %% Add in above and below data and icing level
    [C,ia1,ib1] = intersect(newTime,model.time);
    
    aBelow=model.meltLayer;
    aBelow(aBelow~=10 & aBelow~=20)=nan;
    
    newIce=nan(size(newTime));
    
    for jj=1:length(ib1)        
        newMeltCol=newMelt(:,ia1(jj));
        newMeltCol(isnan(newMeltCol))=aBelow(isnan(newMeltCol),ib1(jj));
        newMelt(:,ia1(jj))=newMeltCol;
        
        newIce(ia1(jj))=model.iceLevel(ib1(jj));
    end

    %% Prepare for plot
    
    disp('Plotting ...');
    
    elevenIndsH=find(model.meltLayer==11);
    twelveIndsH=find(model.meltLayer==12);
    thirteenIndsH=find(model.meltLayer==13);
    fourteenIndsH=find(model.meltLayer==14);
    
    twentyoneIndsH=find(model.meltLayer==21);
    twentytwoIndsH=find(model.meltLayer==22);
    twentythreeIndsH=find(model.meltLayer==23);
    twentyfourIndsH=find(model.meltLayer==24);
    
    elevenIndsL=find(newMelt==11);
    twelveIndsL=find(newMelt==12);
    thirteenIndsL=find(newMelt==13);
    fourteenIndsL=find(newMelt==14);
    
    twentyoneIndsL=find(newMelt==21);
    twentytwoIndsL=find(newMelt==22);
    twentythreeIndsL=find(newMelt==23);
    twentyfourIndsL=find(newMelt==24);
    
    timeMatHigh=repmat(model.time,size(model.meltLayer,1),1);
    timeMatLow=repmat(newTime,size(newMelt,1),1);
    
    pixVec=1:size(model.meltLayer,1);
    pixMatHigh=repmat(pixVec',1,size(model.meltLayer,2));
    pixMatLow=repmat(pixVec',1,size(newMelt,2));
    
    newIndsH=1:100:length(model.time);
    newMeltLayerH=model.meltLayer(:,newIndsH);
    newPixMatHigh=pixMatHigh(:,newIndsH);
    newTimeMatHigh=timeMatHigh(:,newIndsH);
    
    newIndsL=1:10:length(newTime);
    newMeltLayerL=newMelt(:,newIndsL);
    newPixMatLow=pixMatLow(:,newIndsL);
    newTimeMatLow=timeMatLow(:,newIndsL);
    
    %% Plot
    
    close all
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1500,900]);
    
    s1=subplot(3,1,1);
    hold on;
    sub1=surf(newTimeMatHigh,newPixMatHigh,newMeltLayerH,'edgecolor','none');
    s1.Colormap=([1 0 1;1 1 0]);
    view(2);
    scatter(timeMatHigh(twentyoneIndsH),pixMatHigh(twentyoneIndsH),10,'k','filled');
    scatter(timeMatHigh(elevenIndsH),pixMatHigh(elevenIndsH),10,...
        'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
    
    scatter(timeMatHigh(twentyfourIndsH),pixMatHigh(twentyfourIndsH),10,...
        'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
    scatter(timeMatHigh(twentythreeIndsH),pixMatHigh(twentythreeIndsH),10,...
        'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
    scatter(timeMatHigh(twentytwoIndsH),pixMatHigh(twentytwoIndsH),10,...
        'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
    
    scatter(timeMatHigh(fourteenIndsH),pixMatHigh(fourteenIndsH),10,'g','filled');
    scatter(timeMatHigh(thirteenIndsH),pixMatHigh(thirteenIndsH),10,'c','filled');
    scatter(timeMatHigh(twelveIndsH),pixMatHigh(twelveIndsH),10,'b','filled');
    
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim([0 size(timeMatHigh,1)]);
    ylabel('Altitude (km)');
    xlim([model.time(1),model.time(end)]);
    title({['Flight ',num2str(ii),', ',project,', ',...
        datestr(newTime(1),'HH:MM'),' to ',datestr(newTime(end),'HH:MM')];['10hz freezing level']});
    grid on
    
    
    s2=subplot(3,1,2);
    hold on;
    sub2=surf(newTimeMatLow,newPixMatLow,newMeltLayerL,'edgecolor','none');
    s2.Colormap=([1 0 1;1 1 0]);
    view(2);
    scatter(timeMatLow(twentyoneIndsL),pixMatLow(twentyoneIndsL),10,'k','filled');
    scatter(timeMatLow(elevenIndsL),pixMatLow(elevenIndsL),10,...
        'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
    
    scatter(timeMatLow(twentyfourIndsL),pixMatLow(twentyfourIndsL),10,...
        'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
    scatter(timeMatLow(twentythreeIndsL),pixMatLow(twentythreeIndsL),10,...
        'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
    scatter(timeMatLow(twentytwoIndsL),pixMatLow(twentytwoIndsL),10,...
        'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
    
    scatter(timeMatLow(fourteenIndsL),pixMatLow(fourteenIndsL),10,'g','filled');
    scatter(timeMatLow(thirteenIndsL),pixMatLow(thirteenIndsL),10,'c','filled');
    scatter(timeMatLow(twelveIndsL),pixMatLow(twelveIndsL),10,'b','filled');
    
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim([0 size(timeMatHigh,1)]);
    ylabel('Altitude (km)');
    xlim([model.time(1),model.time(end)]);
    title({['Flight ',num2str(ii),', ',project,', ',...
        datestr(newTime(1),'HH:MM'),' to ',datestr(newTime(end),'HH:MM')];['10hz freezing level']});
    grid on
    title('2hz freezing level');
    grid on
    
    s3=subplot(3,1,3);
    hold on;
    plot(model.time,model.iceLevel,'-b','linewidth',3);
    plot(newTime,newIce,'-r','linewidth',1.5);
    xlim([model.time(1),model.time(end)]);
    legend('10hz','2hz');
    title('Icing level')
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'melt_Flight',num2str(ii)],'-dpng','-r0');
    
    %% Loop through HCR data files
    timeModelNum=datenum(newTime);
    
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
        
        modOut.meltLayer=newMelt(:,ib);
        modOut.meltLayer(isnan(modOut.meltLayer))=fillVal;
        modOut.iceLevel=newIce(ib);
        modOut.iceLevel(isnan(modOut.iceLevel))=fillVal;
        
        % Open file
        ncid = netcdf.open(infile,'WRITE');
        netcdf.setFill(ncid,'FILL');
        
        % Get dimensions
        dimtime = netcdf.inqDimID(ncid,'time');
        dimrange = netcdf.inqDimID(ncid,'range');
        
        % Define variables
        netcdf.reDef(ncid);
        varidML = netcdf.defVar(ncid,'MELTING_LAYER','NC_SHORT',[dimrange dimtime]);
        netcdf.defVarFill(ncid,varidML,false,fillVal);
        varidIL = netcdf.defVar(ncid,'ICING_LEVEL','NC_FLOAT',[dimtime]);
        netcdf.defVarFill(ncid,varidIL,false,fillVal);
        netcdf.endDef(ncid);
        
        % Write variables
        netcdf.putVar(ncid,varidML,modOut.meltLayer);
        netcdf.putVar(ncid,varidIL,modOut.iceLevel);
        
        netcdf.close(ncid);
        
        % Write attributes
        ncwriteatt(infile,'MELTING_LAYER','long_name','melting_layer_and_zero_degree_level');
        ncwriteatt(infile,'MELTING_LAYER','standard_name','melting_layer_and_zero_degree_level');
        ncwriteatt(infile,'MELTING_LAYER','units','');
        ncwriteatt(infile,'MELTING_LAYER','flag_values',[10, 11, 12, 13, 14, 20, 21, 22, 23, 24]);
        ncwriteatt(infile,'MELTING_LAYER','flag_meanings',...
            'below_iceLev ERA5_zeroDeg_below_iceLev meltLayer_detected_below/at_iceLev meltLayer_interpolated_below/at_iceLev meltLayer_estimated_below/at_iceLev above_iceLev ERA5_zeroDeg_above_iceLev meltLayer_detected_above_iceLev meltLayer_interpolated_above_iceLev meltLayer_estimated_above_iceLev');
        ncwriteatt(infile,'MELTING_LAYER','is_discrete','true');
        ncwriteatt(infile,'MELTING_LAYER','grid_mapping','grid_mapping');
        ncwriteatt(infile,'MELTING_LAYER','coordinates','time range');
        
        ncwriteatt(infile,'ICING_LEVEL','long_name','icing_level');
        ncwriteatt(infile,'ICING_LEVEL','standard_name','icing_level');
        ncwriteatt(infile,'ICING_LEVEL','units','m');
        ncwriteatt(infile,'ICING_LEVEL','coordinates','time');
        
    end
end