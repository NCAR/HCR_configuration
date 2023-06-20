% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
figdir=['/home/romatsch/plots/HCR/meltingLayer/hourly/',project,'/',freqData,'/dropSondeComp/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 7];

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/',freqData,'/'];
dropsondedir=['/run/media/romatsch/RSF0006/rsf/dropsondes/',project,'/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

% Dropsonde files
if strcmp(project,'otrec')
    dropFormat='nc';
else
    dropFormat='eol';
end

caseList = table2array(readtable(infile));

% Make N by 2 matrix of fieldname + value type
variable_names_types = [["time", "datetime"]; ...
			["sondeAlt", "double"]; ...
			["meltAltMeas", "double"]; ...
			["meltAltInt", "double"]; ...
			["meltAltEst", "double"]; ...
			["zeroDegAlt", "double"]];
% Make table using fieldnames & value types from above
compAlts = table('Size',[0,size(variable_names_types,1)],... 
	'VariableNames', variable_names_types(:,1),...
	'VariableTypes', variable_names_types(:,2));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=startTime;
    
    endTimeIn=datetime(caseList(aa,7:12));
    
    while endTime<endTimeIn
        endTime=endTime+hours(1);
        
        disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
        
        %% Loading dropsonde data
        [dropList,dropTimes]=dropSondeList(startTime,endTime,dropsondedir,dropFormat);
        
        if length(dropList)==0
            continue
        end
        
        [dropAlt,dropT]=getDropData(dropList,dropFormat);
        
        %% Load HCR data
        
        disp('Loading HCR data ...');
        data=[];
        
        if strcmp(freqData,'combined')
            data.HCR_DBZ=[];
        else
            data.DBZ=[];
        end
        data.MELTING_LAYER=[];
        data.ICING_LEVEL=[];
        %data.FLAG=[];
        data.TEMP=[];
        
        dataVars=fieldnames(data);
        
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
        
        dataVars=dataVars(~cellfun('isempty',dataVars));
        
        if strcmp(freqData,'combined')
            data.DBZ=data.HCR_DBZ;
            data=rmfield(data,'HCR_DBZ');
            data.LDR=data.HCR_LDR;
            data=rmfield(data,'HCR_LDR');
            data.VEL_CORR=data.HCR_VEL;
            data=rmfield(data,'HCR_VEL');
        end
        
        if isempty(data.DBZ)
            continue
        end
        
        %% Get indices
        
        elevenInds=find(data.MELTING_LAYER==11);
        twelveInds=find(data.MELTING_LAYER==12);
        thirteenInds=find(data.MELTING_LAYER==13);
        fourteenInds=find(data.MELTING_LAYER==14);
        
        twentyoneInds=find(data.MELTING_LAYER==21);
        twentytwoInds=find(data.MELTING_LAYER==22);
        twentythreeInds=find(data.MELTING_LAYER==23);
        twentyfourInds=find(data.MELTING_LAYER==24);
        
        %% Get melting layer alt and type
        compAltsHourIn=array2table(nan(length(dropAlt),size(variable_names_types,1)-1),...
            'VariableNames', variable_names_types(2:end,1));
        compAltsHour=cat(2,array2table(dropTimes,'VariableNames',{'time'}),compAltsHourIn);
        
        thisMeltAlt=nan;
        allSondeAlts={};
        for jj=1:length(dropAlt)
                        
            % Melting layer altitudes
            [minval sondeTimeInd]=min(abs(etime(datevec(data.time),datevec(dropTimes(jj)))));
            altCol=data.asl(:,sondeTimeInd);
            meltCol=data.MELTING_LAYER(:,sondeTimeInd);
            iceAlt=data.ICING_LEVEL(sondeTimeInd);
            % Melting layer alt and type
            meltIndSonde=find(meltCol==12 | meltCol==13 | meltCol==14);
            meltType=meltCol(meltIndSonde);
            if ~isempty(meltType)
                if meltCol(meltType==12)
                    compAltsHour.meltAltMeas(jj)=min(altCol(meltIndSonde));
                    thisMeltAlt=min(altCol(meltIndSonde));
                elseif meltCol(meltType==13)
                    compAltsHour.meltAltInt(jj)=min(altCol(meltIndSonde));
                    thisMeltAlt=min(altCol(meltIndSonde));
                else
                    compAltsHour.meltAltEst(jj)=min(altCol(meltIndSonde));
                    thisMeltAlt=min(altCol(meltIndSonde));
                end
            end
            % Zero degree alt
            zeroIndSonde=find(meltCol==11 | meltCol==21);
            if ~isempty(zeroIndSonde)
                compAltsHour.zeroDegAlt(jj)=min(altCol(zeroIndSonde));
            end
            
            % Dropsonde altitude
            tempAlt=cat(2,dropT{jj},dropAlt{jj});
            tempAlt(any(isnan(tempAlt),2),:) = [];
            signChT=diff(sign(tempAlt(:,1)));
            sondeAlts=tempAlt(signChT~=0,2);
            if ~isempty(sondeAlts) & ~isnan(thisMeltAlt)
                minDiffAlts=abs(sondeAlts-thisMeltAlt);
                compAltsHour.sondeAlt(jj)=sondeAlts(minDiffAlts==min(minDiffAlts));
            end
            allSondeAlts{end+1}=sondeAlts;
        end
        
        compAlts=cat(1,compAlts,compAltsHour);
        
        %% Plot
        
        timeMat=repmat(data.time,size(data.DBZ,1),1);
        dbzMasked=data.DBZ;
%        dbzMasked(data.FLAG>1)=nan;
        
        close all
        
        if etime(datevec(endTime),datevec(startTime))<=900
            newInds=1:1:length(data.time);
        elseif etime(datevec(endTime),datevec(startTime))<=3600
            newInds=1:10:length(data.time);
        else
            newInds=1:100:length(data.time);
        end
        
        % Resample for plotting
        newDBZ=dbzMasked(:,newInds);
        newASL=data.asl(:,newInds);
        newTEMP=data.TEMP(:,newInds);
        newFindMelt=data.MELTING_LAYER(:,newInds);
        newTime=data.time(newInds);
        
        fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,700]);
        
        subplot(2,1,1)
        hold on;
        sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
        view(2);
        sub1=colMapDBZ(sub1);
        scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
        scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
        
        scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
            'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
        scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
        scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
            'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
        scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
        scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
        scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
        
        % Dropsondes
        for jj=1:length(dropAlt)
            timeVec=repmat(dropTimes(jj),length(dropAlt{jj}),1);
            scatter(timeVec,dropAlt{jj}./1000,20,dropT{jj},'filled');
            set(gca,'clim',[-10 10])
            set(gca,'colormap',jet)
            if ~isempty(allSondeAlts{jj})
                scatter(repmat(dropTimes(jj),1,length(allSondeAlts{jj})),allSondeAlts{jj}/1000,20,'k','filled');
            end
        end
        
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title(['Flight ',num2str(aa),': Reflectivity, melting layer, and dropsonde temperature'])
        grid on
                
        set(colorbar,'visible','off')
        
        subplot(2,1,2)
        hold on;
        sub1=surf(newTime,newASL./1000,newTEMP,'edgecolor','none');
        view(2);
        scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
        scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
        
        scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
            'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
        scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
        scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
            'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
        scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
        scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
        scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
        
        % Dropsondes
        for jj=1:length(dropAlt)
            timeVec=repmat(dropTimes(jj),length(dropAlt{jj}),1);
            scatter(timeVec,dropAlt{jj}./1000,20,dropT{jj},'filled');
            set(gca,'clim',[-10 10])
            set(gca,'colormap',jet)
            if ~isempty(allSondeAlts{jj})
                scatter(repmat(dropTimes(jj),1,length(allSondeAlts{jj})),allSondeAlts{jj}/1000,20,'k','filled');
            end
        end
        
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title(['Melting layer, ERA5 and dropsonde temperature (C)'])
        grid on
                
        colorbar
        
        formatOut = 'yyyymmdd_HHMM';
        set(gcf,'PaperPositionMode','auto')
        print([figdir,'meltRefl_dropsondes_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
        
        startTime=endTime;
    end
end

writetable(compAlts,[figdir,project,'_meltLayer_dropsonde.txt'],'Delimiter',' ');

