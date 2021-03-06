% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/clouds/brightBand/',project,'/hourly/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 8];

indir=HCRdir(project,quality,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=startTime;
    
    endTimeIn=datetime(caseList(aa,7:12));
    
    while endTime<endTimeIn
        endTime=endTime+hours(1);
        
        disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
        data=[];
        
        %% Load data
        
        disp('Loading data ...');
        
        data.DBZ=[];
        data.LDR=[];
        data.VEL_CORR=[];
        data.TEMP=[];
        data.WIDTH=[];
        data.FLAG=[];
        data.TOPO=[];
        data.pitch=[];
        data.roll=[];
        
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
        
        if isempty(data.DBZ)
            continue
        end
        
        %% Find melting layer
        data.dbzMasked=data.DBZ;
        data.dbzMasked(data.FLAG>1)=nan;
        
        findMelt=f_meltLayer(data);
        zeroInds=find(findMelt==0);
        oneInds=find(findMelt==1);
        twoInds=find(findMelt==2);
        threeInds=find(findMelt==3);
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        %% Plot
        close all
        
        % Resample for plotting
        newInds=1:20:length(data.time);
        newDBZ=data.DBZ(:,newInds);
        newLDR=data.LDR(:,newInds);
        newVEL=data.VEL_CORR(:,newInds);
        newASL=data.asl(:,newInds);
        newTEMP=data.TEMP(:,newInds);
        newTime=data.time(newInds);
        
        fig3=figure('DefaultAxesFontSize',11,'position',[100,100,1500,1000]);
        
        s1=subplot(3,1,1);
        hold on;
        sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
        view(2);
        sub1=colMapDBZ(sub1);
        colIndsAll=1:length(data.time);
        scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'c','filled');
        scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'b','filled');
        scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title({['Flight ',num2str(aa),', ',project,', ',...
            datestr(data.time(1),'HH:MM'),' to ',datestr(data.time(end),'HH:MM')];['DBZ and freezing level']});
        grid on
        
        s2=subplot(3,1,2);
        hold on;
        sub2=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
        view(2);
        colorbar
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('LDR')
        grid on
        
        s3=subplot(3,1,3);
        hold on;
        sub2=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
        view(2);
        s3.Colormap=jet;
        colorbar
        caxis([-5 5]);
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('VEL')
        grid on
        
        formatOut = 'yyyymmdd_HHMM';
        set(gcf,'PaperPositionMode','auto')
        print([figdir,'meltRefl_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
        
        %% Calculate difference between zero degree and detected melting layer
        diffInd=[];
        
        for ii=1:size(findMelt,2)
            if data.elevation(ii)<-85
                bbRay=findMelt(:,ii);
                foundInd=min(find(bbRay==1));
                if ~isempty(foundInd)
                    zeroInd=min(find(bbRay==0));
                    diffInd=[diffInd foundInd-zeroInd];
                end
            end
        end
        
        meanDiff=mean(diffInd);
        meanMeters=meanDiff*(data.range(2)-data.range(1));
        
        disp(['Melting layer is on average ',num2str(meanDiff),' gates or ',num2str(meanMeters),' m below the zero degree level.'])
        
        startTime=endTime;
    end
end